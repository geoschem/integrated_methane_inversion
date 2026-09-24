"""Normal (Gaussian) inversion with positivity enforced by a softplus constraint.

The ROI emission elements are reparametrized x = s*log(1 + e^{z/s}) so the posterior stays
positive (no log-barrier, no bounded QP); BC and OH elements stay linear. Levenberg-Marquardt on
the pre-assembled normal-equation products KTinvSoK = K^T So^-1 K, KTinvSoy = K^T So^-1 (y-F(x_A))
and ytinvSoy = (y-F)^T So^-1 (y-F), so it inherits whatever So invert.py built. (The lognormal
transform x = exp(z) is the separate LognormalErrors path.)

With posterior_mean the state is taken at the posterior mean E[softplus(z)] (Gauss-Hermite), not
the mode, matching the lognormal convention (#464). DOFS and the averaging kernel are the linear
ones, not distorted by the transform.
"""

import numpy as np

# Gauss-Hermite nodes/weights for the posterior-mean quadrature.
_GH_NODES, _GH_WEIGHTS = np.polynomial.hermite.hermgauss(24)


def softplus(z, scale):
    """x = scale * log(1 + exp(z/scale)); numerically stable via logaddexp."""
    return scale * np.logaddexp(0.0, z / scale)


def softplus_derivative(z, scale):
    """dx/dz = 1 / (1 + exp(-z/scale)) = logistic(z/scale)."""
    return 1.0 / (1.0 + np.exp(-np.clip(z / scale, -500.0, 500.0)))


def softplus_inverse(x, scale):
    """z such that softplus(z, scale) = x, for x > 0."""
    return scale * np.log(np.expm1(np.clip(x / scale, 1e-12, 700.0)))


def softplus_posterior_mean(z, variance, scale):
    """E[softplus(Z)] for Z ~ N(z, variance), by Gauss-Hermite quadrature."""
    sigma = np.sqrt(np.clip(variance, 0.0, None))
    acc = np.zeros_like(np.asarray(z, dtype=float))
    for node, weight in zip(_GH_NODES, _GH_WEIGHTS):
        acc += weight * softplus(z + np.sqrt(2.0) * sigma * node, scale)
    return acc / np.sqrt(np.pi)


def softplus_posterior_mean_derivative(z, variance, scale):
    """d/dz E[softplus(Z)] for Z ~ N(z, variance)."""
    sigma = np.sqrt(np.clip(variance, 0.0, None))
    acc = np.zeros_like(np.asarray(z, dtype=float))
    for node, weight in zip(_GH_NODES, _GH_WEIGHTS):
        acc += weight * softplus_derivative(z + np.sqrt(2.0) * sigma * node, scale)
    return acc / np.sqrt(np.pi)


def inversion_diagnostics(delta, KTinvSoK, KTinvSoy, ytinvSoy, inv_Sa, gamma, n_obs, n_roi):
    """DOFS, J_A/DOFS and J_O/(m-DOFS) from the linear averaging kernel
    A = (gamma KTinvSoK + Sa^-1)^-1 (gamma KTinvSoK). delta = xhat - x_A;
    J_A = delta_roi^T Sa^-1 delta_roi; J_O = ytinvSoy - 2 delta.KTinvSoy + delta.KTinvSoK.delta."""
    Md = gamma * np.asarray(KTinvSoK, dtype=float)
    A_linear = np.linalg.solve(Md + inv_Sa, Md)
    dofs = float(np.clip(np.diag(A_linear)[:n_roi], 0.0, 1.0).sum())

    d = np.asarray(delta, dtype=float)
    J_A = float(d[:n_roi] @ (inv_Sa[:n_roi, :n_roi] @ d[:n_roi]))
    J_O = float(ytinvSoy - 2.0 * (d @ KTinvSoy) + d @ (KTinvSoK @ d))
    return {
        "DOFS": dofs,
        "J_A": J_A,
        "J_A_over_DOFS": J_A / max(dofs, 1e-9),
        "J_O": J_O,
        "J_O_over_m_minus_DOFS": J_O / max(n_obs - dofs, 1e-9),
    }


def run_softplus(
    KTinvSoK,
    KTinvSoy,
    ytinvSoy,
    inv_Sa,
    inv_Sa_constraint,
    n_obs,
    n_roi,
    x_prior,
    scale=0.1,
    gamma=1.0,
    kappa=10.0,
    posterior_mean=True,
    mean_relaxation=0.5,
    mean_every=20,
    tol=5e-3,
):
    """Softplus-constrained normal inversion in normal-equation space.

    inv_Sa is the linear prior inverse (posterior covariance + diagnostics); inv_Sa_constraint the
    (OH-weighted) prior inverse used in the LM step. The first n_roi elements use the transform, the
    rest (BC, OH) stay linear. x_prior is the full prior (1 on ROI/OH, 0 on BC). scale = softplus
    smoothing s; gamma, kappa = regularization and LM damping. Iterates to convergence (like the
    lognormal solver, no iteration cap).

    Returns xhat, delta, S_post, A, diagnostics, n_iter.
    """
    n = KTinvSoK.shape[0]
    Md = gamma * np.asarray(KTinvSoK, dtype=float)
    vd = gamma * np.asarray(KTinvSoy, dtype=float)

    # Prior in transform space: invert softplus over the ROI, identity elsewhere.
    x_prior = np.asarray(x_prior, dtype=float)
    z_prior = x_prior.copy()
    z_prior[:n_roi] = softplus_inverse(x_prior[:n_roi], scale)

    z = z_prior.copy()
    variance = np.zeros(n)

    def field(state, var):
        scale_factors = state.copy()
        deriv = np.ones(n)
        if posterior_mean:
            scale_factors[:n_roi] = softplus_posterior_mean(state[:n_roi], var[:n_roi], scale)
            deriv[:n_roi] = softplus_posterior_mean_derivative(state[:n_roi], var[:n_roi], scale)
        else:
            scale_factors[:n_roi] = softplus(state[:n_roi], scale)
            deriv[:n_roi] = softplus_derivative(state[:n_roi], scale)
        return scale_factors, deriv

    # No iteration cap (like the lognormal solver): iterate until the ROI scale factors and the
    # posterior-mean variance both stop changing (step, var_change < tol).
    n_iter = 0
    it = 0
    while True:
        n_iter = it + 1
        sf, deriv = field(z, variance)
        data_hessian = Md * np.outer(deriv, deriv)
        gradient = deriv * (vd - Md @ (sf - x_prior))
        z_new = z + np.linalg.solve(
            data_hessian + (1.0 + kappa) * inv_Sa_constraint,
            gradient - inv_Sa_constraint @ (z - z_prior),
        )
        sf_new, _ = field(z_new, variance)
        step = np.max(np.abs(sf_new[:n_roi] - sf[:n_roi]) / np.maximum(sf[:n_roi], 1e-9))
        z = z_new
        # Refresh the posterior-mean variance only periodically -- updating it every iteration
        # chases a moving mean and won't converge on tight priors.
        var_change = 0.0
        if posterior_mean and (it % mean_every == 0 or step < tol):
            var_target = np.clip(np.diag(np.linalg.inv(data_hessian + inv_Sa))[:n_roi], 0.0, None)
            var_prev = variance[:n_roi].copy()
            variance[:n_roi] = (1.0 - mean_relaxation) * variance[:n_roi] + mean_relaxation * var_target
            var_change = float(np.max(np.abs(variance[:n_roi] - var_prev)))
        if it > 0 and step < tol and var_change < tol:
            break
        it += 1

    sf, deriv = field(z, variance)
    data_hessian = Md * np.outer(deriv, deriv)
    S_post = np.linalg.inv(data_hessian + inv_Sa)
    # Linear (normal) averaging kernel -- the data resolution, not distorted by the transform.
    A = np.linalg.solve(Md + inv_Sa, Md)

    xhat = sf.copy()
    delta = xhat - x_prior
    diagnostics = inversion_diagnostics(delta, KTinvSoK, KTinvSoy, ytinvSoy, inv_Sa, gamma, n_obs, n_roi)
    return xhat, delta, S_post, A, diagnostics, n_iter

