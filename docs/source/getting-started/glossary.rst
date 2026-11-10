.. _imi_glossary:

Terminology
==================================

We define here common technical terms used in the IMI documentation.  For an introductory textbook description of inverse 
methods using the IMI terminology, see Chapter 11: Inverse Methods for Atmospheric Chemistry, in  Modeling of Atmospheric 
Chemistry by G.P. Brasseur and D.J. Jacob, Cambridge University Press, 2017. 


.. glossary::
   :sorted:

   Averaging kernel matrix (A)
      A square matrix of dimension (n x n) describing the sensitivity of the **posterior estimate** of the **state vector**
      (:math:`\hat{x}`, dimension n x 1) to its true value x: **A** = :math:`\frac{\delta {\hat{x}}}{\delta x}`. The diagonal terms of 
      the averaging kernel matrix are the **averaging kernel sensitivities** and indicate the ability of the observations to 
      quantify individual state vector elements independently from their **prior estimates** and other state vector elements 
      (1 = fully, 0 = not at all). The off-diagonal terms of the averaging kernel matrix indicate the degree 
      to which inversion results for pairs of state vector elements are independent (0 = fully, ±1 = not at all) .  

   Averaging kernel sensitivity
      Estimates how sensitive the posterior solution for a given state vector (emission) element is to 
      observations as opposed to the prior estimate. An emission element with averaging kernel sensitivity 0 
      is not quantified by the observations at all, and the inversion returns the prior value for that element. 
      An emission element with averaging kernel sensitivity 1 is fully quantified by the observations, and the 
      inversion results for that element are independent of the prior estimate. 

   Boundary Conditions
      3D time-depenent methane concentration fields applied at the edges of the **inversion domain**. They are from a
      global **GEOS-Chem** simulation for the **period of interest** locally corrected with smoothed observations from
      the satellite instrument used in the inversion. See Estrada et al. (2025) for details. The mean boundary
      conditions along the four cardinal edges can be included in the **state vector** to be optimized by the inversion.

   Buffer elements
      The buffer elements in a regional inversion are elements of the **state vector** outside of the **region 
      of interest** but within the **inversion domain**. They are optimized in the IMI to correct for the influence 
      of emissions just outside the region of interest, and as an alternative or supplement to the correction 
      of **boundary conditions**. The user selects the number of buffer elements in the configuration file 
      (default is 8) and the IMI then applies **smart clustering** to aggregate **grid cells** for each buffer element.

   Cost function (J)
      The quantity minimized in the inversion to obtain the **posterior estimate**. It weighs the cost 
      of departing from the **prior estimate** (cost function term :math:`J_A`) and from the observations (cost 
      function term :math:`J_O`) 
      following Bayes’ theorem: :math:`J = J_A + J_O`.

   Degrees of freedom for signal (DOFS)
      This is the trace of the **averaging kernel matrix** (the sum of the **averaging kernel sensitivities**) for
      the region of interest. The DOFS measures the information content of the observations towards optimizing the
      **state vector**, that is  the number of independent pieces of information on the state vector that the 
      observations can quantify independently of the **prior estimate**.

   Dry column mole fraction (X)
      Satellite observations of atmospheric methane are reported as the dry column mole fraction (:math:`X_{CH4}`),
      which is the ratio between the column concentrations [molecules cm -2] of methane (ΩCH4) and dry air (Ωa): 
      XCH4 = ΩCH4/ Ωa. Unit for :math:`X_{CH4}` is parts per billion (ppb).

   Error covariance matrix (S)
      Symmetric square n x n matrix describing the error variance of a vector n x 1.  The diagonal terms are the 
      error variances of the individual vector elements and the off-diagonal terms are the error covariances 
      between vector elements. The **prior error covariance matrix** (:math:`S_A`) describes the error on the **prior 
      estimate** of the **state vector**. The **observational error covariance matrix** (:math:`S_O`) describes the error on the 
      **observation vector** including contributions from the **satellite retrieval** and from the **forward model**. 
      The **posterior error covariance matrix** (:math:`\hat{S}`) describes the error on the **posterior estimate** 
      of the state vector.

   Forward model (F)
      The forward model in the inversion is the chemical transport model (CTM) used to relate emissions to 
      the observed **column concentrations** through atmospheric transport: y = F(x) where y is the **observation 
      vector** and x is the state vector. The IMI uses the **GEOS-Chem CTM** as forward model.

   GEOS-Chem chemical transport model (CTM)
      GEOS-Chem is the **forward model** used in the IMI to relate emissions to observed concentrations 
      through the **observation operator**. GEOS-Chem simulates atmospheric transport with the NASA GEOS 
      assimilated meteorological data on a user-selected rectilinear grid of either 0.125°×0.15625° 
      (≈12x12 km2) or 0.25°x0.3125° (≈25x25 km2).  Coarser resolution is available for global inversions.

   Grid cell
      A rectilinear (latitude x longitude) spatial element of the GEOS-Chem forward model,  
      either 0.125°×0.15625° (≈12x12 km2) or 0.25°x0.3125° (≈25x25 km2) as selected by the user. 
      The IMI regrids the prior emission estimates to the GEOS-Chem grid and returns posterior emission 
      estimates on that grid. The emission state vector elements are also defined on that grid and 
      may be clustered to aggregated grid cells using smart clustering. 

   Inversion domain
      The **region of interest** and a surrounding **buffer region**. The inversion domain is always 
      rectilinear in shape and is defined in the IMI configuration file by latitude and longitude edges.
   
   Inversion ensemble
      The IMI configuration file allows users to automatically carry out a collection of inversion simulations, 
      called an ensemble. This ensemble varies inversion parameters including error variances on 
      **prior emissions**, error variances on **boundary conditions**, and **regularization parameter** γ. 
      As part of the ensemble, users can also select different error probability density functions (pdfs)
      for the prior, including normal, lognormal, and softplus. They can use different subsets of observations 
      such as with versus without glint observations. The IMI discards ensemble members that cause 
      overfit or underfit, based on the prior term of the **cost function** (:math:`J_A`). The range of remaining 
      ensemble members gives a confidence range for the **posterior emission estimate** that is more 
      conservative than the **posterior error covariance matrix**. Large ensembles are inexpensive to 
      compute in the IMI by virtue of the analytical inversion approach.

   Jacobian matrix (K)
      The sensitivity of observed concentrations to emissions through atmospheric transport 
      as simulated by the forward model:  :math:`K = \frac{\delta y}{\delta x}`
      where **y** is the **observation vector** and **x** is the **state vector**. The Jacobian matrix is the 
      derivative of the forward model relative to x:  :math:`K = \frac{\delta y}{\delta x} = \Delta_xF`.
      If the forward model is linear in **x**, then the Jacobian matrix does not depend on **x** and fully describes 
      the forward model: **y = Kx + c** where **c** is a vector of initial conditions. The Jacobian matrix 
      is explicitly constructed in the IMI to allow analytical inversion.

   Kalman filter
      An inversion procedure where the full **period of interest** (for example a year) is partitioned into 
      sub-periods (for example a month) with optimization of emissions for each subperiod. The **posterior 
      estimate** for a subperiod is used as **prior estimate** for the next subperiod, with an optional nudging
      factor back to the original prior estimate. The Kalman filter can also be set up in an open-ended 
      near-real-time mode where emissions are automatically updated at the end of each subperiod.

   Observation operator
      Operator applied to vertical concentration profiles in the forward model to reproduce :math:`X_{CH4}` as would 
      be observed by the satellite instrument. It generally involves regridding of the vertical profile 
      and applying the averaging kernels and prior vertical profiles from the **satellite retrieval**.

   Observation vector (y)
      Ensemble of m observations used in the inversion and arranged as a (1 x m) vector. 
      The IMI uses **super-observations** to construct the observation vector.

   Period of interest
      The period for which the IMI will optimize mean emissions based on TROPOMI observations. 
      It is specified in the configuration file by start and end dates. Only observations made during 
      this period are considered. The period of interest can be split into shorter sub-periods 
      (to be optimized sequentially) through the IMI Kalman filter feature.

   Posterior emission estimate
      Best estimate of emissions after performing the inversion.

   Prior emission estimate
      Best estimate of emissions before performing the inversion, used to guide the optimization
      by the satellite observations. The prior estimate is usually from 2-D gridded bottom-up 
      inventories. The IMI has a library of default bottom-up inventories as well as optional 
      inventories that the user can select through the configuration file.  The user can also 
      provide their own prior emission estimates for the **region of interest**.

   Region of interest
      The region for which optimization of emissions is desired.  This can be specified by the 
      user as a rectilinear domain or as a shapefile. The IMI optimizes emissions within the 
      region of interest at the user-requested grid resolution of either 0.125°×0.15625° (≈12x12 km2) 
      or 0.25°x0.3125° (≈25x25 km2). **Smart clustering** is applied if desired to reduce the size of the
      **state vector**. 

   Regularization parameter (γ)
      A factor γ ≤ 1applied to the observational term of the inversion **cost function** to decrease
      the weight of observations and avoid overfit. Although γ = 1 in straight Bayesian inference, 
      overfit may result if the observation error (including covariances) has been underestimated.  
      If an **inversion ensemble** is requested in the IMI configuration file, different values of γ are 
      tested and the values producing overfit or underfit are rejected.
   
   Satellite retrieval
      The retrieval is the operation by which the backscattered **SWIR** solar spectra measured by the 
      satellite instrument are converted to a **dry column mole fraction** (:math:`X_{CH4}`). It uses a 
      radiative transfer model operating on a vertical grid an where the column concentration 
      is retrieved with a certain vertical sensitivity (averaging kernel vector) and prior vertical 
      profile, both of which are generally reported together with :math:`X_{CH4}` and are necessary for applying 
      the **observation operator** of the **forward model**.

   Shortwave infrared (SWIR)
      The solar spectrum radiation range used for satellite observvation of methane. 

   Smart clustering: 
      Reducing the number of **state vector** elements by grouping together grid cells within 
      the **region of interest**, The user specifies the desired number of state vector elements 
      in the configuration file and the IMI then clusters weakly emitting and/or poorly observed 
      grid cells in order to preserve high resolution in areas of high emissions and dense observations.

   State vector
      The collection of elements (state vector elements) to be optimized in the inversion. 
      This includes 2D gridded emission elements in the **inversion domain** and may also include 
      **boundary conditions** and the OH sink as specified by the user in the configuration file.
   
   Super-observations
      Observations averaged over the **GEOS-Chem** grid and hourly time step, and used to construct 
      the **observation vector** for the inversion.

   
   