.. _imi-known-bugs:

#####################
Known bugs and issues
#####################

Please see our `Issue tracker on GitHub
<https://github.com/geoschem/integrated_methane_inversion/issues>`_ for a list of recent
bugs and fixes.

===================
Current bug reports
===================

These `bug reports (on GitHub)
<https://github.com/geoschem/integrated_methane_inversion/issues?q=state%3Aopen%20label%3Abug>`_
are currently unresolved. We hope to fix these in future releases.

=======================================
Other issues that you should know about
=======================================

.. _gc-known-bugs-gcc12:

GCC 12.2.0 is discontinued in Spack v1.0.0
------------------------------------------

As of Spack v1.0, `spack-packages <https://packages.spack.io/>`_ has
been split off into its own separate repository. This change includes
the unfortunate deprecation of the :program:`GNU Compiler Collection
(GCC)` version 12.2.0. It appears that only the most recent minor
release in each major release is now treated as stable. These
deprecations are updated promptly for example, GCC 12.4.0 is already
marked as deprecated just 10 days after the release of GCC 12.5.0.

Deprecated GCC versions are no longer listed with the :command:`spack
info` command, so rather than warning users about deprecation, Spack
simply fails with an unhelpful error message about not being able to
satisfy the request.

For the time being, we recommend that you use `Spack release v0.23.1
<https://github.com/spack/spack/releases/tag/v0.23.1>`_ which still
supports GCC 12.2.0 and related libraries.  Please see our
:ref:`spackguide` Supplemental Guide for an updated Spack
installation workflow.

Discontinuity in GEOS-FP convection at 01 Jun 2020
--------------------------------------------------

The convection scheme used to generate archived GEOS-FP meteorology
files changed from RAS to Grell-Freitas starting 01 June 2020 with
impact on vertical transport. Discussion and analysis of the impact is
available at https://github.com/geoschem/geos-chem/issues/1409.

In addition, there is a bug in convective precipitation flux following
the switch where all values are zero. While this bug is automatically
fixed by calling different convection schemes in GEOS-Chem, the
convection scheme called is based only on run start date. This means
that using meteorology for a year different than simulation year may
result in choosing the wrong convection scheme. It also means that
simulations which span 01 June 2020 will incorrectly use the same
convection scheme for the entire run.

Due to these issues we recommend splitting up GEOS-FP runs in time
such that a single simulation does not run across 01 June 2020.
Instead. set one run to stop on 01 June 2020 and then restart a new
run from there. If you wish to use a GEOS-FP meteorology year
different from your simulation year please create a GEOS-Chem GitHub
issue for assistance.

============================
Bugs that have been resolved
============================

These `bugs (reported on GitHub) <https://github.com/geoschem/integrated_methane_inversion/pulls?q=label%3Abugfix+>`_ have been resolved.
