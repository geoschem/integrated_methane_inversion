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
