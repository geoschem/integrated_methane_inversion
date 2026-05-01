Integrated Methane Inversion (IMI)
==================================
.. raw:: html

   <p>
     <a href="https://github.com/geoschem/integrated_methane_inversion/releases"><img src="https://img.shields.io/github/v/release/geoschem/integrated_methane_inversion?label=Latest%20Stable%20Release" alt="Latest release"></a>
     <a href="https://github.com/geoschem/integrated_methane_inversion/releases/"><img src="https://img.shields.io/github/release-date/geoschem/integrated_methane_inversion" alt="Release date"></a>
     <a href="https://doi.org/10.5281/zenodo.6081933"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.6081933.svg" alt="DOI"></a><br/>
     <a href="https://github.com/geoschem/integrated_methane_inversion/blob/main/LICENSE.md"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License"></a>
   <a href="https://imi.readthedocs.io/en/latest/"><img src="https://img.shields.io/readthedocs/imi?label=ReadTheDocs"></a>
   </p>

The Integrated Methane Inversion (IMI) workflow is a cloud-computing tool for quantifying methane emissions
by inversion of satellite observations from the TROPOspheric Monitoring Instrument (TROPOMI). 
It uses `GEOS-Chem <http://geos-chem.org>`_ as forward model for the inversion and infers methane emissions 
at up to 12 x 12 km\ :sup:`2`\  resolution.

This site provides instructions for using the IMI, including launching an AWS compute instance, 
configuring and running an inversion, and analyzing the results with a ready-made jupyter notebook.

Some instructions are specific to the Amazon Web Services (AWS) cloud, but the IMI can also be run on a 
local compute cluster by either manually building the environment or using a docker container. 


.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   getting-started/quick-start.rst
   getting-started/best-practices.rst
   getting-started/glossary.rst
   getting-started/faqs.rst

.. toctree::
   :maxdepth: 1
   :caption: User Guide

   user-guide/imi-config-file.rst
   user-guide/imi-preview.rst
   user-guide/imi-input-data.rst
   user-guide/imi-output.rst
   user-guide/ami-specifications.rst
   user-guide/setting-up-jupyter.rst

.. toctree::
   :maxdepth: 1
   :caption: Advanced Features
   
   advanced/custom-state-vector.rst
   advanced/using-clustering-options.rst
   advanced/custom-prior-emissions-hemco.rst
   advanced/custom-region.rst
   advanced/inversion-ensemble.rst
   advanced/kalman-filter-mode.rst
   advanced/running-with-tmux.rst
   advanced/local-cluster.rst
   advanced/imi-docker-container.rst

.. toctree::
   :maxdepth: 1
   :caption: Other

   other/common-configurations.rst
   other/minimizing-cost-tips.rst

.. toctree::
   :maxdepth: 1
   :caption: Help & Reference

   help-and-reference/version-history.rst	
   help-and-reference/known-bugs.rst
   help-and-reference/SUPPORT.md
   help-and-reference/CONTRIBUTING.md
   geos-chem-shared-docs/editing_these_docs.rst
