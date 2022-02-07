Integrated Methane Inversion (IMI)
==================================
.. raw:: html

   <p>
   <a href="https://imi.readthedocs.io/en/latest/"><img src="https://img.shields.io/readthedocs/imi?label=ReadTheDocs"></a>
   <a href="https://github.com/ACMG-CH4/CH4_inversion_workflow/releases"><img src="https://img.shields.io/github/v/release/ACMG-CH4/CH4_inversion_workflow?include_prereleases&label=Latest%20Pre-Release"></a>
   </p>

.. important:: This is a prerelease of the Integrated Methane Inversion (IMI) user guide.
   These pages contain the most up-to-date and accurate instructions for the IMI, but they
   are still a work in progress. 
   
   Contributions (e.g., suggestions, edits, revisions) would be greatly appreciated. See
   :ref:`editing this guide <editing_this_user_guide>` and our :doc:`contributing guidelines <reference/CONTRIBUTING>`. 
   If you find something hard to understand, let us know!

The Integrated Methane Inversion (IMI) workflow is a cloud-computing tool for quantifying methane emissions
by inversion of satellite observations from the TROPOspheric Monitoring Instrument (TROPOMI). 
It uses `GEOS-Chem <http://geos-chem.org>`_ as forward model for the inversion and infers methane emissions 
at 25 × 25 km\ :sup:`2`\  resolution.

This site provides instructions for using the IMI, including launching an AWS compute instance with pre-installed IMI 
software, configuring and running an inversion, and analyzing the results with a ready-made jupyter notebook.

Some instructions are specific to the Amazon Web Services (AWS) cloud, but the IMI can also be run on a local compute cluster. 


.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   getting-started/quick-start.rst
   getting-started/imi-config-file.rst
   getting-started/running-with-tmux.rst
   getting-started/ami-specifications.rst
   getting-started/minimizing-cost-tips.rst

.. toctree::
   :maxdepth: 1
   :caption: Advanced Features

   advanced/custom-state-vector.rst
   advanced/custom-prior-emissions-hemco.rst
   advanced/inversion-ensemble.rst

.. toctree::
   :maxdepth: 1
   :caption: Help & Reference

   reference/known-bugs.rst
   reference/SUPPORT.md
   reference/CONTRIBUTING.md
   geos-chem-shared-docs/editing_these_docs.rst
   reference/glossary.rst
