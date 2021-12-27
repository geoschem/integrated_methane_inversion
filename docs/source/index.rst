Integrated Methane Inversion Workflow
=====================================
.. raw:: html

   <p>
   <a href="https://imi.readthedocs.io/en/latest/"><img src="https://img.shields.io/readthedocs/imi?label=ReadTheDocs"></a>
   <a href="https://github.com/ACMG-CH4/CH4_inversion_workflow/releases"><img src="https://img.shields.io/github/v/release/ACMG-CH4/CH4_inversion_workflow?include_prereleases&label=Latest%20Pre-Release"></a>
   </p>

.. important:: This is a prerelease of the Integrated Methane Inversion Workflow user guide.
   These pages are the most up-to-date and accurate instructions for the IMI Workflow, but they
   are still a work in progress. 
   
   Contributions (e.g., suggestions, edits, revisions) would be greatly appreciated. See
   :ref:`editing this guide <editing_this_user_guide>` and our :doc:`contributing guidelines <reference/CONTRIBUTING>`. 
   If you find something hard to understand---let us know!

The Integrated Methane Inversion (IMI) workflow is a cloud-based (or local cluster) tool for quantifying methane emissisions with 25 × 25 km\ :sup:`2`\  resolution by inversion of satellite observations from the TROPOspheric Monitoring Instrument (TROPOMI) and uses `GEOS-Chem <http://geos-chem.org>`_ as the forward model.

This site provides instructions for setting up and running the Integrated Methane Inversion workflow. Some instructions are specific to running the workflow on the AWS Cloud, but most steps in the workflow do not rely on using a cloud environment. We provide instructions for launching an EC2 instance from our custom AMI, configuring and running the inversion workflow, and analyzing the results with a ready-made jupyter notebook.


.. toctree::
   :maxdepth: 1
   :caption: Getting Started


   getting-started/quick-start.rst
   getting-started/imi-config-file.rst
   getting-started/ami-specifications.rst
   getting-started/manual-running.rst
   getting-started/running-with-tmux.rst
   getting-started/minimizing-cost-tips.rst

.. toctree::
   :maxdepth: 1
   :caption: Help & Reference

   reference/known-bugs.rst
   reference/SUPPORT.md
   reference/CONTRIBUTING.md
   geos-chem-shared-docs/editing_these_docs.rst
   reference/glossary.rst