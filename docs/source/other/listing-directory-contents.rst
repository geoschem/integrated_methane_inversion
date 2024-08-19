IMI directory contents 
======================

This page describes the contents of various file directories generated and populated by the IMI in the course of an inversion.

Inversion directory
-------------------

The inversion directory is where the IMI computes the Jacobian, obtains the optimal estimate of emissions, and saves the results.

It is located at ``/home/ubuntu/imi_output_dir/{YourRunName}/inversion``.

In addition to a shell script and several Python scripts used in the inversion, you will find
the following items in the inversion directory after completing an inversion:

.. list-table::
   :widths: 30, 70
   :class: tight-table
  
   * - ``data_converted/``
     - | Directory of Python ``.pkl`` files containing
       
         - TROPOMI observations
         - virtual TROPOMI observations of the GEOS-Chem reference simulation 
         - elements of the Jacobian matrix
         
       | for each TROPOMI orbit relevant to the inversion.
       | 
       | All quantities have been "converted" to 1D fields indexed by latitude and longitude.
   * - ``data_converted_posterior/``
     - | Directory of Python ``.pkl`` files containing
       
         - TROPOMI observations
         - virtual TROPOMI observations of the GEOS-Chem posterior simulation
         
       | for each TROPOMI orbit relevant to the inversion.
       |
       | All quantities have been "converted" to 1D fields indexed by latitude and longitude.
   * - ``data_geoschem/``
     - | Directory of ``.nc`` files containing daily GEOS-Chem ``SpeciesConc`` output from the reference simulation. 
       |
       | These files are used to generate virtual TROPOMI observations for comparison with the true observations.
   * - ``data_geoschem_posterior/``
     - | Directory of ``.nc`` files containing daily GEOS-Chem ``SpeciesConc`` output from the posterior simulation. 
       |
       | These files are used to generate virtual TROPOMI observations for comparison with the true observations.
   * - ``data_sensitivities/``
     - | Directory of ``.nc`` files containing daily 4-D GEOS-Chem sensitivities to perturbations in the 
         state variables of the inversion (i.e., in the emission elements being optimized). 
       |
       | The data have dimensions ``(element, lev, lat, lon)``, where ``element`` is the emission element id
         (state variable id) and ``lev`` is the vertical dimension. 
       |
       | These files are used to compute the Jacobian matrix by application of the TROPOMI operator.
   * - ``inversion_result.nc``
     - | File containing the raw output of the inversion (``invert.py``) as vectors (posterior emission
         estimate) and matrices (posterior error covariance matrix, averaging kernel matrix).
   * - ``gridded_posterior.nc``
     - | File containing the posterior emission estimate, posterior error covariance matrix, and averaging
         kernel matrix projected onto the 2-D inversion grid.
   * - ``visualization_notebook.ipynb``
     - | Jupyter notebook for quickly visualizing key results of the inversion.
   * - ``visualization_notebook.html``
     - | For convenience, we run the Jupyter notebook as part of the IMI and output the resultant notebook 
         to an html file, which can be viewed and rendered by your browser.
