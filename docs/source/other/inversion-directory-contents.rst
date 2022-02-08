Contents of the inversion directory
===================================

The inversion directory is where the IMI computes the Jacobian, obtains the optimal estimate of emissions, and saves the results.
This page describes the contents of the directory.

After completing an inversion, navigate to the inversion directory and display its contents::

    $ cd /home/ubuntu/CH4_Workflow/{YourRunName}/inversion
    $ ll

In addition to a shell script and several Python scripts used in the inversion, you will see the following items:

.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``data_converted``
     - temp
   * - ``data_converted_posterior``
     - temp
   * - ``data_GC``
     - temp
   * - ``Sensi``
     - temp
   * - ``inversion_result.nc``
     - temp
   * - ``gridded_posterior.nc``
     - temp
   * - ``visualization_notebook.ipynb``
     -