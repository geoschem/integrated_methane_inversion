IMI directory contents 
======================

This page describes the contents of various file directories generated and populated by the IMI in the course of an inversion.

Inversion directory
-------------------

The inversion directory is where the IMI computes the Jacobian, obtains the optimal estimate of emissions, and saves the results.

It is located at ``/home/ubuntu/CH4_Workflow/{YourRunName}/inversion``.

In addition to a shell script and several Python scripts used in the inversion, 
the following items will appear in the inversion directory after completing an inversion:

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
     - temp