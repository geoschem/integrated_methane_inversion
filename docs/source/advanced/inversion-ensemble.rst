Constructing an inversion ensemble
==================================

After performing an inversion, you can use the IMI to create a low-cost ensemble of sensitivity inversions with different 
inversion parameters. This is because the Jacobian matrix computed in the first inversion can easily be reused.

See the `Common configurations page <../other/common-configurations.html#running-a-sensitivity-inversion>`__ 
for instructions on how to re-configure the IMI to use a pre-computed Jacobian. Then modify
the values of ``PriorError``, ``ObsError``, and/or ``Gamma`` in the configuration file and re-run the inversion.

.. note::
    Make sure to archive the final results of the original inversion (``inversion_result.nc`` and ``gridded_posterior.nc``) 
    before running the sensitivity inversion. Those files will be overwritten.

If you want to run a sensitivity inversion with updated prior emission inventories, the pre-computed Jacobian needs
to be scaled according to the differences between the original and updated inventories. Instructions for this to come.