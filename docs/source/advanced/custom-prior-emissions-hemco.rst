Modifying prior emission estimates
==================================

To modify the default prior emission inventories, first generate the template run directory following
the `configuration instructions for modifying emissions <../other/common-configurations.html#modifying-prior-emission-estimates>`__.

Once the template run directory is ready, you will need to modify the emission inventories via 
`HEMCO <http://wiki.seas.harvard.edu/geos-chem/index.php/HEMCO>`_.

Start by transferring your custom emission inventory to EC2::

    $ scp -i /local/path/to/my-key-pair.pem /local/path/to/my-inventory.nc ubuntu@my-instance-public-dns-name:/path/to/my-inventory.nc

The emissions need to be defined in a netcdf file formatted as HEMCO expects. See the 
`HEMCO documentation for preparing data files <http://wiki.seas.harvard.edu/geos-chem/index.php/Preparing_data_files_for_use_with_HEMCO>`_
for details on how to format your custom emission inventory for use with HEMCO.

Once your inventory has been properly formatted, you can include it as an emission field via HEMCO. To do this, navigate to the template 
run directory and open the HEMCO configuration file with vim (``vi``) or emacs::

    $ cd /home/ubuntu/imi_output_dir/{YourRunName}/template_run
    $ emacs HEMCO_Config.rc

Follow `instructions in the HEMCO User's Guide <http://wiki.seas.harvard.edu/geos-chem/index.php/The_HEMCO_User%27s_Guide#Base_emissions>`_
to add a new emission field.

You can run the :doc:`IMI preview <../getting-started/imi-preview>` to quickly check that the updated emissions are working as expected. 