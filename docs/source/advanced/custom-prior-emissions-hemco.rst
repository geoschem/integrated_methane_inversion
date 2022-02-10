Modifying prior emission estimates
==================================

To modify the default prior emission inventories, first generate the template run directory following
the `configuration instructions for modifying emissions <../other/common-configurations.html#modifying-prior-emission-estimates>`__.

Once the template run directory is ready, you will need to modify the emission inventories via 
`HEMCO <http://wiki.seas.harvard.edu/geos-chem/index.php/HEMCO>`_.

Start by transferring your custom emission inventory to EC2::

    $ scp -i /local/path/to/my-key-pair.pem /local/path/to/my-inventory.nc ubuntu@my-instance-public-dns-name:/path/to/my-inventory.nc

The emissions need to be defined in a netcdf file formatted as HEMCO expects. 

TODO