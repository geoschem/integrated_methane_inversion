Creating a custom state vector file
===================================

By default the IMI uses `latitude/longitude bounds <../getting-started/imi-config-file.html#region-of-interest>`__ 
to automatically create a gridded state vector file for a rectilinear region of interest with surrounding buffer elements.

The state vector file is located at ::

    /home/ubuntu/imi_output_dir/{YourRunName}/StateVector.nc

It contains the state variable labels for every grid cell in the inversion domain. For example, if the region of interest 
contains 200 emission elements and the IMI is configured to use 8 additional buffer elements, then the total number of state 
variables is 208 and the state vector file will assign a number between 1 and 208 to every grid cell in the inversion domain.

Instead of the default rectilinear region of interest, you may want to use an irregular region as was done for the Permian 
Basin by Varon et al. (2022; link_TODO). To do so you will need to generate the state vector file yourself.

The easiest way to do this is by using a shapefile for the region of interest in conjunction with the
``statevector_from_shapefile.ipynb`` jupyter notebook. The notebook is located at ::

    /home/ubuntu/integrated_methane_inversion/src/notebooks/statevector_from_shapefile.ipynb

First upload a shapefile for the custom region of interest to your EC2 instance::

    $ scp -i /local/path/to/my-key-pair.pem /local/path/to/my-shapefile ubuntu@my-instance-public-dns-name:/path/to/my-shapefile

Next, open the configuration file and insert the path to your shapefile in the
`custom state vector section <../getting-started/imi-config-file.html#custom-pre-generated-state-vector>`__.
Also provide the latitude/longitude bounds for the desired inversion domain, which will include both the irregular region of
interest and the additional coarse buffer elements.

Next, follow `these short instructions <https://docs.aws.amazon.com/dlami/latest/devguide/setup-jupyter.html>`_ to set up and connect to
a jupyter notebook server on AWS. Once connected to the server, open ``statevector_from_shapefile.ipynb`` and run its contents to
generate a state vector file from your shapefile.

If no shapefile is available, you will need to construct the custom state vector file manually. You may want to start from an 
automatically generated rectilinear state vector file. 