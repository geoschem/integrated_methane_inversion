IMI preview
===========

The IMI preview feature allows users to estimate the quality and information content of a proposed inversion 
without actually performing the inversion.

Under the default configuration, the IMI performs only the preview and then stops. This is to prevent 
accidental initiation of low-quality but potentially expensive inversions. For more details about
the preview (default) configuration, see the 
`Common configurations page <../other/common-configurations.html#default-preview-configuration>`__.

To run the preview after selecting a region and time period of interest in the :doc:`configuration file <imi-config-file>` 
(and modifying any other configurable settings), simply `run the IMI <quick-start.html#run-the-imi>`__ with the ``DoPreview``
configuration option set to ``true``.

The IMI preview provides the following information for users to assess their proposed inversion (as it is 
described in the configuration file):

  - Map of mean TROPOMI observations for the region and period of interest
  - Map of prior emission estimates to be used in the inversion
  - Map of observation density for the region and period of interest
  - Map of mean SWIR albedo for the region and period of interest
  - Total number of observations available in the region of interest during the inversion period
  - Rough estimate of degrees of freedom for signal (DOFS) for the inversion
  - Rough estimate of USD financial cost of the inversion

This information is generated as a ``.txt`` file and collection of ``.png`` files in the preview directory, 
which is located at::

    /home/ubuntu/imi_output_dir/{YourRunName}/preview_run/

The ``.txt`` file can be viewed directly in the terminal. To view the ``.png`` files, first download them from
EC2 to your local computer using::

    $ scp -i /local/path/to/my-key-pair.pem ubuntu@my-instance-public-dns-name:/path/to/my-file.png /local/path/to/my-file.png

For more informaton on this command, see the 
`AWS Documentation <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AccessingInstancesLinux.html>`_.

**Tips for interpreting the IMI preview results:**

  - Inspect the maps of XCH4 observations and observation density to evaluate TROPOMI coverage for the 
    region and period of interest.
  - Compare the maps of XCH4 observations and prior emission estimates to evaluate spatial correspondence 
    between the two datasets. 
  - Compare the maps of XCH4 observations and SWIR albedo to confirm that there are no obvious albedo-related 
    artifacts in the methane retrieval field, which would be diagnosed by high spatial correlation between 
    XCH4 and albedo.
  - DOFS > 1 is a bare minimum to achieve any solid information about emissions. 
  - DOFS < 2 is marginal for most applications.
  - If there is an obvious mismatch between the XCH4 observations and the prior emissions (for example due 
    to severe bias in the prior inventory) OR if the expected DOFS is low, consider: (a) using an improved 
    prior inventory, (b) increasing the inversion period to incorporate more observations, and/or 
    (c) increasing the prior error estimate.
  - If there is indication of albedo-related artifacts in the XCH4 field, consider removing the affected
    observations. This can be done by modifying the TROPOMI data filters via the ``filter_tropomi()``
    function in ``/home/ubuntu/integrated_methane_inversion/src/inversion_scripts/utils.py``.