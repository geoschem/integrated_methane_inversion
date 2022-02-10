IMI preview
===========

The IMI preview feature allows users to project the quality and information content of a proposed inversion 
before actually performing the inversion. It takes approximately 10 minutes to run for a 1-year inversion
performed on a ``c5.9xlarge`` instance with 36 CPUs.

Under the default configuration, the IMI performs only the preview and then stops. This is to prevent 
the accidental initiation of low-quality but potentially expensive inversions. For more details about
the preview (default) configuration, see `<../other/common-configurations.html#default-(preview)-configuration>`__.

The IMI preview provides the following information for users to assess their proposed inversion (as described
in the configuration file):

  - Map of mean TROPOMI observations for the region and period of interest
  - Map of prior emissions estimates to be used in the inversion
  - Map of observation density for the region and period of interest
  - Map of mean SWIR albedo for the region and period of interest
  - Total number of observations available in the region of interest during the inversion period
  - Rough estimate of degrees of freedom for signal (DOFS) for the inversion
  - Rough estimate of USD financial cost of the inversion

This information is generated as a `.txt` file and collection of `.png` files in the preview directory, 
which is located at::

    $ /home/ubuntu/imi_output_dir/{YourRunName}/preview_run/

The `.txt` file can be viewed directly in the terminal. To view the `.png` files, first download them from
your EC2 instance to your local computer using::

    $ scp -i /local/path/to/my-key-pair.pem ubuntu@my-instance-public-dns-name:/path/to/my-file.png /local/path/to/my-file.png

Some tips for interpreting the IMI preview results:

  - Inspect the maps of XCH4 observations and observation density to evaluate TROPOMI coverage for the 
    region and period of interest.
  - Compare the maps of XCH4 observations and prior emission estimates to evaluate spatial correspondence 
    between the two datasets. 
  - Compare the maps of XCH4 observations and SWIR albedo to confirm that there are no obvious albedo-related 
    artifacts in the methane retrieval field, which would be diagnosed by high spatial correlation between 
    the XCH4 and albedo fields.
  - DOFS > 1 is a bare minimum to achieve any solid information about emissions. 
  - DOFS < 2 is marginal for most applications.
  - If there is an obvious mismatch between the XCH4 observations and the prior emissions (for example due 
    to severely bias in the prior inventory) OR if the expected DOFS is low, consider: (a) using an improved 
    prior inventory, (b) increasing the inversion period to incorporate more observations, and/or 
    (c) increasing the prior error estimate.