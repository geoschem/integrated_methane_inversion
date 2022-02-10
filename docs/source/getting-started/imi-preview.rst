IMI preview
===========

The IMI preview feature allows users to project the quality and information content of a proposed inversion 
before actually performing the inversion. It takes approximately 10 minutes to run for a 1-year inversion
performed on a ``c5.9xlarge`` instance with 36 CPUs.

Under the default configuration, the IMI performs only the preview and then stops. This is to prevent 
the accidental initiation of low-quality but potentially expensive inversions. For more details about
the preview (default) configuration, see `<../other/common-configurations.html#default-(preview)-configuration>`__.

The IMI preview provides the following information for users to assess their proposed inversion (as configured
in the configuration file):

  - Map of mean TROPOMI observations for the region and period of interest
  - Map of prior emissions estimates to be used in the inversion
  - Map of observation density for the region and period of interest
  - Map of mean albedo for the region and period of interest
  - Total number of observations available in the region of interest during the inversion period
  - Rough estimate of degrees of freedom for signal (DOFS) for the inversion
  - Rough estimate of USD financial cost of the inversion

This information is generated as a collection of `.png` and `.txt` files in the preview directory located at::

    $ /home/ubuntu/imi_output_dir/{YourRunName}/preview_run/

