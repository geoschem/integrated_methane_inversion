Common IMI configurations
=========================

This page provides examples of how to configure the IMI `setup modules <../getting-started/imi-config-file.html#setup-modules>`__ 
and `run modules <../getting-started/imi-config-file.html#run-modules>`__ to accomplish some common tasks.

Default (preview) configuration
-------------------------------

By default the IMI will download the TROPOMI data for the period of interest, set up the template GEOS-Chem run directory, 
run the preview, and then stop. ::
    
    ## Setup modules
    SetupTemplateRundir: true
    SetupSpinupRun: false
    SetupJacobianRuns: false
    SetupInversion: false
    SetupPosteriorRun: false
    
    ## Run modules
    RunSetup: true
    DoSpinup: false
    DoJacobian: false
    DoInversion: false
    DoPosterior: false
    
    ## IMI preview
    DoPreview: true

If the results of the preview are satisfactory, try the next configuration on this page to run the inversion.
If they are not satisfactory, modify the configuration file (e.g., the region and/or period of interest) and try again.


Inversion after successful preview
----------------------------------

The preview is complete and the results are satisfactory. Now proceed with the inversion (without re-running the preview). ::

    ## Inversion
    PrecomputedJacobian: false

    ## Setup modules
    SetupTemplateRundir: false
    SetupSpinupRun: true
    SetupJacobianRuns: true
    SetupInversion: true
    SetupPosteriorRun: true
    
    ## Run modules
    RunSetup: true
    DoSpinup: true
    DoJacobian: true
    DoInversion: true
    DoPosterior: true
    
    ## IMI preview
    DoPreview: false


Sensitivity inversion
---------------------

An initial inversion is complete. Now re-run the inversion with modified prior error (``PriorError``), 
observational error (``ObsError``), or regularization parameter (``Gamma``). ::

    ## Inversion
    PrecomputedJacobian: true

    ## Setup modules
    SetupTemplateRundir: false
    SetupSpinupRun: false
    SetupJacobianRuns: false
    SetupInversion: false
    SetupPosteriorRun: false
    
    ## Run modules
    RunSetup: false
    DoSpinup: false
    DoJacobian: false
    DoInversion: true
    DoPosterior: false
    
    ## IMI preview
    DoPreview: false

Note that the final results of the original inversion (``inversion_result.nc`` and ``gridded_posterior.nc``) 
will be overwritten if not archived before running the sensitivity inversion.