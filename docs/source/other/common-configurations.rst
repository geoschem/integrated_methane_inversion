Common IMI configurations
=========================

This page provides examples of how to configure the IMI `setup modules <../getting-started/imi-config-file.html#setup-modules>`__ 
and `run modules <../getting-started/imi-config-file.html#run-modules>`__ to accomplish some common tasks.

Default (preview) configuration
-------------------------------

By default the IMI will download the TROPOMI data for the period of interest, set up the template GEOS-Chem run directory, 
compute the prior emissions with a HEMCO standalone run, run the preview, and then stop. ::
    
    ## Setup modules
    ##   Turn on/off different steps in setting up the inversion 
    RunSetup: true
    SetupTemplateRundir: true
    SetupSpinupRun: false
    SetupJacobianRuns: false
    SetupInversion: false
    SetupPosteriorRun: false

    ## Run modules
    ##   Turn on/off different steps in performing the inversion
    DoHemcoPriorEmis: true
    DoSpinup: false
    DoJacobian: false
    ReDoJacobian: true
    DoInversion: false
    DoPosterior: false

    ## IMI preview
    ##   NOTE: RunSetup must be true to run preview
    DoPreview: true

If the results of the preview are satisfactory, you can try the next configuration on this page to run the inversion.
If they are not satisfactory, modify the configuration file (e.g., the region and/or period of interest) and try again.


Running an inversion after the preview
--------------------------------------

If the preview is complete and the results are satisfactory, you can proceed with the inversion (without re-running the preview or computing prior emissions). ::

    ## Inversion
    PrecomputedJacobian: false

    ## Setup modules
    ##   Turn on/off different steps in setting up the inversion 
    RunSetup: true
    SetupTemplateRundir: false
    SetupSpinupRun: true
    SetupJacobianRuns: true
    SetupInversion: true
    SetupPosteriorRun: true

    ## Run modules
    ##   Turn on/off different steps in performing the inversion
    DoHemcoPriorEmis: false
    DoSpinup: true
    DoJacobian: true
    ReDoJacobian: true
    DoInversion: true
    DoPosterior: true

    ## IMI preview
    ##   NOTE: RunSetup must be true to run preview
    DoPreview: false


Running a sensitivity inversion
-------------------------------

You've completed an initial inversion. Use the following configuration to run a new inversion with modified prior error (``PriorError``), 
observational error (``ObsError``), or regularization parameter (``Gamma``). ::

    ## Inversion
    PrecomputedJacobian: true
    ReferenceRunDir: "/path/to/your/run/dir"

    ## Setup modules
    ##   Turn on/off different steps in setting up the inversion 
    RunSetup: false
    SetupTemplateRundir: false
    SetupSpinupRun: false
    SetupJacobianRuns: false
    SetupInversion: false
    SetupPosteriorRun: false

    ## Run modules
    ##   Turn on/off different steps in performing the inversion
    DoHemcoPriorEmis: false
    DoSpinup: false
    DoJacobian: false
    ReDoJacobian: false
    DoInversion: true
    DoPosterior: false

    ## IMI preview
    ##   NOTE: RunSetup must be true to run preview
    DoPreview: false

Note that the final results of the original inversion (``inversion_result.nc`` and ``gridded_posterior.nc``) 
will be overwritten if not archived before running the sensitivity inversion.


Running an inversion without the preview
----------------------------------------

We generally don't recommend doing this, but if you wish to perform an inversion without manually inspecting the results 
of the IMI preview, use the following configuration to run the IMI from end to end, with a threshold on the expected degrees of
freedom for signal (DOFS) to cancel the inversion; if the expected DOFS are below the threshold, the IMI will exit with a warning. ::

    ## Setup modules
    ##   Turn on/off different steps in setting up the inversion 
    RunSetup: true
    SetupTemplateRundir: true
    SetupSpinupRun: true
    SetupJacobianRuns: true
    SetupInversion: true
    SetupPosteriorRun: true

    ## Run modules
    ##   Turn on/off different steps in performing the inversion
    DoHemcoPriorEmis: true
    DoSpinup: true
    DoJacobian: true
    ReDoJacobian: true
    DoInversion: true
    DoPosterior: true

    ## IMI preview
    ##   NOTE: RunSetup must be true to run preview
    DoPreview: true
    DOFSThreshold: {insert-threshold-value}


Modifying prior emission estimates
----------------------------------

**Set up the template run directory**
::

    ## Setup modules
    ##   Turn on/off different steps in setting up the inversion 
    RunSetup: true
    SetupTemplateRundir: true
    SetupSpinupRun: false
    SetupJacobianRuns: false
    SetupInversion: false
    SetupPosteriorRun: false

    ## Run modules
    ##   Turn on/off different steps in performing the inversion
    DoHemcoPriorEmis: false
    DoSpinup: false
    DoJacobian: false
    ReDoJacobian: false
    DoInversion: false
    DoPosterior: false

    ## IMI preview
    ##   NOTE: RunSetup must be true to run preview
    DoPreview: false

**Run the preview**

After :doc:`modifying the prior emission inventories <../advanced/custom-prior-emissions-hemco>`,
run the preview without setting up the template run directory. ::

    ## Setup modules
    ##   Turn on/off different steps in setting up the inversion 
    RunSetup: true
    SetupTemplateRundir: false
    SetupSpinupRun: false
    SetupJacobianRuns: false
    SetupInversion: false
    SetupPosteriorRun: false

    ## Run modules
    ##   Turn on/off different steps in performing the inversion
    DoHemcoPriorEmis: true
    DoSpinup: false
    DoJacobian: false
    ReDoJacobian: false
    DoInversion: false
    DoPosterior: false

    ## IMI preview
    ##   NOTE: RunSetup must be true to run preview
    DoPreview: true

If satisfied with the preview results, continue with one of the above configurations to run the inversion.
