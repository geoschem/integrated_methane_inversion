IMI Configuration File
======================
This page documents and explains the various settings for the configuration file (``config.yml``).

General
~~~~~~~
- ``isAWS``: Boolean value of whether you are running the workflow on AWS (``true``) or a local cluster (``false``).
- ``RunName``: Name for your workflow that will be used for directory names and prefixes.
- ``UseSlurm``: Boolean value of whether you are running the end to end script as a batch job with ``sbatch`` (eg. ``sbatch run_ch4_inversion.sh``) instead of interactively (eg. ./run_ch4_inversion.sh)

Period of Interest
~~~~~~~~~~~~~~~~~~
- ``StartDate``: Beginning of the inversion time frame, using ``YYYYMMDD`` format.
- ``EndDate``: End of the inversion time frame, using ``YYYYMMDD`` format.
- ``SpinupMonths``: Number of months for the spinup simulation.

Region of Interest
~~~~~~~~~~~~~~~~~~
- ``LonMin``: Minimum longitude edge of your domain.
- ``LonMax``: Maximum longitude edge of your domain.
- ``LatMin``: Minimum latitude edge of your domain.
- ``LatMax``: Maximum latitude edge of your domain.

Inversion
~~~~~~~~~
- ``PriorError``: TODO
- ``ObsError``: TODO
- ``Gamma``: TODO

Grid
~~~~
- ``Res``: Resolution for simulation. Options are "4x5", "2x2.5", "0.5x0.625", and "0.25x0.3125".
- ``Met``: Meteorology to use for the simulation. Options are "merra2" or "geosfp".
- ``HalfPolar``: Set to "T" for global simulations to use half polar boxes. Set to "F" for nested grid simulations.
- ``Levs``: Set to 47 to use the reduced 47-level grid recommended for CH4 simulations.
- ``NestedGrid``: Set to "F" for global simulations or "T" for nested simulations.
- ``REGION``: Set to "" for global or "NA", "AS", "CH", "EU" for default domains.
- ``Buffer``: Set to "0 0 0 0" for global simulations. For nested simulations, the recommendation is to use "3 3 3 3" to use 3 grid cells along the nested-grid domain for your buffer zone.

Setup modules
~~~~~~~~~~~~~
These settingd give options to turn on/off different steps in setting up GEOS-Chem and the inversion.

- ``CreateStateVectorFile``: Create a netCDF cluster file containing the gridboxes or regions where emissions will be perturbed. If this is set to false, the ClusterFile option must be specified further down.
- ``SetupTemplateRundir``: Copy run directory files from GEOS-Chem and replace text in input files according to settings in ``config.yml``.
- ``SetupSpinupRun``: Create a run directory to spinup a new restart file representative of your model setup.
- ``SetupJacobianRuns``: Setup run directories for each of your perturbation clusters. The output from these simulations will be used to construct the Jacobian.
- ``SetupInversion``: Copy scripts used to post-process GEOS-Chem data, build the Jacobian, and run the inversion.
- ``SetupPosteriorRun``: Create a run directory to submit a posterior run.

Run modules
~~~~~~~~~~~
Turn on/off different steps in performing the inversion.

- ``RunSetup``: Run the setup script to create all necessary run directories for the workflow.
- ``DoSpinup``: Set to true to run a spinup simulation to generate a new restart file representative of your model setup.
- ``DoJacobian``: Set to true to submit the jacobian runs.
- ``DoInversion``: Set to true to run the inversion.
- ``DoPosterior``: Set to true to run the posterior simulation.

State vector settings 
~~~~~~~~~~~~~~~~~~~~~
Only used if ``CreateStateVectorFile`` is set to true

- ``BufferDeg``: TODO
- ``nBufferClusters``: TODO
- ``LandThreshold``: TODO

Custom Statevector File
~~~~~~~~~~~~~~~~~~~~~~~
Only used if CreateStateVectorFile is set to false.

- ``StateVectorFile``: Path to your custom state vector file
- ``ShapeFile``: Path to your shapefile.
- ``LonMinCustomStateVector``: Minimum longitude edge of your state vector file.
- ``LonMaxCustomStateVector``: Maximum longitude edge of your state vector file.
- ``LatMinCustomStateVector``: Minimum latitude edge of your state vector file.
- ``LatMaxCustomStateVector``: Maximum latitude edge of your state vector file.

Compute Resources to Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Note: Specific to Harvard's cannon cluster. Not used for cloud runs.

- ``nCPUs``: Number of cpus to use in sbatch scripts
- ``partition``: Name of the cluster partition to use with sbatch (eg. "huce_cascade")
