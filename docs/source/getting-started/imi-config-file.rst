IMI configuration file
======================
This page documents settings in the IMI configuration file (``config.yml``).

General
~~~~~~~
.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``RunName``
     - Name for this inversion; will be used for directory names and prefixes.
   * - ``isAWS``
     - Boolean for running the IMI on AWS (``true``) or a local cluster (``false``).
   * - ``UseSlurm``
     - Boolean for running the IMI as a batch job with ``sbatch`` instead of interactively.
       Select ``true`` to run the IMI with ``sbatch run_imi.sh``.
       Select ``false`` to run the IMI with ``./run_imi.sh`` (:doc:`via tmux <../advanced/running-with-tmux>`).
   * - ``SafeMode``
     - Boolean for running in safe mode to prevent overwriting existing files.

Period of interest
~~~~~~~~~~~~~~~~~~
.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``StartDate``
     - Beginning of the inversion period in ``YYYYMMDD`` format (this date is included in the inversion, 0-24h UTC).
   * - ``EndDate``
     - End of the inversion period in ``YYYYMMDD`` format (this date is excluded from the inversion, 0-24h UTC).
   * - ``SpinupMonths``
     - Number of months for the spinup simulation. 

TROPOMI data type
~~~~~~~~~~~~~~~~~~
.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``BlendedTROPOMI``
     - Boolean for if the Blended TROPOMI+GOSAT data should be used (``true``) or if the operational data should be used (``false``).

Region of interest
~~~~~~~~~~~~~~~~~~
.. list-table::
   :widths: 30, 70
   :class: tight-table 

   * - ``LonMin``
     - Minimum longitude edge of the region of interest (only used if ``CreateAutomaticRectilinearStateVectorFile`` is ``true``).
   * - ``LonMax``
     - Maximum longitude edge of the region of interest (only used if ``CreateAutomaticRectilinearStateVectorFile`` is ``true``).
   * - ``LatMin``
     - Minimum latitude edge of the region of interest (only used if ``CreateAutomaticRectilinearStateVectorFile`` is ``true``).
   * - ``LatMax``
     - Maximum latitude edge of the region of interest (only used if ``CreateAutomaticRectilinearStateVectorFile`` is ``true``).
   * - ``NestedGrid``
     - Boolean for using the GEOS-Chem nested grid simulation. Must be
       ``true`` for IMI regional inversions.
   * - ``NestedRegion``
     - Nesting domain for the inversion. Select ``AF`` for Africa, ``AS`` for Asia, ``EU`` for Europe, ``ME`` for the Middle East, ``NA`` for North America, ``OC`` for Oceania, ``RU`` for Russia, or ``SA`` for South America. For global met fields set this option to ``""`` See the `GEOS-Chem horizontal grids <http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_horizontal_grids>`_ documentation for details about the available nested-grid domains.

State vector 
~~~~~~~~~~~~
.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``CreateAutomaticRectilinearStateVectorFile``
     - Boolean for whether the IMI should automatically create a rectilinear state vector for the inversion. If ``false``, a custom/pre-generated state vector netcdf file must be provided below.
   * - ``nBufferClusters``
     - Number of buffer elements (clusters of GEOS-Chem grid cells lying outside the region of interest) to add to the state vector of emissions being optimized in the inversion. Default value is ``8``.
   * - ``BufferDeg``
     - Width of the buffer elements, in degrees; will not be used if ``CreateAutomaticRectilinearStateVectorFile`` is ``false``. Default is ``5`` (~500 km).
   * - ``LandThreshold``
     - Land-cover fraction below which to exclude GEOS-Chem grid cells from the state vector when creating the state vector file. Default value is ``0.25``.
   * - ``OffshoreEmisThreshold``
     - Offshore GEOS-Chem grid cells with oil/gas emissions above this threshold will be included in the state vector. Default value is ``0``.

Clustering Options
^^^^^^^^^^^^^^^^^^
For more information on using the clustering options take a look at the `clustering options page <../advanced/using-clustering-options.html>`__.

.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``ReducedDimensionStateVector``
     - Boolean for whether to reduce the dimension of the statevector from the native resolution version by clustering elements. If ``false`` the native state vector is used with no dimension reduction.
   * - ``ClusteringMethod``
     - Clustering method to use for state vector reduction. (eg. "kmeans" or "mini-batch-kmeans")
   * - ``NumberOfElements``
     - Number of elements in the reduced dimension state vector. This is only used if ``ReducedDimensionStateVector`` is ``true``.
   * - ``ForcedNativeResolutionElements``
     - yaml list of of coordinates that you would like to force as native resolution state vector elements [lat, lon]. This is useful for ensuring hotspot locations are at the highest available resolution. 

Custom/pre-generated state vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These settings are only used if ``CreateAutomaticRectilinearStateVectorFile`` is ``false``. Use them to :doc:`create a custom state vector file <../advanced/custom-state-vector>` from a shapefile in conjunction with the ``statevector_from_shapefile.ipynb`` jupyter notebook located at::

  $ /home/ubuntu/integrated_methane_inversion/src/notebooks/statevector_from_shapefile.ipynb

.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``StateVectorFile``
     - Path to the custom or pre-generated state vector netcdf file. File will be saved here if generating it from a shapefile.
   * - ``ShapeFile``
     - Path to the shapefile.

Note: To setup a remote Jupyter notebook check out the quick start guide `visualize results with python <../getting-started/quick-start.html#visualize-results-with-python>`__ section.

Inversion
~~~~~~~~~
.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``PriorError``
     - Error in the prior estimates (1-sigma; relative). Default is ``0.5`` (50%) error.
   * - ``ObsError``
     - Observational error (1-sigma; absolute; ppb). Default value is ``15`` ppb error.
   * - ``Gamma``
     - Regularization parameter; typically between 0 and 1. Default value is ``1.0``.
   * - ``PrecomputedJacobian``
     - Boolean for whether the Jacobian matrix has already been computed (``true``) or not (``false``). Default value is ``false``.

Grid
~~~~
.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``Res``
     - Resolution for inversion. Options are ``"0.25x0.3125"`` and ``"0.5x0.625"``.
   * - ``Met``
     - Meteorology to use for the inversion. Options are ``"geosfp"`` (for ``Res: "0.25x0.3125"``) and ``"merra2"`` (for ``Res: "0.5x0.625"``).

Setup modules
~~~~~~~~~~~~~
These settings turn on/off (``true`` / ``false``) different steps for setting up the IMI.

.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``SetupTemplateRundir``
     - Boolean to create a GEOS-Chem run directory and modify it with settings from ``config.yml``.
   * - ``SetupSpinupRun``
     - Boolean to set up a run directory for the spinup-simulation by copying the template run directory and modifying the start/end dates, restart file, and diagnostics.
   * - ``SetupJacobianRuns``
     - Boolean to set up run directories for N+1 simulations (one reference simulation, plus N sensitivity simulations for the N state vector elements) by copying the template run directory and modifying the start/end dates, restart file, and diagnostics. Output from these simulations will be used to construct the Jacobian.
   * - ``SetupInversion``
     - Boolean to set up the inversion directory containing scripts needed to perform the inverse analysis; inversion results will be saved here.
   * - ``SetupPosteriorRun``
     - Boolean to set up the run directory for the posterior simulation by copying the template run directory and modifying the start/end dates, restart file, and diagnostics.

Run modules
~~~~~~~~~~~
These settings turn on/off (``true`` / ``false``) different steps for running the inversion.

.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``RunSetup``
     - Boolean to run the setup script (``setup_imi.sh``), including selected setup modules above.
   * - ``DoSpinup``
     - Boolean to run the spin-up simulation.
   * - ``DoJacobian``
     - Boolean to run the reference and sensitivity simulations.
   * - ``DoInversion``
     - Boolean to run the inverse analysis code.
   * - ``DoPosterior``
     - Boolean to run the posterior simulation.

SLURM Resource Allocation
~~~~~~~~~~~~~~~~~~~~~~~~~
These settings are used to allocate resources (CPUs and Memory) to the different simulations needed to run the inversion.
Note: some python scripts are also deployed using slurm and default to using the ``SimulationCPUs`` and ``SimulationMemory`` settings.

.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``RequestedTime``
     - Max amount of time to allocate to each sbatch job (eg. "0-6:00")
   * - ``SimulationCPUs``
     - Number of cores to allocate to each in series simulation.
   * - ``SimulationMemory``
     - Amount of memory to allocate to each in series simulation (in MB).
   * - ``JacobianCPUs``
     - Number of cores to allocate to each jacobian simulation (run in parallel).
   * - ``JacobianMemory``
     - Amount of memory to allocate to each jacobian simulation (in MB).
   * - ``SchedulerPartition``
     - Name of the partition(s) you would like all slurm jobs to run on (eg. "debug,huce_intel,seas_compute,etc").
   
IMI preview
~~~~~~~~~~~
.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``DoPreview``
     - Boolean to run the :doc:`IMI preview <imi-preview>` (``true``) or not (``false``).
   * - ``DOFSThreshold``
     - Threshold for estimated DOFS below which the IMI should automatically exit with a warning after performing the preview.
       Default value ``0`` prevents exit.

Advanced settings: GEOS-Chem options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These settings are intended for advanced users who wish to modify additional GEOS-Chem options.

.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``PerturbValue``
     - Value to perturb emissions by in each sensitivity simulation. Default value is ``1.5``.
   * - ``UseEmisSF``
     - Boolean to apply emissions scale factors derived from a previous inversion. This file should be provided as a netCDF file and specified in HEMCO_Config.rc. Default value is ``false``.
   * - ``UseOHSF``
     - Boolean to apply OH scale factors derived from a previous inversion. This file should be provided as a netCDF file and specified in HEMCO_Config.rc. Default value is ``false``.
   * - ``HourlyCH4``
     - Boolean to save out hourly diagnostics from GEOS-Chem. This output is used in satellite operators via post-processing. Default value is ``true``.
   * - ``PLANEFLIGHT``
     - Boolean to save out the planeflight diagnostic in GEOS-Chem. This output may be used to compare GEOS-Chem against planeflight data. The path to those data must be specified in input.geos. See the `planeflight diagnostic <http://wiki.seas.harvard.edu/geos-chem/index.php/Planeflight_diagnostic>`_ documentation for details. Default value is ``false``.
   * - ``GOSAT``
     - Boolean to turn on the GOSAT observation operator in GEOS-Chem. This will save out text files comparing GEOS-Chem to observations, but has to be manually incorporated into the IMI. Default value is ``false``.
   * - ``TCCON``
     - Boolean to turn on the TCCON observation operator in GEOS-Chem. This will save out text files comparing GEOS-Chem to observations, but has to be manually incorporated into the IMI. Default value is ``false``.
   * - ``AIRS``
     - Boolean to turn on the AIRS observation operator in GEOS-Chem. This will save out text files comparing GEOS-Chem to observations, but has to be manually incorporated into the IMI. Default value is ``false``.

Advanced settings: Local cluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These settings are intended for advanced users who wish to (:doc:`run
the IMI on a local cluster<../advanced/local-cluster>`).

.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``OutputPath``
     - Path for IMI runs and output.
   * - ``DataPath``
     - Path to GEOS-Chem input data.
   * - ``DataPathTROPOMI``
     - Path to TROPOMI input data.
   * - ``CondaFile``
     - Path to file containing Conda environment settings.
   * - ``CondaEnv``
     - Name of conda environment.
   * - ``RestartDownload``
     - Boolean for downloading an initial restart file from AWS S3. Default value is ``true``.
   * - ``RestartFilePrefix``
     - Path to initial GEOS-Chem restart file plus file prefix (e.g. ``GEOSChem.BoundaryConditions.`` or ``GEOSChem.Restart.``). The date string and file extension (``YYYYMMDD_0000z.nc4``) will be appended. This file will be used to initialize the spinup simulation.
   * - ``RestartFilePreviewPrefix``
     - Path to initial GEOS-Chem restart file plus file prefix (e.g. ``GEOSChem.BoundaryConditions.`` or ``GEOSChem.Restart.``). The date string and file extension (``YYYYMMDD_0000z.nc4``) will be appended. This file will be used to initialize the preview simulation.
   * - ``BCpath``
     - Path to GEOS-Chem boundary condition files (for nested grid simulations).
   * - ``BCversion``
     - Version of TROPOMI smoothed boundary conditions to use (e.g. ``v2023-04``). Note: this will be appended onto BCpath as a subdirectory.
   * - ``PreviewDryRun``
     - Boolean to download missing GEOS-Chem data for the preview run. Default value is ``true``.
   * - ``SpinupDryRun``
     - Boolean to download missing GEOS-Chem data for the spinup simulation. Default value is ``true``.
   * - ``ProductionDryRun``
     - Boolean to download missing GEOS-Chem data for the production (i.e. Jacobian) simulations. Default value is ``true``.
   * - ``PosteriorDryRun``
     - Boolean to download missing GEOS-Chem data for the posterior simulation. Default value is ``true``.
   * - ``BCDryRun``
     - Boolean to download missing GEOS-Chem data for the preview run. Default value is ``true``.
   * - ``PreviewDryRun``
     - Boolean to download missing GEOS-Chem boundary condition files. Default value is ``true``.

Note for ``*DryRun`` options: If you are running on AWS, you will be charged if your ec2 instance is not in the us-east-1 region. If running on a local cluster you must have AWS CLI enabled or you can modify the ``./download_data.py`` commands in ``setup_imi.sh`` to use ``washu`` instead of ``aws``. See the `GEOS-Chem documentation <https://geos-chem.readthedocs.io/en/latest/inputs/dry-run.html>`_ for more details.
