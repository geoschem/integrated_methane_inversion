IMI configuration file
======================
This page documents settings in the IMI configuration file (``config.yml``).

General
~~~~~~~
.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``isAWS``
     - Boolean for running the IMI on AWS (``true``) or a local cluster (``false``).
   * - ``RunName``
     - Name for this inversion; will be used for directory names and prefixes.
   * - ``UseSlurm``
     - Boolean for running the IMI as a batch job with ``sbatch`` instead of interactively.
       Select ``true`` to run the IMI with ``sbatch run_imi.sh``.
       Select ``false`` to run the IMI with ``./run_imi.sh`` (:doc:`via tmux <../advanced/running-with-tmux>`).

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
   * - ``NestedRegion``
     - Nesting domain for the inversion. 
       Select ``AS`` for Asia, ``EU`` for Europe, or ``NA`` for North America.
       See the `GEOS-Chem horizontal grids <http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_horizontal_grids>`_ documentation
       for details about the available nesting domains.

State vector 
~~~~~~~~~~~~
.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``CreateAutomaticRectilinearStateVectorFile``
     - Boolean for whether the IMI should automatically create a rectilinear state vector for the inversion. 
       If ``false``, a custom/pre-generated state vector netcdf file must be provided below.
   * - ``nBufferClusters``
     - Number of buffer elements (clusters of GEOS-Chem grid cells lying outside the region of interest) to add to the state vector 
       of emissions being optimized in the inversion.
   * - ``BufferDeg``
     - Width of the buffer elements, in degrees; will not be used if ``CreateAutomaticRectilinearStateVectorFile`` is ``false``.
   * - ``LandThreshold``
     - Land-cover fraction below which to exclude GEOS-Chem grid cells from the state vector when creating the state vector file.

Custom/pre-generated state vector
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These settings are only used if ``CreateAutomaticRectilinearStateVectorFile`` is ``false``. 
Use them to :doc:`create a custom state vector file <../advanced/custom-state-vector>` from a shapefile 
in conjunction with the ``statevector_from_shapefile.ipynb`` jupyter notebook located at::

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
     - Error in the prior estimates (1-sigma; relative). E.g., ``0.5`` means 50% error.
   * - ``ObsError``
     - Observational error (1-sigma; absolute; ppb). E.g., ``15`` means 15 ppb error.
   * - ``Gamma``
     - Regularization parameter; typically between 0 and 1.
   * - ``PrecomputedJacobian``
     - Boolean for whether the Jacobian matrix has already been computed (``true``) or not (``false``).

Grid
~~~~
.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``Res``
     - Resolution for inversion. Options are ``"0.25x0.3125"`` and ``"0.5x0.625"``.
   * - ``Met``
     - Meteorology to use for the inversion. Options are ``"geosfp"`` (for ``Res: "0.25x0.3125"``) and ``"merra2"`` (for ``Res: "0.5x0.625"``).
   * - ``HalfPolar``
     - Set to ``"F"`` for nested grid simulations. 
   * - ``Levs``
     - Set to ``"47"`` to use the reduced 47-level grid recommended for CH4 simulations.
   * - ``NestedGrid``
     - ``"T"`` for nested simulations.
   * - ``Buffer``
     - Width of GEOS-Chem meteorological buffer zone (different from state vector buffer clusters above). 
       For nested simulations, the recommended input is ``"3 3 3 3"`` to use 3 grid cells along the edges of the nested-grid domain as buffer zone.

Setup modules
~~~~~~~~~~~~~
These settings turn on/off (``true`` / ``false``) different steps for setting up the IMI.

.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``SetupTemplateRundir``
     - Copy run directory files from GEOS-Chem and modify with information from ``config.yml``.
   * - ``SetupSpinupRun``
     - Setup the run directory for the spin-up simulation.
   * - ``SetupJacobianRuns``
     - Setup run directories for N+1 simulations (one reference simulation, plus N sensitivity simulations for the N state vector elements). 
       Output of these simulations will be used to construct the Jacobian.
   * - ``SetupInversion``
     - Setup the inversion directory containing scripts needed to perform the inverse analysis; inversion results will be saved here.
   * - ``SetupPosteriorRun``
     - Setup the run directory for the posterior simulation.

Run modules
~~~~~~~~~~~
These settings turn on/off (``true`` / ``false``) different steps of the inversion.

.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``RunSetup``
     - Run the setup script (``setup_imi.sh``), including selected setup modules above.
   * - ``DoSpinup``
     - Run the spin-up simulation.
   * - ``DoJacobian``
     - Run the reference and sensitivity simulations.
   * - ``DoInversion``
     - Run the inverse analysis code.
   * - ``DoPosterior``
     - Run the posterior simulation.

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

Compute resources to request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These settings are specific to Harvard's Cannon compute cluster. Not used for cloud runs.

.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``nCPUs``
     - Number of cpus to use in ``sbatch`` scripts.
   * - ``partition``
     - Name of the cluster partition to use with ``sbatch`` (eg. ``"huce_cascade"``).
