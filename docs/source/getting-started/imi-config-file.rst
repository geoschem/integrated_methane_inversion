IMI Configuration File
======================
This page documents settings in the IMI configuration file (``config.yml``).

General
~~~~~~~
.. list-table::
   :widths: 40, 60
   :class: tight-table

   * - ``isAWS``
     - Boolean for running the IMI on AWS (``true``) or a local cluster (``false``).
   * - ``RunName``
     - Name for this inversion; will be used for directory names and prefixes.
   * - ``UseSlurm``
     - Boolean for running the IMI as a batch job with ``sbatch`` instead of interactively.
       Select ``true`` to run the IMI with ``sbatch run_ch4_inversion.sh``.
       Select ``false`` to run the IMI with ``./run_ch4_inversion.sh`` (via tmux_link_TODO).

Period of Interest
~~~~~~~~~~~~~~~~~~
.. list-table::
   :widths: 40, 60
   :class: tight-table

   * - ``StartDate``
     - Beginning of the inversion period in ``YYYYMMDD`` format (this date is included in the inversion, 0-24h UTC).
   * - ``EndDate``
     - End of the inversion period in ``YYYYMMDD`` format (this date is excluded from the inversion, 0-24h UTC).
   * - ``SpinupMonths``
     - Number of months for the spinup simulation. 

Region of Interest
~~~~~~~~~~~~~~~~~~
.. list-table::
   :widths: 40, 60
   :class: tight-table 

   * - ``LonMin``
     - Minimum longitude edge of the region of interest (only used if ``CreateStateVectorFile`` is ``true``).
   * - ``LonMax``
     - Maximum longitude edge of the region of interest (only used if ``CreateStateVectorFile`` is ``true``).
   * - ``LatMin``
     - Minimum latitude edge of the region of interest (only used if ``CreateStateVectorFile`` is ``true``).
   * - ``LatMax``
     - Maximum latitude edge of the region of interest (only used if ``CreateStateVectorFile`` is ``true``).
   * - ``REGION``
     - Nesting domain for the inversion. 
       Select ``AS`` for Asia, ``EU`` for Europe, or ``NA`` for North America.
       See the `GEOS-Chem horizontal grids <http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_horizontal_grids>`_ documentation
       for details about the available nesting domains.

State vector 
~~~~~~~~~~~~
.. list-table::
   :widths: 40, 60
   :class: tight-table

   * - ``CreateStateVectorFile``
     - Boolean for whether the IMI should automatically create a rectilinear state vector for the inversion. 
       If ``false``, a custom/pre-generated state vector netcdf file must be provided below.
   * - ``nBufferClusters``
     - Number of buffer elements (clusters of GEOS-Chem grid cells lying outside the region of interest) to add to the state vector 
       of emissions being optimized in the inversion.
   * - ``BufferDeg``
     - Width of the buffer elements, in degrees; will not be used if ``CreateStateVectorFile`` is ``false``.
   * - ``LandThreshold``
     - Land-cover fraction below which to exclude GEOS-Chem grid cells from the state vector when creating the state vector file.

Custom/pre-generated state vector file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These settings are only used if ``CreateStateVectorFile`` is ``false``.

.. list-table::
   :widths: 40, 60
   :class: tight-table

   * - ``StateVectorFile``
     - Path to the custom or pre-generated state vector netcdf file.

To create a custom state vector file from a shapefile, use the following settings in conjunction with the ``statevector_from_shapefile.ipynb`` jupyter notebook.

.. list-table::
   :widths: 40, 60
   :class: tight-table

   * - ``ShapeFile``
     - Path to the shapefile.
   * - ``LonMinCustomStateVector``
     - Minimum longitude edge of the desired inversion domain (includes region of interest and buffer elements).
   * - ``LonMaxCustomStateVector``
     - Maximum longitude edge of the desired inversion domain (includes region of interest and buffer elements).
   * - ``LatMinCustomStateVector``
     - Minimum latitude edge of the desired inversion domain (includes region of interest and buffer elements).
   * - ``LatMaxCustomStateVector``
     - Maximum latitude edge of the desired inversion domain (includes region of interest and buffer elements).

Inversion
~~~~~~~~~
.. list-table::
   :widths: 40, 60
   :class: tight-table

   * - ``PriorError``
     - Error in the prior estimates (1-sigma; relative). E.g., ``0.5`` means 50% error.
   * - ``ObsError``
     - Observational error (1-sigma; absolute; ppb). E.g., ``15`` means 15 ppb error.
   * - ``Gamma``
     - Regularization parameter; typically between 0 and 1.

Grid
~~~~
.. list-table::
   :widths: 40, 60
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
   :widths: 40, 60
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
   :widths: 40, 60
   :class: tight-table

   * - ``RunSetup``
     - Run the setup script, including selected setup modules above.
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
   :widths: 40, 60
   :class: tight-table

   * - ``DoPreview``
     - Boolean to run the IMI preview (``true``) or not (``false``).
   * - ``DOFSThreshold``
     - Threshold for estimated DOFS below which the IMI should automatically exit with a warning after performing the preview.
       Default value ``0`` prevents exit.

Compute Resources to Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
These settings are specific to Harvard's Cannon compute cluster. Not used for cloud runs.

.. list-table::
   :widths: 40, 60
   :class: tight-table

   * - ``nCPUs``
     - Number of cpus to use in ``sbatch`` scripts.
   * - ``partition``
     - Name of the cluster partition to use with ``sbatch`` (eg. ``"huce_cascade"``).
