IMI input data 
==============

Input data for the IMI is stored on AWS and automatically accessed when running the IMI on AWS. When
running the IMI on a local compute cluster, you will need to download the required input data.

Input data for GEOS-Chem
------------------------

The forward model in the IMI is `GEOS-Chem chemical transport model <https://geoschem.github.io/>`_. 
The input files needed to run GEOS-Chem within the IMI can be accessed from the  `GEOS-Chem Nested Input Data` 
<https://registry.opendata.aws/geoschem-nested-input-data/>`_ portal (aka `s3://gcgrid <https://gcgrid.s3.amazonaws.com/index.html>`_).

The GEOS-Chem Input Data portal is part of the `AWS Open Data Sponsorship Program <https://aws.amazon.com/opendata/open-data-sponsorship-program/>`_.
As a result, **the data is completely free to use**.  You will NOT incur any data egress fees when downloading data from the
`s3://gcgrid <https://gcgrid.s3.amazonaws.com/index.html>`_ bucket.  

Input data for GEOS-Chem include:

- `Meteorological fields <https://geos-chem.readthedocs.io/en/stable/geos-chem-shared-docs/doc/gcid-data-on-aws.html#gcid-data-org-met>`_ (GEOS-FP or MERRA-2)
- `Emissions inventories <https://geos-chem.readthedocs.io/en/latest/geos-chem-shared-docs/doc/gcid-data-on-aws.html#gcid-data-org-emis-inputs>`_
- Chemistry input data (e.g. archived OH fields)
- `Initial conditions for starting GEOS-Chem simulations <https://geos-chem.readthedocs.io/en/latest/geos-chem-shared-docs/doc/gcid-data-on-aws.html#gcid-data-org-init-cond>`_

To automatically download these data for your inversion domain and time period, we recommend setting the ``*DryRun`` options 
in the :doc:`IMI configuration file <imi-config-file>` to true. This will execute a `GEOS-chem dry-run simulation`<https://geos-chem.readthedocs.io/en/latest/gcclassic-user-guide/dry-run.html>_
to identify and download the necessary giles.

You may also download these files manually using `AWS CLI <https://aws.amazon.com/cli/`_. See `this tutorial <https://geos-chem.readthedocs.io/en/latest/geos-chem-shared-docs/doc/gcid-awscli-tutorial.html>`_ for instructions.
Alternatively, you can access the data via `AWS S3 Explorer <https://s3.amazonaws.com/gcgrid/index.html>`_.

Satellite data
--------------

The IMI currently supports the following datasets for use in the inversion:

- ``TROPOMI``: The operational `TROPOMI <https://www.tropomi.eu/>`_ retrieval product developed by the SRON Netherlands Institute for Space Research.
- ``blendedTROPOMI``: The Blended TROPOMI+GOSAT retrieval product developed by `Balasus et al. (2023) <https://amt.copernicus.org/articles/16/3787/2023/>`_ to mitigate retrieval artifacts in the operational product.


Boundary conditions
-------------------

The IMI uses gridded 3D boundary conditions saved out from global 2°x2.5° GEOS-Chem simulation to define the inflow/outflow of methane at the edges of a regional simulation domain, essential for accurate regional inversions.
These files are further smoothed using either TROPOMI or the blended TROPOMI+GOSAT retrieval products to remove any systematic biases. These boundary condition files archive hourly species concentrations
and may also be used as initial conditions for the forward model simulations.

The IMI boundary condition files are stored on AWS in `s3://imi-boundary-conditions <https://imi-boundary-conditions.s3.amazonaws.com/index.html>`_ and may be downloaded using AWS CLI as described above.


