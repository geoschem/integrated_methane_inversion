Running the IMI on a local cluster
==================================

The IMI is setup to run on AWS by default. However, if you have a
local cluster available to you, you may choose to run the IMI
there. This option requires some manual changes and is therefore only
recommended for advanced users.

You must first ensure you have the proper `hardware <https://geos-chem.readthedocs.io/en/latest/starting/hardware-requirements.html>`__ and `software <https://geos-chem.readthedocs.io/en/latest/starting/software-requirements.html>`__ requirements for running GEOS-Chem.

When logged onto your local cluster, navigate to the path where you want to download the IMI repository and type the following command:

.. code-block:: console

    $ git clone https://github.com/geoschem/integrated_methane_inversion.git

This will clone the IMI code into a local folder named ``integrated methane_inversion``.

.. tip:: If you wish, you can clone the IMI repository into a
	 different local folder by supplying the name of the folder at
	 the end of the :command:`git clone` command. For example:
         ::
            git clone https://github.com/geoschem/integrated_methane_inversion.git imi-v1.0

Next, download the GEOS-Chem source code and its submodules within the
IMI folder using these commands:

.. code-block:: console

    $ cd integrated_methane_inversion
    $ git clone https://github.com/geoschem/GCClassic.git
    $ cd GCClassic
    $ git submodule update --init --recursive

See `Downloading rhe GEOS-Chem source code <https://geos-chem.readthedocs.io/en/latest/building-gc/download-source-code.html>`__ for more details.

Navigate back to the top-level IMI folder and view the contents:

.. code-block:: console

    $ cd ..
    $ ls 

Within the IMI is a subfolder called ``envs`` that constains files for
running the IMI on different systems. By default, files are provided
for AWS and Harvard's Cannon cluster. We recommend you add a subfolder
within ``envs`` for your own system to easily access files needed for
the IMI. In this directory, we recommend storing any environment files
needed to load the libraries for GEOS-Chem (e.g. fortran compiler, netcdf,
openmpi, cmake), a conda environment file, and an IMI configuration
file modified for your system. See the files in
``envs/Harvard-Cannon`` for examples.

Within the copied IMI configuration file, you will need to modify the
settings in the section labeled "Settings for running on your local
cluster." If you already have the GEOS-Chem input data on your system,
you may set the ``*DryRun`` options to ``false``.

It is recommended that you set up and run the IMI in stages when
running on a local cluster to ensure that each stage works
properly. You can do this by modifying the settings under ``Setup
modules`` and ``Run modules`` and turning them on one or a few at a
time. You may find that you need to manually edit some files. For
example, after creating the template run directory you should open
``ch4_run.template`` in a text editor and modify as needed for your
system (by default this script is set up to submit to a SLURM scheduler).





    
