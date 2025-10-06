# README: Environment files for the MSU-HPC Orion

 ## Environment files for the IMI

 In this directory we provide a sample **environment file**.  This bash script will load all of the necessary software for running the IMI into your login environment.

   $ source gcclassic.intel2022.orion.env

This will load the following modules:

  1) intel-oneapi-compilers/2022.2.1   4) libszip/2.1.1   7) hdf5/1.12.2
  2) intel-oneapi-mpi/2021.7.1         5) zlib/1.2.13     8) netcdf-c/4.9.0
  3) intel-oneapi-mkl/2022.2.1         6) hdf4/4.2.16     9) netcdf-fortran/4.6.0


## Python environment for the IMI

 Use the file `imi_env.yml` to create a Conda environment for the IMI on your system.  This Conda environment will load all of the packages needed to execute the scripts in the repository.

 Instructions for creating a Conda environment from `imi_env.yml`

 ```console
 conda env create --file imi_env.yml
 ```

 To activate:

 ```console
 conda activate imi_env
 ```

 To deactivate:

 ```console
 conda deactivate
 ```

 For more information, see https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#


## Configuration file for the IMI

 File config.msu-hpc_orion.global_inv.yml is the IMI configuration file for a global inversion. A regional inversion file for the IMI on MSU-HPC Orion is not yet developed. You may create one by adapting envs/Harvard-Cannon/config.harvard-cannon.yml, using config.msu-hpc_orion.global_inv.yml as a guide for data paths.


IMPORTANT NOTES:

1. MSU-HPC Orion compute nodes do not have internet access. This causes a problem if running preview in the IMI. We therefore recommend finding a work-around to using compute nodes when running preview.

2. Some data needed for the IMI may not be available on MSU-HPC Orion. You will need to download required data. See config.msu-hpc_orion.global_inv.yml for the location of input data on the cluster.

3. If issues are encountered using the IMI you may open a GitHub issue at https://github.com/geoschem/integrated_methane_inversion.