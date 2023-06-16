# README: Environment files for the Harvard Cannon cluster

 ## Environment files for the IMI

 In this directory we provide sample **environment files**.  These are bash scripts that will load all of the necessary software for running the IMI into your login environment.


 | Environment file                    | With compilers            | Software included | 
 | ----------------                    | --------------            | ----------------- |
 | `gcclassic.rocky+gnu10.env`         | gcc, g++, gfortran 10.2.0 | FASRC + Spack     |
 | `gcclassic.gocky+gnu10.minimal.env` | gcc, g++, gfortran 10.2.0 | FASRC             |
 | `gcclassic.rocky+gnu12.env`         | gcc, g++, gfortran 12.2.0 | FASRC + Spack     |
 | `gcclassic.gocky+gnu12.minimal.env` | gcc, g++, gfortran 12.2.0 | FASRC             |

 NOTES:

 1. Use the `gcclassic.rocky+gnu10.env` or `gcclassic.rocky+gnu12.env` environment files in interactive login sessions. These environment files will load the complete set of software packages into your Linux environment.

 2. Use the `gcclassic.rocky+gnu10.minimal.env` or `gcclassic.rocky+gnu12.minimal.env` environment files in run scripts that will be submitted to Cannon queues. These environment files will only load the essential software packages needed to run GEOS-Chem Classic for the IMI, without any of the interactive utilities that have been built with Spack.

 3. We recommend that you copy or link these files to a convenient location (such as your home directory).

 To apply the settings contained in an environment file, type

 ```console
 source <env-file-name>
 ```
 where `<env-file-name>` refers to one of the files listed above.

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
