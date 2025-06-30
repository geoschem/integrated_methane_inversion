## Directions to generate the boundary conditions.

1. Change settings in `config_boundary_conditions.yml`.
   - `startDate`       - first day you want boundary conditions for.
   - `endDate`         - last day you want boundary conditions for.
   - `workDir`         - where boundary conditions will be generated.
   - `tropomiDir`      - where TROPOMI files are located.
   - `blendedDir`      - where Blended TROPOMI+GOSAT files are located.
   - `CondaEnv`        - conda environment to use for the Python script.
   - `GEOSChemEnv`     - environment file for GEOS-Chem.
   - `Partition`       - which partition to run the jobs on.
   - `restartFilePath` - restart file for GEOS-Chem.
      - if your simulation starts on 1 April 2018, this won't be used (`GEOSChem.Restart.20180401_0000z.nc4` will).
         - this file comes from a CH4 simulation by Todd Mooring that is constrained by NOAA surface observations.
      - if your simulation starts on another date, it should use a restart file from the simulation that started on 1 April 2018.
      - this is accommodated by the fact that the simulation run here is setup to write daily restart files.
      - to determine what day your restart file should be for, subtract 30 days from `startDate`.
         - for example, if `startDate: 20230601`, your restart file should be `GEOSChem.Restart.20230502_0000z.nc4`.
         - this allows for a 15 day previous average (using 18 May 2023 to 1 June 2023) and the 1 day at the start where no first hour will be written by GEOS-Chem.
         - in rare cases of low data density, a 30 day previous average is used (which would use 3 May 2023 to 1 June 2023).
   - `debug`           - whether or not to delete the `debug.log`.
      - this includes information about your environment file and the build of GEOS-Chem.
      - all important information and errors are written to `boundary_conditions.log`.
2. Run `sbatch -p huce_ice run_boundary_conditions.sh`.
   - GEOS-Chem will be run first (2.0 x 2.5, GEOS-FP, CH4, 47 L, daily restart files).
   - Bias-corrected boundary conditions will be written via `write_boundary_conditions.py` (description in file).
   - In `workDir`, the folders `tropomi-boundary-conditions` and `blended-boundary-conditions` will be populated.

## Directions for doing this operationally at Harvard
1. Run from `20180401` until the last day you have both satellite data and met fields.
   - example: `startDate: 20180401`, `endDate: 20230531`.
2. **Before deleting your `workDir`**, in `workDir/gc_run/Restarts/`, copy the restart file from > 30 days before the next day you will need boundary conditions to a persistent storage location.
   - example: `GEOSChem.Restart.20230401_0000z.nc4`.
3. When the satellite data and met fields become available, generate boundary conditions up until your new end date, but start with a little overlap to check for consistency.
   - example: `startDate: 20230501`, `endDate: 20230630` (using `GEOSChem.Restart.20230401_0000z.nc4`), then check consistency between your new and previously generated boundary conditions for `20230501` until `20230531` before using the new `20230601` to `20230630` boundary conditions.
