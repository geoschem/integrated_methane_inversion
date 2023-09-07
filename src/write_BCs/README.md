## Directions to generate the boundary conditions.

1. Change settings in `config_boundary_conditions.yml`.
   - You need to provide a restart file for 15 days before the beginning of the day you want BCs for (unless 1 April 2018 is your start). In that case, the file `GEOSChem.Restart.20180401_0000z.nc4` is used, which comes from a simulation constrained by NOAA surface observations.
   - This is okay though because daily restarts are output by the GEOS-Chem run.
2. Run `sbatch -p huce_cascade run_boundary_conditions.sh`.
   - GEOS-Chem will be run first (2.0 x 2.5, GEOS-FP, CH4, 47 L, daily restart files).
   - Bias-corrected boundary conditions will be written via `write_boundary_conditions.py` (description in file).