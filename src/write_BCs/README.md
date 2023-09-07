## Directions to generate the boundary conditions.

1. Change settings in `config_boundary_conditions.yml`.
2. Run `sbatch run_boundary_conditions.sh`.
   - GEOS-Chem will be run first (2.0 x 2.5, GEOS-FP, CH4, 47 L, daily restart files).
   - Bias-corrected boundary conditions will be written via `write_boundary_conditions.py` (description in file).