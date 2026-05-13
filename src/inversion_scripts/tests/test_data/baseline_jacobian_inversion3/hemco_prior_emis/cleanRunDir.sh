#!/bin/bash

#============================================================================
# cleanRunDir.sh: Removes files created by HEMCO standalone from a run dir.
#
# Usage:
# ------
# $ ./cleanRunDir.sh          # Removes model output files in the run
#                             # directory.  Also prompts the user before
#                             # removing diagnostic output files from
#                             # from OutputDir/.
#
# $ ./cleanRunDir.sh --force  # Removes model output files in the run
#                             # directory, but will remove diagnostic
#                             # output files without prompting first.
#                             # USE WITH CAUTION!
#============================================================================

# Clean model output files in the run directory
rm -fv *~
rm -fv HEMCO.log
rm -fv log*
rm -fv slurm-*
rm -fv core.*
rm -fv fort.*

#----------------------------------------------------------------------------
# Clean data files in OutputDir.
# These are netCDF files (*.nc) and KPP standalone interface files (*.txt).
#----------------------------------------------------------------------------
if [[ "x${1}" == "x" ]]; then      # User confirmation required
    rm -Iv ./OutputDir/*.nc*
    rm -Iv ./OutputDir/*.txt
else                               # User Confirmation not required
    rm -fv ./OutputDir/*.nc*
    rm -fv ./OutputDir/*.txt*
fi
