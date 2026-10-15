#!/bin/bash

# Print the status of every Jacobian run directory for the target faces set
# in gchp_imi_common.sh, then a summary. Usage:
#   ./check_gchp_imi_jacobian.sh <end_date YYYYMMDD> <max_days>
#
# Status of each run (see run_status in gchp_imi_common.sh):
#   DONE            cap_restart has reached the end date
#   ACTIVE          a job for the run is pending or running in Slurm
#   READY           can be submitted
#   FAILED          the last job wrote FAILED to job_status. See the run log.
#   NO_CAP_RESTART  cap_restart is missing. Investigate the run directory.
#
# Exits with an error if a READY run would need more than max_days to reach
# the end date, since that usually means a typo in the arguments.

set -e
cd "$(dirname "$0")"
source ./gchp_imi_common.sh

#------------------------------------------------------
# Check arguments
#------------------------------------------------------
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <end_date YYYYMMDD> <max_days>"
    exit 1
fi
end_date=$1
max_days=$2

if [[ ! -d "${test_dir}" ]]; then
    echo "ERROR: Test directory ${test_dir} does not exist"
    exit 1
elif [[ ! "${end_date}" =~ ^[0-9]{8}$ ]] || ! date -d "${end_date}" >/dev/null 2>&1; then
    echo "ERROR: End date '${end_date}' is not a valid date YYYYMMDD"
    exit 1
elif [[ ! "${max_days}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Maximum days '${max_days}' is not a valid integer"
    exit 1
fi

#------------------------------------------------------
# Loop over target faces and their Jacobian runs
#------------------------------------------------------
active=$(active_run_dirs)
num_done=0; num_active=0; num_ready=0; num_failed=0; num_other=0; num_too_long=0

echo "Status of Jacobian runs in ${test_dir}, end date ${end_date}:"
for ((i_face=start_face; i_face<=num_faces; i_face++)); do
    face_name=$(face_name_for ${i_face})
    for run_dir in "${test_dir}/${face_name}/jacobian_runs/${face_name}_"*; do
        [[ -d "${run_dir}" ]] || continue
        status=$(run_status "${run_dir}" "${end_date}" "${active}")
        details=""
        case "${status}" in
            DONE)   num_done=$((num_done + 1)) ;;
            ACTIVE) num_active=$((num_active + 1)) ;;
            READY)
                num_ready=$((num_ready + 1))
                start_date=$(cut -c 1-8 "${run_dir}/cap_restart")
                days=$(( ( $(date -d "${end_date}" +%s) - $(date -d "${start_date}" +%s) ) / 86400 ))
                details="starts ${start_date}, ${days} days"
                if (( days > max_days )); then
                    details="${details}, MORE THAN ${max_days} DAYS"
                    num_too_long=$((num_too_long + 1))
                fi
                ;;
            FAILED)
                num_failed=$((num_failed + 1))
                details=$(grep "^FAILED" "${run_dir}/job_status" | tail -n 1)
                ;;
            *)      num_other=$((num_other + 1)) ;;
        esac
        printf "  %-36s %-15s %s\n" "$(basename "${run_dir}")" "${status}" "${details}"
    done
done

#------------------------------------------------------
# Print summary
#------------------------------------------------------
echo ""
echo "Summary for end date ${end_date}:"
echo "  DONE:   ${num_done}"
echo "  ACTIVE: ${num_active}"
echo "  READY:  ${num_ready}"
echo "  FAILED: ${num_failed}"
echo "  Other (e.g. NO_CAP_RESTART): ${num_other}"
echo "  Total:  $((num_done + num_active + num_ready + num_failed + num_other))"

if (( num_too_long > 0 )); then
    echo ""
    echo "ERROR: ${num_too_long} READY runs need more than ${max_days} days. Check the arguments."
    exit 1
fi
