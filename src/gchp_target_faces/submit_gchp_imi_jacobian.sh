#!/bin/bash

# Submit Jacobian runs that have not reached the end date and have no job in
# the Slurm queue. Runs check_gchp_imi_jacobian.sh first, then asks whether
# to continue and how many runs to submit. Usage:
#   ./submit_gchp_imi_jacobian.sh <end_date YYYYMMDD> <max_days>
#
# Each run is submitted from its own directory as
#   sbatch -p ${run_partition} -t ${run_time} -J <run dir name> gchp_ch4_run_spot.run <end_date> <max_days>
# To limit how many jobs run at once, set the maximum number of nodes for
# the partition in the cluster config.

set -e
cd "$(dirname "$0")"
source ./gchp_imi_common.sh

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <end_date YYYYMMDD> <max_days>"
    exit 1
fi
end_date=$1
max_days=$2

#------------------------------------------------------
# Check status first. This stops if the arguments are bad.
#------------------------------------------------------
./check_gchp_imi_jacobian.sh "${end_date}" "${max_days}"

echo ""
read -p "Proceed to submitting runs? [y/N] " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Exiting."
    exit 0
fi
read -p "Also resubmit FAILED runs? [y/N] " -n 1 -r
echo
retry_failed=false
if [[ $REPLY =~ ^[Yy]$ ]]; then
    retry_failed=true
fi
read -p "Maximum number of runs to submit: " max_runs
if [[ ! "${max_runs}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: '${max_runs}' is not a valid integer."
    exit 1
fi

log_file="${test_dir}/logs/submit.log"
echo -e "\n$(date): end date ${end_date}, max days ${max_days}, max runs ${max_runs}, retry failed ${retry_failed}" >> "${log_file}"

#------------------------------------------------------
# Loop over target faces and their Jacobian runs and submit
#------------------------------------------------------

# Nodes and cores for each run, from the base IMI config
run_nodes=$(grep "^NUM_NODES:" "${template_config}" | awk '{print $2}')
run_cores=$(grep "^TOTAL_CORES:" "${template_config}" | awk '{print $2}')

# Read the queue again in case it changed while answering the questions
active=$(active_run_dirs)
num_submitted=0
for ((i_face=start_face; i_face<=num_faces; i_face++)); do
    face_name=$(face_name_for ${i_face})
    for run_dir in "${test_dir}/${face_name}/jacobian_runs/${face_name}_"*; do
        [[ -d "${run_dir}" ]] || continue
        if (( num_submitted >= max_runs )); then
            echo "Reached the limit of ${max_runs} submitted runs."
            break 2
        fi

        status=$(run_status "${run_dir}" "${end_date}" "${active}")
        if [[ "${status}" == READY || ( "${status}" == FAILED && "${retry_failed}" == true ) ]]; then
            run_name=$(basename "${run_dir}")
	    msg=$(cd "${run_dir}" && sbatch -p "${run_partition}" \
					    -t "${run_time}" \
					    -N "${run_nodes}" \
					    -n "${run_cores}" \
					    -J "${run_name}" "${run_name}.run" "${end_date}" "${max_days}")
            echo "  ${run_name} (${status}): ${msg}" | tee -a "${log_file}"
            num_submitted=$((num_submitted + 1))
        fi
    done
done

echo "Done. ${num_submitted} jobs submitted." | tee -a "${log_file}"
