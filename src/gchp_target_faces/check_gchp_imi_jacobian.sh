#!/bin/bash

# Use this bash script to submit N number of Jacobian run slurm jobs for a
# target end date. If end date is reached for all Jacobian run directories then
# no more jobs are submitted. To control how many jobs are run at once, set up
# the pcluster partition (queue) to specify max number of jobs.

#------------------------------------------------------------------
# Assign arguments
#-----------------------------------------------------------------

if [ "$#" -ne 3 ]; then
    echo -e "\nERROR: Exactly 3 arguments are required. You passed $#.\n"
    echo "Usage: $0 <YYYYMMDD> <max_duration_days> <test_name>"
    exit 1
fi

end_date=$1
max_duration=$2
test_name=$3  # e.g. Instance_tests_400tracer_2025_C36S10

#-----------------------------------------
# Set Directories
#-----------------------------------------

# Must be set up prior to running this script, using setup_gchp_imi_jacobian.sh
main_dir=/home/ec2-user/fsx/${test_name}

# Expected directory name prefix
face_dir_prefix=Global_1yr_2025_C36S10

#------------------------------------------------------
# Do argument checks
#------------------------------------------------------

# Check that test directory exists
if [ ! -d "$main_dir" ]; then
    echo -e "\nERROR: Test directory $main_dir does not exist\n"
elif [[ ! "${end_date}" =~ ^[0-9]{8}$ ]]; then
    echo -e "\nERROR: First argument '${end_date}' must be an 8-digit integer (YYYYMMDD).\n" >&2
    exit 1
elif ! date -d "${end_date}" >/dev/null 2>&1; then
    echo -e "\nERROR: First argument '${end_date}' is not a valid calendar date.\n" >&2
    exit 1
elif [[ ! "${max_duration}" =~ ^-?[0-9]+$ ]]; then
    echo -e "\nERROR: Second argument '${max_duration}' is not a valid integer.\n" >&2
    exit 1
fi

#-----------------------------------------------
# Do individual run directory checks in loop
#-----------------------------------------------

# Initialize stats to track
num_active_runs=0
num_runs_ready_to_submit=0
num_runs_already_done=0
num_runs=0

# Option to check only the first N runs
N=957

# Loop over target faces
echo -e "\nChecking status of all jacobian run directories"
cd ${main_dir}
for face_dir in ${face_dir_prefix}_*; do
    face_id="${face_dir##*_}"
    cd ${main_dir}/${face_dir}/jacobian_runs

    # Loop over jacobian runs for this target face
    for jac_dir in ${face_dir_prefix}_*; do

	# Break out if reached N runs (configured before loop)
	((num_runs++))
	if (( num_runs > N )); then
            echo "Limit of $N runs reached. Exiting run directory loop."
            break 2
	fi
	
	run_id="${jac_dir##*_}"
	cd ${main_dir}/$face_dir/jacobian_runs/${jac_dir}
	echo "   ID ${face_id}_${run_id}: "

        #--------------------------------------------------------------------------------
        # If latest restart file date >= end date, increment num_runs_already_done
	#--------------------------------------------------------------------------------
	latest_file=$(ls -1 ./Restarts/GEOSChem.Restart.*_0000z.c36.nc4 2>/dev/null | tail -n 1)
	if [ -n "$latest_file" ]; then
	    filename=$(basename "$latest_file")
	    file_date="${filename:17:8}"
	    if [ "${file_date}" -ge "${end_date}" ]; then
		echo "        Run not needed (latest restart file is $filename)"
		((num_runs_already_done++))
		cd ..
		continue
	    fi
	fi

	#--------------------------------------------------------------------------------
        # If run is in progress, increment num_active_runs
	#--------------------------------------------------------------------------------
	job_running=false
        if [ -f "job_status" ]; then
	    last_line=$(tail -n 1 job_status)
	    if [[ "${last_line}" == *"Submitted batch job"* ]]; then
		echo "        Job currently in progress"
		((num_active_runs++))
		cd ..
		continue
	    fi
        fi
	
	#------------------------------------------------------------------------------------
	# Additional checks to determine if run is ready to submit
	#------------------------------------------------------------------------------------

	# Exit with error if cap_restart does not exist
	if [ ! -f cap_restart ]; then
	    echo -e "        ERROR: File cap_restart does not exist for ${face_id}_${run_id}. Investigate problem directory.\n"
	    exit 1
	fi

	start_date=$(cut -c 1-8 cap_restart)
	start_sec=$(date -d "${start_date}" +%s)
	end_sec=$(date -d "${end_date}" +%s)

	# Exit with error if start date is after end date
	if [ "${start_sec}" -gt "${end_sec}" ]; then
	    echo -e "        ERROR: Start date (${start_date}) cannot be after end date (${end_date}). Edit end date passed to script.\n"
	    exit 1
	fi

	# If start is 0 to max duration days before end, increment num_runs_ready_to_submit
	diff_days=$(( (end_sec - start_sec) / 86400 ))
	if [ "$diff_days" -le "${max_duration}" ]; then
	    echo "        Run can be submitted. Will start on ${start_date} and run for ${diff_days} days"
	    ((num_runs_ready_to_submit++))
	    cd ..
	    continue
	else
	    # Exit with error if submitted run would exceed max duration
	    echo -e "        ERROR: Run length (${diff_days} days) would exceed max allowed (${max_duration} days). Edit arguments passed to script.\n"
	    exit 1
	fi

    done
done

#------------------------------------------------------------------------------------
# Print summary
#------------------------------------------------------------------------------------
echo -e "\nSummary: "
echo "   End date: ${end_date}"
echo "   Total runs already run through end date: ${num_runs_already_done}"
echo "   Total runs not yet there but are actively running: ${num_active_runs}"
echo "   Total runs ready to submit: ${num_runs_ready_to_submit}"
echo -e "\n   Sanity check total runs counted: $((num_runs_already_done + num_active_runs + num_runs_ready_to_submit))\n"

exit 0
