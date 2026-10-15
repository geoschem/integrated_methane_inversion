#!/bin/bash

# This script looks at a given Restarts directory, determines first, last,
# and second to last dates of the files from the filenames, and deletes
# everything except files first, last, second to last, or Jan 1 (0101_) in
# the filename.

# Keeping this as loop over directories simply to make it more versatile in
# case want to do recursively from higher level

base_path=.
for dir in "$base_path"/Restarts; do
    if [ -d "$dir" ]; then

        echo -e "\nDirectory $dir\n"
	
	# Get all restart files, with either prefix
	files=($(ls $dir/GEOSChem.Restart.*.nc4 $dir/gcchem_internal_checkpoint.*.nc4 2>/dev/null))

	# Exit if not more than 3
	if [[ ${#files[@]} -lt 3 ]]; then
	    echo "Fewer than 3 files. No need to delete."
	    exit 0
	fi

	# Extract all unique YYYYMMDD dates across all files and sort them
	dates=($(for f in "${files[@]}"; do
	    echo "$f" | grep -oE '[0-9]{8}'
	done | sort -u))

	first_date="${dates[0]}"
        second_to_last_date="${dates[-2]}"
	last_date="${dates[-1]}"

	# Keep any file (either prefix) matching first, last, second to last, or
	# contains Jan 1
	for file in "${files[@]}"; do
	    if [[ "$file" == *"${first_date}"* || "$file" == *"${last_date}"* || "$file" == *"${second_to_last_date}"* || "$file" == *"0101_"* ]]; then
		echo "--> Keeping  $file"
	    else
		echo "--> Deleting $file"
		rm "$file"
	    fi
	done
    fi
done

exit 0
