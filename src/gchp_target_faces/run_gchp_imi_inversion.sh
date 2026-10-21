#!/bin/bash

# Run satellite diagnostics and the IMI inversion stage for all target
# faces set in gchp_imi_common.sh, for one segment, max_parallel_inversions
# faces at a time. See inversion_face.sh for what is done for each face.
#
# Run this after all Jacobian runs have reached the segment end
# (check_gchp_imi_jacobian.sh shows all runs DONE). It waits for all faces
# to finish, so run it on the head node inside tmux:
#   tmux new -s inversion
#   ./run_gchp_imi_inversion.sh <segment_end YYYYMMDD>
#
# Faces already DONE for this segment end are skipped, so this can be
# rerun after fixing failed faces.

set -e
cd "$(dirname "$0")"
source ./gchp_imi_common.sh

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <segment_end YYYYMMDD>"
    exit 1
fi
segment_end=$1
if [[ ! "${segment_end}" =~ ^[0-9]{8}$ ]] || ! date -d "${segment_end}" >/dev/null 2>&1; then
    echo "ERROR: Segment end '${segment_end}' is not a valid date YYYYMMDD"
    exit 1
fi

#---------------------------------------------------------------------------
# Function to print the last line of a face's inversion_status, or nothing
#---------------------------------------------------------------------------
# Argument: face name
last_status() {
    tail -n 1 "${test_dir}/$1/inversion_status" 2>/dev/null || true
}

mkdir -p "${test_dir}/logs"
echo "Processing faces ${start_face} to ${num_faces} through ${segment_end}, ${max_parallel_inversions} at a time"

#-----------------------------------
# Loop over target faces
#-----------------------------------
for ((i_face=start_face; i_face<=num_faces; i_face++)); do
    face_name=$(face_name_for ${i_face})
    if [[ "$(last_status "${face_name}")" == "DONE ${segment_end} "* ]]; then
        echo "  ${face_name}: already DONE for ${segment_end}, skipping"
        continue
    fi

    # Wait until fewer than max_parallel_inversions faces are running
    while (( $(jobs -rp | wc -l) >= max_parallel_inversions )); do
        wait -n || true
    done

    echo "  ${face_name}: started at $(date)"
    bash inversion_face.sh "${face_name}" "${segment_end}" \
        > "${test_dir}/logs/inversion_${face_name}.${segment_end}.log" 2>&1 &
done

#-----------------------------------
# Wait for the last faces to finish
#-----------------------------------
wait
echo ""
echo "Done. Status of each face:"
for ((i_face=start_face; i_face<=num_faces; i_face++)); do
    face_name=$(face_name_for ${i_face})
    echo "  ${face_name}: $(last_status "${face_name}")"
done
