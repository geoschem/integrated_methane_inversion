#!/bin/bash

# Run the IMI inversion stage for all target faces set in
# gchp_imi_common.sh, max_parallel_inversions faces at a time. This waits
# for all faces to finish, so run it on the head node inside tmux:
#   tmux new -s inversion
#   ./run_gchp_imi_inversion.sh
#
# Faces whose inversion_status ends with DONE are skipped, so this can be
# rerun after fixing failed faces.

set -e
cd "$(dirname "$0")"
source ./gchp_imi_common.sh

mkdir -p "${test_dir}/logs"
echo "Running inversion stage for faces ${start_face} to ${num_faces}, ${max_parallel_inversions} at a time"

for ((i_face=start_face; i_face<=num_faces; i_face++)); do
    face_name=$(face_name_for ${i_face})
    status_file="${test_dir}/${face_name}/inversion_status"
    if [[ -f "${status_file}" && "$(tail -n 1 "${status_file}")" == DONE* ]]; then
        echo "  ${face_name}: already DONE, skipping"
        continue
    fi

    # Wait until fewer than max_parallel_inversions faces are running
    while (( $(jobs -rp | wc -l) >= max_parallel_inversions )); do
        wait -n || true
    done

    echo "  ${face_name}: started at $(date)"
    bash inversion_face.sh "${face_name}" > "${test_dir}/logs/inversion_${face_name}.log" 2>&1 &
done

# Wait for the last faces to finish
wait
echo ""
echo "Done. Status of each face:"
for ((i_face=start_face; i_face<=num_faces; i_face++)); do
    face_name=$(face_name_for ${i_face})
    echo "  ${face_name}: $(tail -n 1 "${test_dir}/${face_name}/inversion_status" 2>/dev/null)"
done
