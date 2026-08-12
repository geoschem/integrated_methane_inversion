May 2026

The jacobian integration tests use external test data that is not checked into
this repository. Before running the tests, set `JACOBIAN_TEST_DATA_ROOT` to
either:

- the directory containing `baseline_jacobian_inversion3`, or
- the `baseline_jacobian_inversion3` directory itself.

For example:

`export JACOBIAN_TEST_DATA_ROOT=/path/to/test_data`

or:

`export JACOBIAN_TEST_DATA_ROOT=/path/to/test_data/baseline_jacobian_inversion3`

Then run:

`pytest src/inversion_scripts/tests/test_jacobian.py`

The test copies the data into pytest's temporary directory and rewrites the
copied config `OutputPath` so test runs do not write generated files into the
repository.

The outputs in `baseline_jacobian_inversion3/inversion/data_converted` and
`data_visualization` were generated with this command, replacing `$TEST_DATA`
with the external directory containing `baseline_jacobian_inversion3`:

`python -u src/inversion_scripts/jacobian.py $TEST_DATA/baseline_jacobian_inversion3/inversion/ $TEST_DATA/baseline_jacobian_inversion3/config_baseline_jacobian_inversion3.yml "20180505" "20180507" "-104" "-103" "31" "32" "27" "CH4" $TEST_DATA/baseline_jacobian_inversion3/satellite_data_may/ "BlendedTROPOMI" "false" "false" "1" "true" "false"`

The same command is run in `test_jacobian.py::test_jacobian_end_to_end_baseline`
to verify that the actual outputs of `jacobian.py` are the same as the expected
outputs from a fixed run. If you make any changes to the observation operator
or the jacobian construction that change jacobian entries by more than O(1 ppb),
the `test_jacobian.py` tests will fail. Once you've verified that your changes
are correct, re-run the command above to regenerate the expected output files.
