"""Unit and integration tests for jacobian.py."""

import pickle
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


class TestEndToEndIntegration:
    """End-to-end integration tests running jacobian.py with real test data."""

    @pytest.fixture
    def baseline_test_data_dir(self):
        """Path to baseline test data directory."""
        test_data_dir = Path(__file__).parent / "test_data" / "baseline_jacobian_inversion3"
        if not test_data_dir.exists():
            pytest.skip(f"Test data not found at {test_data_dir}")
        return test_data_dir

    @pytest.fixture
    def temp_inversion_workspace(self, temp_dir, baseline_test_data_dir):
        """
        Create a temporary inversion workspace by copying baseline test data.
        Clones baseline_jacobian_inversion3 structure into temp directory.
        """

        # Copy the entire baseline directory to temp location (note this includes golden files)
        temp_workspace = temp_dir / "inversion_workspace"
        temp_workspace.mkdir(exist_ok=True)

        # Copy inversion subdirectory and all its contents
        src_inversion_dir = baseline_test_data_dir / "inversion"
        dest_inversion_dir = temp_workspace / "inversion"

        if src_inversion_dir.exists():
            shutil.copytree(src_inversion_dir, dest_inversion_dir)
            # NOTE: this line removes the data_converted directory that 
            # contains the golden output files, so that the test can generate 
            # new output files and compare against the golden files without 
            # accidentally using the pre-existing files in the baseline directory.
            shutil.rmtree(dest_inversion_dir / "data_converted")
            shutil.rmtree(dest_inversion_dir / "data_visualization")

        # Copy satellite data
        src_satellite_dir = baseline_test_data_dir / "satellite_data_may"
        dest_satellite_dir = temp_workspace / "satellite_data_may"

        if src_satellite_dir.exists():
            shutil.copytree(src_satellite_dir, dest_satellite_dir)

        # Create data directories that will be populated by jacobian.py
        (dest_inversion_dir / "data_geoschem").mkdir(exist_ok=True)
        (dest_inversion_dir / "data_converted").mkdir(exist_ok=True)
        (dest_inversion_dir / "data_visualization").mkdir(exist_ok=True)

        # Copy GEOS-Chem data if it exists in baseline
        gc_src = baseline_test_data_dir / "data_geoschem"
        if gc_src.exists():
            gc_dest = dest_inversion_dir / "data_geoschem"
            if not gc_dest.exists():
                shutil.copytree(gc_src, gc_dest)

        return temp_workspace

    @pytest.fixture
    def avg_tropomi_expected_outputs_dir(self, baseline_test_data_dir):
        """Path to expected output files of average tropomi operator."""
        return baseline_test_data_dir / "inversion" / "data_converted"

    @pytest.fixture
    def regular_tropomi_expected_outputs_dir(self, baseline_test_data_dir):
        """Path to expected output files of regular tropomi operator."""
        return baseline_test_data_dir / "inversion" / "data_visualization"

    def test_jacobian_end_to_end_baseline(
        self,
        temp_inversion_workspace,
        baseline_test_data_dir,
        avg_tropomi_expected_outputs_dir,
        regular_tropomi_expected_outputs_dir,
    ):
        """
        End-to-end test: Run jacobian.py with baseline test data and verify outputs.
        
        This test:
        1. Copies baseline test data to temp directory
        2. Calls jacobian.py with test parameters
        3. Compares generated pickle outputs with expected reference files
        4. Verifies output structure and numerical accuracy
        """

        workdir = temp_inversion_workspace / "inversion"
        config_file = baseline_test_data_dir / "config_baseline_jacobian_inversion3.yml"
        satellite_cache = baseline_test_data_dir / "satellite_data_may"

        # Test parameters (matching baseline_jacobian_inversion3 setup)
        startday = "20180505"
        endday = "20180507"
        lonmin = "-104"
        lonmax = "-103"
        latmin = "31"
        latmax = "32"
        n_elements = "27"
        species = "CH4"
        satellite_product = "BlendedTROPOMI"
        use_water_obs = "false"
        is_post = "false"
        period_i = "1"
        build_jacobian = "true"
        viz_prior = "false"

        # Build command
        cmd = [
            sys.executable,
            "-u",
            "src/inversion_scripts/jacobian.py",
            str(workdir),
            str(config_file),
            startday,
            endday,
            lonmin,
            lonmax,
            latmin,
            latmax,
            n_elements,
            species,
            str(satellite_cache),
            satellite_product,
            use_water_obs,
            is_post,
            period_i,
            build_jacobian,
            viz_prior,
        ]

        # Run jacobian.py
        env = {**os.environ}
        workspace_root = Path(__file__).parent.parent.parent.parent
        env["PYTHONPATH"] = f"{workspace_root}:{env.get('PYTHONPATH', '')}"

        result = subprocess.run(
            cmd,
            cwd=workspace_root,
            env=env,
            capture_output=True,
            text=True,
            timeout=300,
        )

        # Check that jacobian.py ran successfully
        if result.returncode != 0:
            pytest.fail(
                f"jacobian.py exited with code {result.returncode}\n"
                f"STDOUT:\n{result.stdout}\n"
                f"STDERR:\n{result.stderr}"
            )

        # Verify avg tropomi operator output pickle files were created
        output_dir = workdir / "data_converted"
        output_pkl_files = list(output_dir.glob("*_GCtoSatellite.pkl"))

        assert len(output_pkl_files) > 0, (
            f"No output pickle files found in {output_dir}. "
            f"Expected *_GCtoSatellite.pkl files"
        )

        # Compare outputs with expected files
        if avg_tropomi_expected_outputs_dir.exists():
            expected_pkl_files = list(avg_tropomi_expected_outputs_dir.glob("*_GCtoSatellite.pkl"))

            for expected_file in expected_pkl_files:
                output_file = output_dir / expected_file.name
                if output_file.exists():
                    # Load both pickle files
                    with open(expected_file, "rb") as f:
                        expected_output = pickle.load(f)

                    with open(output_file, "rb") as f:
                        actual_output = pickle.load(f)

                    # Verify structure matches
                    assert set(expected_output.keys()) == set(actual_output.keys()), (
                        f"Output keys don't match for {expected_file.name}. "
                        f"Expected: {set(expected_output.keys())}, "
                        f"Got: {set(actual_output.keys())}"
                    )

                    # Compare numerical arrays with tolerance
                    for key in expected_output.keys():
                        if isinstance(expected_output[key], np.ndarray):
                            assert isinstance(actual_output[key], np.ndarray), (
                                f"Array {key} has different type: "
                                f"expected ndarray, got {type(actual_output[key])}"
                            )

                            assert expected_output[key].shape == actual_output[key].shape, (
                                f"Array {key} has different shape: "
                                f"expected {expected_output[key].shape}, "
                                f"got {actual_output[key].shape}"
                            )

                            # Use relative tolerance for floating point comparisons
                            if expected_output[key].size > 0:
                                assert np.allclose(
                                    expected_output[key],
                                    actual_output[key],
                                    rtol=1e-5,
                                    atol=1e-8,
                                    equal_nan=True,
                                ), (
                                    f"Array {key} values don't match within tolerance. "
                                    f"Max difference: "
                                    f"{np.max(np.abs(expected_output[key] - actual_output[key]))}"
                                )

        # Verify regular tropomi operator output pickle files were created
        output_dir = workdir / "data_visualization"
        output_pkl_files = list(output_dir.glob("*_GCtoSatellite.pkl"))

        assert len(output_pkl_files) > 0, (
            f"No output pickle files found in {output_dir}. "
            f"Expected *_GCtoSatellite.pkl files"
        )

        # Compare outputs with expected files
        if regular_tropomi_expected_outputs_dir.exists():
            expected_pkl_files = list(regular_tropomi_expected_outputs_dir.glob("*_GCtoSatellite.pkl"))

            for expected_file in expected_pkl_files:
                output_file = output_dir / expected_file.name
                if output_file.exists():
                    # Load both pickle files
                    with open(expected_file, "rb") as f:
                        expected_output = pickle.load(f)

                    with open(output_file, "rb") as f:
                        actual_output = pickle.load(f)

                    # Verify structure matches
                    assert set(expected_output.keys()) == set(actual_output.keys()), (
                        f"Output keys don't match for {expected_file.name}. "
                        f"Expected: {set(expected_output.keys())}, "
                        f"Got: {set(actual_output.keys())}"
                    )

                    # Compare numerical arrays with tolerance
                    for key in expected_output.keys():
                        if isinstance(expected_output[key], np.ndarray):
                            assert isinstance(actual_output[key], np.ndarray), (
                                f"Array {key} has different type: "
                                f"expected ndarray, got {type(actual_output[key])}"
                            )

                            assert expected_output[key].shape == actual_output[key].shape, (
                                f"Array {key} has different shape: "
                                f"expected {expected_output[key].shape}, "
                                f"got {actual_output[key].shape}"
                            )

                            # Use relative tolerance for floating point comparisons
                            if expected_output[key].size > 0:
                                assert np.allclose(
                                    expected_output[key],
                                    actual_output[key],
                                    rtol=1e-5,
                                    atol=1e-8,
                                    equal_nan=True,
                                ), (
                                    f"Array {key} values don't match within tolerance. "
                                    f"Max difference: "
                                    f"{np.max(np.abs(expected_output[key] - actual_output[key]))}"
                                )


    # NOTE: this test is a bit extraneous, it's literally just testing whether the .pkl outputs of jacobian.py have obs_GC and K as expected
    def test_jacobian_output_structure_validation(
        self, temp_inversion_workspace, baseline_test_data_dir
    ):
        """
        Verify that jacobian.py output has correct structure and data types.
        """

        workdir = temp_inversion_workspace / "inversion"
        config_file = baseline_test_data_dir / "config_baseline_jacobian_inversion3.yml"
        satellite_cache = temp_inversion_workspace / "satellite_data_may"

        cmd = [
            sys.executable,
            "-u",
            "src/inversion_scripts/jacobian.py",
            str(workdir),
            str(config_file),
            "20180505",
            "20180507",
            "-104",
            "-103",
            "31",
            "32",
            "27",
            "CH4",
            str(satellite_cache),
            "BlendedTROPOMI",
            "false",
            "false",
            "1",
            "true",
            "false",
        ]

        env = {**os.environ}
        workspace_root = Path(__file__).parent.parent.parent.parent
        env["PYTHONPATH"] = f"{workspace_root}:{env.get('PYTHONPATH', '')}"

        result = subprocess.run(
            cmd,
            cwd=workspace_root,
            env=env,
            capture_output=True,
            text=True,
            timeout=300,
        )

        assert result.returncode == 0, (
            f"jacobian.py failed: {result.stderr}"
        )

        output_dir = workdir / "data_converted"
        output_files = list(output_dir.glob("*_GCtoSatellite.pkl"))

        for pkl_file in output_files:
            with open(pkl_file, "rb") as f:
                output = pickle.load(f)

            required_keys = [
                "obs_GC",
                "K"
            ]
            for key in required_keys:
                assert key in output, f"Missing required key: {key}"
                assert isinstance(output[key], np.ndarray), (
                    f"Key {key} should be numpy array, got {type(output[key])}"
                )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
