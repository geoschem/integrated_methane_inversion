import importlib.util
from pathlib import Path
import sys

import numpy as np
import netCDF4 as nc
import pytest
import xarray as xr

from src.inversion_scripts.operators.superobservation import (
    OPTIONAL_SUPEROBSERVATION_VARIABLES,
    structured_superobservations_to_dataset,
    validate_superobservation_dataset,
)
from src.inversion_scripts.operators.msat_funcs import (
    average_methanesat_observations,
    coordinate_edges,
)
from src.inversion_scripts.operators.satellite_operator import (
    _ravel_rectilinear_grid_indices,
)
from src.inversion_scripts.operators.superobservation import imi_superobservation_dtype
from GOOPy.parsers import read_MSAT as read_goopy_msat


def _observations():
    dtype = imi_superobservation_dtype(species="CH4", n_pressure_edges=3, n_layers=2)
    observations = np.zeros(1, dtype=dtype)
    observations["iGC"] = 2
    observations["jGC"] = 3
    observations["lat"] = 40.0
    observations["lon"] = -105.0
    observations["lat_sat"] = 40.01
    observations["lon_sat"] = -105.01
    observations["CH4"] = 1900.0
    observations["time"] = "20240101_12"
    observations["p_sat"] = [[1000.0, 500.0, 0.1]]
    observations["dry_air_subcolumns"] = [[2.0, 1.0]]
    observations["apriori"] = [[4.0e-6, 2.0e-6]]
    observations["avkern"] = [[0.5, 1.0]]
    observations["observation_count"] = 4.0
    return observations


def test_canonical_dataset_normalizes_units_and_vertical_fields():
    dataset = structured_superobservations_to_dataset(
        _observations(), "CH4", "input.nc"
    )

    assert dataset.attrs["format_name"] == "IMI_superobservation"
    np.testing.assert_allclose(dataset["column"], [1.9e-6])
    np.testing.assert_allclose(dataset["pressure_edges"], [[100000, 50000, 10]])
    np.testing.assert_allclose(dataset["pressure_weight"], [[2 / 3, 1 / 3]])
    np.testing.assert_allclose(dataset["prior_profile"], [[2e-6, 2e-6]])


@pytest.mark.parametrize("invalid_value", [0.0, -1.0, np.nan, np.inf])
def test_canonical_dataset_rejects_invalid_dry_air_subcolumns(invalid_value):
    observations = _observations()
    observations["dry_air_subcolumns"][0, 0] = invalid_value

    with pytest.raises(
        ValueError, match="dry-air subcolumn must be finite and strictly positive"
    ):
        structured_superobservations_to_dataset(
            observations, "CH4", "input.nc"
        )


def test_canonical_validator_accepts_dataset_without_optional_provenance():
    dataset = structured_superobservations_to_dataset(
        _observations(), "CH4", "input.nc"
    )
    optional_fields = OPTIONAL_SUPEROBSERVATION_VARIABLES.intersection(
        dataset.variables
    )

    validate_superobservation_dataset(dataset.drop_vars(optional_fields))


def test_converter_accepts_observations_without_optional_coordinates_or_indices():
    observations = _observations()
    retained_fields = [
        name for name in observations.dtype.names
        if name not in {"lat_sat", "lon_sat", "iGC", "jGC"}
    ]

    dataset = structured_superobservations_to_dataset(
        observations[retained_fields], "CH4", "input.nc"
    )

    validate_superobservation_dataset(dataset)
    assert "satellite_latitude" not in dataset
    assert "model_i" not in dataset


def test_coordinate_edges_support_nonuniform_centers():
    np.testing.assert_allclose(
        coordinate_edges([0.0, 1.0, 3.0, 6.0]),
        [-0.5, 0.5, 2.0, 4.5, 7.5],
    )


def test_goopy_parser_reads_canonical_dataset(tmp_path):
    path = tmp_path / "superobs.nc"
    structured_superobservations_to_dataset(
        _observations(), "CH4", "input.nc"
    ).to_netcdf(path)

    parser_path = Path(__file__).parents[3] / "GOOPy" / "parsers.py"
    sys.path.insert(0, str(parser_path.parent))
    spec = importlib.util.spec_from_file_location("goopy_parsers", parser_path)
    parsers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parsers)
    parsed = parsers.read_IMI_superobservations(path, {})

    assert parsed.sizes == {"N_OBS": 1, "N_CENTERS": 2, "N_EDGES": 3}
    np.testing.assert_allclose(parsed["PRESSURE_EDGES"], [[1000, 500, 0.1]])
    np.testing.assert_allclose(parsed["SATELLITE_COLUMN"], [1.9e-6])


def test_methanesat_averager_chunks_and_uses_joint_mask(tmp_path):
    state_vector_path = tmp_path / "StateVector.nc"
    xr.Dataset(
        {"StateVector": (("lat", "lon"), np.ones((2, 2)))},
        coords={"lat": [0.5, 1.5], "lon": [0.5, 1.5]},
    ).to_netcdf(state_vector_path)

    msat_path = tmp_path / "msat.nc"
    base_time = 1_730_000_000.0
    with nc.Dataset(msat_path, "w") as dataset:
        dataset.createDimension("lat", 4)
        dataset.createDimension("lon", 4)
        dataset.createVariable("lat", "f8", ("lat",))[:] = [.25, .75, 1.25, 1.75]
        dataset.createVariable("lon", "f8", ("lon",))[:] = [.25, .75, 1.25, 1.75]
        xch4 = dataset.createVariable("xch4", "f4", ("lat", "lon"), fill_value=1e36)
        time = dataset.createVariable("time", "f8", ("lat", "lon"), fill_value=0.0)
        time.units = "seconds since 1970-1-1 0:0:0, in UTC"
        xch4[:] = np.arange(16, dtype=np.float32).reshape(4, 4) + 1800
        time[:] = base_time
        apriori = dataset.createGroup("apriori_data")
        pressure = apriori.createVariable(
            "surface_pressure", "f4", ("lat", "lon"), fill_value=1e36
        )
        pressure[:] = 1000.0

        # Three different invalid fields leave only pixel [1, 1] in cell [0, 0].
        xch4[0, 0] = 1e36
        time[0, 1] = 0.0
        pressure[1, 0] = 1e36

    observations = average_methanesat_observations(
        msat_path,
        state_vector_path,
        gc_startdate="2024-10-26",
        gc_enddate="2024-10-27",
        row_chunk_size=1,
    )

    assert len(observations) == 4
    first = observations[(observations["jGC"] == 0) & (observations["iGC"] == 0)][0]
    assert first["observation_count"] == 1
    np.testing.assert_allclose(first["CH4"], 1805.0)
    np.testing.assert_allclose(first["surface_pressure"], 1000.0)
    assert np.all(np.diff(first["p_sat"]) < 0)

    target_observations = average_methanesat_observations(
        msat_path,
        state_vector_path,
        gc_startdate="2024-10-26",
        gc_enddate="2024-10-27",
        target_lats=np.array([0.5, 1.5]),
        target_lons=np.array([0.5, 1.0]),
    )
    assert target_observations["iGC"].max() < 2
    assert set(target_observations["lon"]) <= {0.5, 1.0}

    parsed = read_goopy_msat(
        msat_path,
        {
            "N_OBS": "none", "N_EDGES": "none", "N_CENTERS": "none",
            "PRESSURE_EDGES": "none", "PRESSURE_WEIGHT": "none",
            "LATITUDE": "lat", "LONGITUDE": "lon", "TIME": "time",
            "AVERAGING_KERNEL": "none", "PRIOR_PROFILE": "none",
            "SATELLITE_COLUMN": "xch4", "QUALITY_FLAG": "none",
        },
    )
    assert np.issubdtype(parsed["TIME"].dtype, np.datetime64)


def test_rectilinear_grid_index_validation_reports_mismatch():
    with np.testing.assert_raises_regex(
        ValueError,
        r"iGC=\[6, 26\].*GEOS-Chem shape=\(31, 26\).*invalid observations=1",
    ):
        _ravel_rectilinear_grid_indices(
            np.array([2, 15]),
            np.array([6, 26]),
            (31, 26),
        )
