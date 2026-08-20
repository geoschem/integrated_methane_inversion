import importlib.util
from pathlib import Path
import sys

import numpy as np
import netCDF4 as nc
import xarray as xr

from src.inversion_scripts.operators.superobservation import (
    structured_superobservations_to_dataset,
)
from src.inversion_scripts.operators.msat_funcs import (
    average_methanesat_observations,
)
from GOOPy.parsers import read_MSAT as read_goopy_msat


def _observations():
    dtype = [
        ("iGC", "i4"), ("jGC", "i4"),
        ("lat_sat", "f4"), ("lon_sat", "f4"),
        ("CH4", "f4"), ("time", "U13"),
        ("p_sat", "f4", (3,)),
        ("surface_pressure", "f4"),
        ("nir_albedo", "f4"), ("swir_albedo", "f4"),
        ("dry_air_subcolumns", "f4", (2,)),
        ("apriori", "f4", (2,)),
        ("avkern", "f4", (2,)),
        ("layer", "f4", (2,)),
        ("observation_count", "f4"),
        ("lat", "f4"), ("lon", "f4"),
    ]
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
