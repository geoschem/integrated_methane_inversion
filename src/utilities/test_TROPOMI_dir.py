"""Validate an IMI observation directory for the selected product."""

import re
import sys
from pathlib import Path

import netCDF4 as nc

from src.inversion_scripts.utils import extract_observation_date


def check_for_duplicate_orbit_numbers(observation_dir):
    """Retain the processor-version/orbit check used by TROPOMI products."""
    files = sorted(Path(observation_dir).glob("*.nc"))
    orbit_numbers = []
    for path in files:
        matches = re.findall(r"_(\d{5})_", path.name)
        if len(matches) != 1:
            raise ValueError(
                "Please check the TROPOMI filename format: " f"{path}"
            )
        orbit_numbers.append(matches[0])
    if len(set(orbit_numbers)) != len(files):
        raise ValueError(
            "Duplicate TROPOMI orbit numbers found; retain only one processor "
            "version per orbit"
        )


def check_methanesat_files(observation_dir):
    """Check filenames and the minimal L3 variables used by the MSAT averager."""
    files = sorted(Path(observation_dir).glob("*.nc"))
    filenames = [path.name for path in files]
    if len(set(filenames)) != len(filenames):
        raise ValueError("Duplicate MethaneSAT filenames found")

    for path in files:
        date = extract_observation_date(path).strftime("%Y%m%d")
        with nc.Dataset(path) as dataset:
            missing = {"lat", "lon", "time", "xch4"}.difference(dataset.variables)
            if missing:
                raise ValueError(f"{path} is missing MSAT variables: {sorted(missing)}")
            apriori = dataset.groups.get("apriori_data")
            if apriori is None or "surface_pressure" not in apriori.variables:
                raise ValueError(f"{path} is missing apriori_data/surface_pressure")


def validate_observation_directory(observation_dir, satellite_product):
    files = list(Path(observation_dir).glob("*.nc"))
    if not files:
        raise ValueError(f"No NetCDF observation files found in {observation_dir}")
    if satellite_product in {"TROPOMI", "BlendedTROPOMI"}:
        check_for_duplicate_orbit_numbers(observation_dir)
    elif satellite_product == "MSAT":
        check_methanesat_files(observation_dir)
    else:
        for path in files:
            extract_observation_date(path)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(
            "Usage: test_TROPOMI_dir.py OBSERVATION_DIR SATELLITE_PRODUCT"
        )
    validate_observation_directory(sys.argv[1], sys.argv[2])
