"""Validate an IMI observation directory for the selected product."""

import sys
from pathlib import Path

from src.inversion_scripts.satellite_products import get_satellite_product


def check_for_duplicate_orbit_numbers(observation_dir):
    """Validate TROPOMI filenames and orbit uniqueness."""
    get_satellite_product("TROPOMI").validate_directory(Path(observation_dir))


def check_methanesat_files(observation_dir):
    """Validate the raw files needed by the MSAT averager."""
    get_satellite_product("MSAT").validate_directory(Path(observation_dir))


def validate_observation_directory(observation_dir, satellite_product):
    directory = Path(observation_dir)
    files = list(directory.glob("*.nc"))
    if not files:
        raise ValueError(f"No NetCDF observation files found in {observation_dir}")
    get_satellite_product(satellite_product).validate_directory(directory)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(
            "Usage: validate_observation_directory.py OBSERVATION_DIR SATELLITE_PRODUCT"
        )
    validate_observation_directory(sys.argv[1], sys.argv[2])
