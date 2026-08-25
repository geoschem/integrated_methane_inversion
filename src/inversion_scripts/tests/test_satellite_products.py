"""Tests for centralized satellite-product registration."""

import datetime

import numpy as np

import pytest

from src.inversion_scripts.satellite_products import (
    ObservationRequest,
    get_satellite_product,
    supported_satellite_products,
)


def test_builtin_satellite_products_are_registered():
    assert supported_satellite_products() == (
        "BlendedTROPOMI",
        "MSAT",
        "TROPOMI",
    )


@pytest.mark.parametrize(
    "name,goopy_name,visualization_source",
    [
        ("TROPOMI", "TROPOMI", "raw"),
        ("BlendedTROPOMI", "TROPOMI_blended", "raw"),
        ("MSAT", "MSAT", "superobservation"),
    ],
)
def test_product_metadata(name, goopy_name, visualization_source):
    product = get_satellite_product(name)
    assert product.goopy_raw_product_name == goopy_name
    assert product.visualization_source == visualization_source


def test_unknown_satellite_product_lists_supported_names():
    with pytest.raises(ValueError, match="Unsupported satellite product 'UNKNOWN'"):
        get_satellite_product("UNKNOWN")


def test_msat_preview_returns_datetime_values(monkeypatch):
    observations = np.zeros(
        1,
        dtype=[
            ("lat_sat", "f4"),
            ("lon_sat", "f4"),
            ("CH4", "f4"),
            ("time", "U13"),
            ("observation_count", "f4"),
        ],
    )
    observations["time"] = "20241026_21"
    monkeypatch.setattr(
        "src.inversion_scripts.satellite_products.products."
        "average_methanesat_observations",
        lambda *args, **kwargs: observations,
    )
    request = ObservationRequest(
        filename="msat.nc",
        species="CH4",
        start_date=np.datetime64("2024-10-26"),
        end_date=np.datetime64("2024-10-27"),
        xlim=(-1.0, 1.0),
        ylim=(-1.0, 1.0),
        state_vector_path="StateVector.nc",
    )

    preview = get_satellite_product("MSAT").preview(request)

    assert preview["time"] == [datetime.datetime(2024, 10, 26, 21)]
