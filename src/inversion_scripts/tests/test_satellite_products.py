"""Tests for centralized satellite-product registration."""

import pytest

from src.inversion_scripts.satellite_products import (
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
