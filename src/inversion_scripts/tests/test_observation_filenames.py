from datetime import datetime

import pytest

from src.inversion_scripts.satellite_products import get_satellite_product


MSAT_NAME = (
    "regrid_t56_20241026_c03490380_p959_45m_MSAT_L3_45m_c03490380_"
    "p959_v01010000_20241026T210413Z_210445Z.nc"
)
CANONICAL_MSAT_NAME = (
    "MSAT_L3_45m_c04970420_p3385_v01014000_20241124T094110Z_094142Z.nc"
)
TROPOMI_NAME = (
    "S5P_OFFL_L2__CH4____20210704T174744_20210704T192914_19300_"
    "03_020400_20221126T123158.nc"
)
BLENDED_TROPOMI_NAME = (
    "S5P_BLND_L2__CH4____20210704T174744_20210704T192914_19300_"
    "03_020400_20221126T123158.nc"
)


@pytest.mark.parametrize("name", [MSAT_NAME, f"{MSAT_NAME[:-3]}_GCtoSatellite.pkl"])
def test_extract_observation_date_from_msat_names(name):
    assert get_satellite_product("MSAT").observation_date(name) == datetime(2024, 10, 26)


def test_extract_observation_date_from_canonical_msat_name():
    assert get_satellite_product("MSAT").observation_date(
        CANONICAL_MSAT_NAME
    ) == datetime(2024, 11, 24)


def test_extract_observation_date_from_tropomi_name():
    assert get_satellite_product("TROPOMI").observation_date(TROPOMI_NAME) == datetime(2021, 7, 4)


def test_extract_observation_date_from_tropomi_derived_name():
    name = f"{BLENDED_TROPOMI_NAME[:-3]}_GCtoSatellite.pkl"
    assert get_satellite_product("BlendedTROPOMI").observation_date(name) == datetime(
        2021, 7, 4
    )


def test_tropomi_products_reject_each_others_filenames():
    with pytest.raises(ValueError, match="TROPOMI naming convention"):
        get_satellite_product("TROPOMI").observation_date(BLENDED_TROPOMI_NAME)
    with pytest.raises(ValueError, match="BlendedTROPOMI naming convention"):
        get_satellite_product("BlendedTROPOMI").observation_date(TROPOMI_NAME)


def test_extract_observation_date_rejects_unknown_name():
    with pytest.raises(ValueError, match="MSAT naming convention"):
        get_satellite_product("MSAT").observation_date("observation.nc")


def test_tropomi_parser_does_not_accept_an_unstructured_date():
    with pytest.raises(ValueError, match="TROPOMI naming convention"):
        get_satellite_product("TROPOMI").observation_date("observation_20210704.nc")
