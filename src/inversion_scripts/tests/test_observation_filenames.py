from datetime import datetime

import pytest

from src.inversion_scripts.utils import extract_observation_date


MSAT_NAME = (
    "regrid_t56_20241026_c03490380_p959_45m_MSAT_L3_45m_c03490380_"
    "p959_v01010000_20241026T210413Z_210445Z.nc"
)
TROPOMI_NAME = (
    "S5P_OFFL_L2__CH4____20210704T174744_20210704T192914_19300_"
    "03_020400_20221126T123158.nc"
)


@pytest.mark.parametrize("name", [MSAT_NAME, f"{MSAT_NAME[:-3]}_GCtoSatellite.pkl"])
def test_extract_observation_date_from_msat_names(name):
    assert extract_observation_date(name) == datetime(2024, 10, 26)


def test_extract_observation_date_from_tropomi_name():
    assert extract_observation_date(TROPOMI_NAME) == datetime(2021, 7, 4)


def test_extract_observation_date_rejects_unknown_name():
    with pytest.raises(ValueError, match="acquisition date"):
        extract_observation_date("observation.nc")
