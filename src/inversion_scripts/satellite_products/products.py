"""Built-in IMI satellite-product strategies."""

from __future__ import annotations

import datetime
import os
import re
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pygeohash as pgh

from GOOPy.parsers import read_MSAT
from src.inversion_scripts.operators.msat_funcs import average_methanesat_observations
from src.inversion_scripts.operators.operator_utilities import get_gc_lat_lon
from src.inversion_scripts.satellite_products.base import (
    ObservationRequest,
    SatelliteProduct,
    SuperobservationResult,
)
from src.inversion_scripts.utils import (
    extract_observation_date,
    filter_MSAT,
    filter_blended,
    filter_tropomi,
    read_blended,
    read_tropomi,
)


def _preview_from_raw(
    satellite: dict, indices: np.ndarray, species: str
) -> dict:
    """Convert filtered raw observations to the common preview dictionary."""
    preview = {
        "lat": [],
        "lon": [],
        species: [],
        "swir_albedo": [],
        "time": [],
        "observation_count": [],
    }
    for lat_index, lon_index in zip(*indices):
        preview["lat"].append(satellite["latitude"][lat_index, lon_index])
        preview["lon"].append(satellite["longitude"][lat_index, lon_index])
        preview[species].append(satellite[species][lat_index, lon_index])
        preview["swir_albedo"].append(
            satellite["swir_albedo"][lat_index, lon_index]
        )
        preview["time"].append(satellite["time"][lat_index, lon_index])
        preview["observation_count"].append(1.0)
    return preview


class TropomiFamilyProduct(SatelliteProduct):
    """Shared implementation for operational and blended TROPOMI products."""

    reader = staticmethod(read_tropomi)
    observation_filter = staticmethod(filter_tropomi)

    def validate_directory(self, directory: Path) -> None:
        files = sorted(directory.glob("*.nc"))
        orbit_numbers = []
        for path in files:
            matches = re.findall(r"_(\d{5})_", path.name)
            if len(matches) != 1:
                raise ValueError(f"Please check the TROPOMI filename format: {path}")
            orbit_numbers.append(matches[0])
        if len(set(orbit_numbers)) != len(files):
            raise ValueError(
                "Duplicate TROPOMI orbit numbers found; retain only one "
                "processor version per orbit"
            )

    def read_and_filter(
        self, request: ObservationRequest
    ) -> tuple[dict, np.ndarray] | None:
        satellite = self.reader(request.filename)
        if satellite is None:
            print(f"Skipping {request.filename} due to file processing issue.")
            return None
        indices = self.observation_filter(
            satellite,
            request.xlim,
            request.ylim,
            request.start_date,
            request.end_date,
            request.use_water_observations,
        )
        return satellite, indices

    def create_superobservations(
        self, request: ObservationRequest
    ) -> SuperobservationResult | None:
        # Local imports avoid coupling low-level product definitions to the
        # high-level operator during module initialization.
        from src.inversion_scripts.operators.satellite_operator import (
            average_satellite_observations,
            average_satellite_observations_to_CSgrid,
        )

        result = self.read_and_filter(request)
        if result is None:
            return None
        satellite, indices = result
        observation_count = len(indices[0])
        if observation_count == 0:
            print(f"No satellite observations found in {request.filename}. Skipping.")
            return None
        print("Found", observation_count, "satellite observations.")

        if request.config is None or request.gc_cache is None or request.time_threshold is None:
            raise ValueError("Superobservation generation requires config, gc_cache, and time_threshold")

        config = request.config
        if config["UseGCHP"]:
            if config["STRETCH_GRID"]:
                stretch = f"{config['STRETCH_FACTOR']:.2f}".replace(".", "d")
                target = pgh.encode(config["TARGET_LAT"], config["TARGET_LON"])
                gridspec_path = (
                    f"c{config['CS_RES']}_s{stretch}_t{target}_gridspec.nc"
                )
            else:
                gridspec_path = f"c{config['CS_RES']}_gridspec.nc"
            grid_shape = (6, config["CS_RES"], config["CS_RES"])
            cs_grid_dir = os.path.join(
                os.path.expandvars(config["OutputPath"]),
                config["RunName"],
                "CS_grids",
            )
            observations = average_satellite_observations_to_CSgrid(
                satellite,
                request.species,
                str(request.filename),
                indices,
                request.time_threshold,
                cs_grid_dir,
                gridspec_path,
                grid_shape,
            )
        else:
            gc_lat_lon = get_gc_lat_lon(request.gc_cache, request.start_date)
            observations = average_satellite_observations(
                satellite,
                request.species,
                gc_lat_lon,
                indices,
                request.time_threshold,
            )
            grid_shape = (len(gc_lat_lon["lat"]), len(gc_lat_lon["lon"]))
        return SuperobservationResult(observations, grid_shape)

    def preview(self, request: ObservationRequest) -> dict | None:
        result = self.read_and_filter(request)
        if result is None:
            return None
        return _preview_from_raw(*result, request.species)


class TropomiProduct(TropomiFamilyProduct):
    name = "TROPOMI"
    goopy_raw_product_name = "TROPOMI"

    def download(self, start_date: str, end_date: str, destination: Path) -> None:
        from src.utilities.download_TROPOMI import download_operational_TROPOMI

        start = datetime.datetime.strptime(start_date, "%Y%m%d")
        end = datetime.datetime.strptime(end_date, "%Y%m%d")
        download_operational_TROPOMI(start, end, str(destination))


class BlendedTropomiProduct(TropomiFamilyProduct):
    name = "BlendedTROPOMI"
    goopy_raw_product_name = "TROPOMI_blended"
    reader = staticmethod(read_blended)
    observation_filter = staticmethod(filter_blended)

    def download(self, start_date: str, end_date: str, destination: Path) -> None:
        from src.utilities.download_blended_TROPOMI import download_blended

        start = datetime.datetime.strptime(start_date, "%Y%m%d")
        end = datetime.datetime.strptime(end_date, "%Y%m%d")
        download_blended(start, end, str(destination))


class MethaneSatProduct(SatelliteProduct):
    name = "MSAT"
    goopy_raw_product_name = "MSAT"
    visualization_source = "superobservation"

    def validate_directory(self, directory: Path) -> None:
        for path in sorted(directory.glob("*.nc")):
            extract_observation_date(path)
            with nc.Dataset(path) as dataset:
                missing = {"lat", "lon", "time", "xch4"}.difference(
                    dataset.variables
                )
                if missing:
                    raise ValueError(
                        f"{path} is missing MSAT variables: {sorted(missing)}"
                    )
                apriori = dataset.groups.get("apriori_data")
                if apriori is None or "surface_pressure" not in apriori.variables:
                    raise ValueError(
                        f"{path} is missing apriori_data/surface_pressure"
                    )

    def read_and_filter(
        self, request: ObservationRequest
    ) -> tuple[dict, np.ndarray] | None:
        satellite = read_MSAT(request.filename)
        if satellite is None:
            print(f"Skipping {request.filename} due to file processing issue.")
            return None
        indices = filter_MSAT(
            satellite,
            request.xlim,
            request.ylim,
            request.start_date,
            request.end_date,
            request.use_water_observations,
        )
        return satellite, indices

    def create_superobservations(
        self, request: ObservationRequest
    ) -> SuperobservationResult:
        if request.config is not None and request.config["UseGCHP"]:
            raise NotImplementedError(
                "MethaneSAT superobservations currently require a rectilinear grid"
            )
        if request.state_vector_path is None:
            raise ValueError("MSAT superobservations require the state-vector path")
        if request.gc_cache is None:
            raise ValueError("MSAT superobservations require gc_cache")
        observations = average_methanesat_observations(
            request.filename,
            request.state_vector_path,
            species=request.species,
            time_threshold=request.time_threshold,
            gc_startdate=request.start_date,
            gc_enddate=request.end_date,
        )
        gc_lat_lon = get_gc_lat_lon(request.gc_cache, request.start_date)
        grid_shape = (len(gc_lat_lon["lat"]), len(gc_lat_lon["lon"]))
        return SuperobservationResult(observations, grid_shape)

    def preview(self, request: ObservationRequest) -> dict:
        if request.state_vector_path is None:
            raise ValueError("MSAT preview requires the state-vector path")
        observations = average_methanesat_observations(
            request.filename,
            request.state_vector_path,
            species=request.species,
            gc_startdate=request.start_date,
            gc_enddate=request.end_date,
        )
        count = len(observations)
        times = [
            datetime.datetime.strptime(value, "%Y%m%d_%H")
            for value in observations["time"]
        ]
        return {
            "lat": observations["lat_sat"].tolist(),
            "lon": observations["lon_sat"].tolist(),
            request.species: observations[request.species].tolist(),
            "swir_albedo": [np.nan] * count,
            "time": times,
            "observation_count": observations["observation_count"].tolist(),
        }
