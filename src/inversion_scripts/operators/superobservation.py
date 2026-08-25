"""Canonical, product-neutral IMI superobservation files."""

import numpy as np
import pandas as pd
import xarray as xr


FORMAT_NAME = "IMI_superobservation"
FORMAT_VERSION = "1.0"


def imi_superobservation_dtype(
    species: str,
    n_pressure_edges: int,
    n_layers: int,
) -> np.dtype:
    """Return the canonical IMI rectilinear superobservation NumPy dtype.

    Args:
        species: Name used for the retrieved column field, such as ``CH4``.
        n_pressure_edges: Number of vertical pressure edges per observation.
        n_layers: Number of layer-resolved values per observation.

    Returns:
        Structured NumPy dtype shared by satellite-specific rectilinear
        averagers and the canonical NetCDF writer.
    """
    if n_pressure_edges < 0 or n_layers < 0:
        raise ValueError("Superobservation vertical dimensions cannot be negative")
    return np.dtype([
        ("iGC", "i4"), ("jGC", "i4"),
        ("lat_sat", "f4"), ("lon_sat", "f4"),
        (species, "f4"), ("time", "U13"),
        ("p_sat", "f4", (n_pressure_edges,)),
        ("surface_pressure", "f4"),
        ("nir_albedo", "f4"), ("swir_albedo", "f4"),
        ("dry_air_subcolumns", "f4", (n_layers,)),
        ("apriori", "f4", (n_layers,)),
        ("avkern", "f4", (n_layers,)),
        ("layer", "f4", (n_layers,)),
        ("observation_count", "f4"),
        ("lat", "f4"), ("lon", "f4"),
    ])


def structured_superobservations_to_dataset(
    observations, species, source_file, source_product="unknown"
) -> xr.Dataset:
    """Convert the IMI in-memory representation to the canonical NetCDF schema."""
    if len(observations) == 0:
        raise ValueError("Cannot create a superobservation file with no observations")

    n_obs = len(observations)
    n_layers = observations["avkern"].shape[1]
    n_edges = observations["p_sat"].shape[1]
    if n_edges != n_layers + 1:
        raise ValueError(
            "Superobservation pressure_edges must have one more element than "
            f"the layer fields; got {n_edges} edges and {n_layers} layers"
        )

    dry_air = np.asarray(observations["dry_air_subcolumns"], dtype=np.float64)
    dry_air_total = dry_air.sum(axis=1, keepdims=True)
    if np.any(~np.isfinite(dry_air_total)) or np.any(dry_air_total <= 0):
        raise ValueError("Every superobservation must have a positive dry-air column")

    # The legacy IMI structure stores the prior as a partial column (mol m-2).
    prior_profile = np.divide(
        np.asarray(observations["apriori"], dtype=np.float64),
        dry_air,
        out=np.full_like(dry_air, np.nan),
        where=dry_air > 0,
    )
    pressure_weight = dry_air / dry_air_total

    times = pd.to_datetime(
        observations["time"], format="%Y%m%d_%H"
    ).to_numpy(dtype="datetime64[ns]")

    field_names = set(observations.dtype.names)
    latitude = observations["lat"] if "lat" in field_names else observations["lat_sat"]
    longitude = observations["lon"] if "lon" in field_names else observations["lon_sat"]
    data_vars = {
            "latitude": ("observation", latitude.astype(np.float64)),
            "longitude": ("observation", longitude.astype(np.float64)),
            "satellite_latitude": (
                "observation", observations["lat_sat"].astype(np.float64)
            ),
            "satellite_longitude": (
                "observation", observations["lon_sat"].astype(np.float64)
            ),
            "time": ("observation", times),
            "column": (
                "observation", observations[species].astype(np.float64) * 1e-9
            ),
            # Canonical files use Pa. The legacy IMI pressure arrays use hPa.
            "pressure_edges": (
                ("observation", "edge"),
                observations["p_sat"].astype(np.float64) * 100.0,
            ),
            "pressure_weight": (
                ("observation", "layer"), pressure_weight
            ),
            "averaging_kernel": (
                ("observation", "layer"),
                observations["avkern"].astype(np.float64),
            ),
            "prior_profile": (
                ("observation", "layer"), prior_profile
            ),
            "dry_air_subcolumn": (
                ("observation", "layer"), dry_air
            ),
            "observation_count": (
                "observation", observations["observation_count"].astype(np.float64)
            ),
        }
    if {"iGC", "jGC"}.issubset(field_names):
        data_vars.update({
            "model_i": ("observation", observations["iGC"].astype(np.int32)),
            "model_j": ("observation", observations["jGC"].astype(np.int32)),
        })
    elif {"nfi", "Ydimi", "Xdimi"}.issubset(field_names):
        data_vars.update({
            "model_face": ("observation", observations["nfi"].astype(np.int32)),
            "model_y": ("observation", observations["Ydimi"].astype(np.int32)),
            "model_x": ("observation", observations["Xdimi"].astype(np.int32)),
        })
    else:
        raise ValueError("Superobservations lack model grid indices")

    dataset = xr.Dataset(
        data_vars=data_vars,
        coords={
            "observation": np.arange(n_obs, dtype=np.int64),
            "layer": np.arange(n_layers, dtype=np.int32),
            "edge": np.arange(n_edges, dtype=np.int32),
        },
        attrs={
            "format_name": FORMAT_NAME,
            "format_version": FORMAT_VERSION,
            "source_file": str(source_file),
            "source_product": str(source_product),
            "species": species,
            "title": "Grid-cell-averaged satellite superobservations",
        },
    )

    units = {
        "latitude": "degrees_north",
        "longitude": "degrees_east",
        "satellite_latitude": "degrees_north",
        "satellite_longitude": "degrees_east",
        "column": "mol mol-1",
        "pressure_edges": "Pa",
        "pressure_weight": "1",
        "averaging_kernel": "1",
        "prior_profile": "mol mol-1",
        "dry_air_subcolumn": "mol m-2",
        "observation_count": "1",
    }
    for variable, unit in units.items():
        dataset[variable].attrs["units"] = unit

    validate_superobservation_dataset(dataset)
    return dataset


def validate_superobservation_dataset(dataset):
    """Validate the canonical fields used by the GOOPy observation operator."""
    required = {
        "latitude", "longitude", "time", "column", "pressure_edges",
        "pressure_weight", "averaging_kernel", "prior_profile",
        "observation_count",
    }
    missing = required.difference(dataset.variables)
    if missing:
        raise ValueError(f"Missing canonical superobservation fields: {sorted(missing)}")
    if dataset.attrs.get("format_name") != FORMAT_NAME:
        raise ValueError(f"format_name must be {FORMAT_NAME!r}")
    if dataset.sizes["edge"] != dataset.sizes["layer"] + 1:
        raise ValueError("The edge dimension must be one longer than layer")
    if not np.allclose(dataset["pressure_weight"].sum("layer"), 1.0):
        raise ValueError("pressure_weight must sum to one for each observation")
    if np.any(dataset["observation_count"].values <= 0):
        raise ValueError("observation_count must be positive")
