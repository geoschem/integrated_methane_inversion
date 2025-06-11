#!/usr/bin/env python3
"""
Classify observations into a cubed-sphere grid based on geographic location.

"""

import numpy as np
import xarray as xr
from shapely.geometry import Point, Polygon
from scipy.spatial import cKDTree


def latlon_to_cartesian(lat, lon):
    """Convert latitude and longitude to Cartesian coordinates on a unit sphere."""
    lat_rad = np.radians(lat)
    lon_rad = np.radians(lon)
    x = np.cos(lat_rad) * np.cos(lon_rad)
    y = np.cos(lat_rad) * np.sin(lon_rad)
    z = np.sin(lat_rad)
    return np.vstack((x, y, z)).T


def build_kdtree(lats, lons):
    """Build KDTree from flattened simulation grid."""
    lat_flat = lats.reshape(-1)
    lon_flat = lons.reshape(-1)
    cart_coords = latlon_to_cartesian(lat_flat, lon_flat)
    return cKDTree(cart_coords), lats.shape


def precompute_polygons(corner_lats, corner_lons):
    """Precompute simulation grid cell polygons from corner coordinates."""
    corner_lons[corner_lons>180] -= 360
    nf, YCdim, XCdim = corner_lats.shape
    Ydim = YCdim - 1
    Xdim = XCdim - 1
    polygons = np.empty((nf, Ydim, Xdim), dtype=object)

    for f in range(nf):
        for j in range(Ydim):
            for i in range(Xdim):
                try:
                    polygons[f, j, i] = Polygon([
                        (corner_lons[f, j, i],     corner_lats[f, j, i]),
                        (corner_lons[f, j, i + 1], corner_lats[f, j, i + 1]),
                        (corner_lons[f, j + 1, i + 1], corner_lats[f, j + 1, i + 1]),
                        (corner_lons[f, j + 1, i], corner_lats[f, j + 1, i])
                    ])
                except Exception:
                    polygons[f, j, i] = None
    return polygons


def classify_obs_to_cs_grid(
    obs_dict: dict,
    grid_path: str,
    k: int = 1
):
    """
    Classify observation points into a cubed-sphere grid.

    Parameters:
    - obs_dict: filtered TROPOMI observations in dict format
    - grid_path: Path to cubed-sphere grid NetCDF file
    - k: Number of KDTree neighbors to consider for polygon check

    Returns:
    - obs_super: classified observation in CS grids in dict format
    """
    # Load datasets
    grid_ds = xr.open_dataset(grid_path)

    # Read in values
    obs_lat = np.array(obs_dict['lat'])
    obs_lon = np.array(obs_dict['lon'])
    obs_lon[obs_lon>180] -= 360

    # Read simulation grid
    lats = np.array(grid_ds["lats"])
    lons = np.array(grid_ds["lons"])
    corner_lats = np.array(grid_ds["corner_lats"])
    corner_lons = np.array(grid_ds["corner_lons"])

    # Build KDTree and polygons
    kdtree, shape = build_kdtree(lats, lons)
    polygons = precompute_polygons(corner_lats, corner_lons)

    # Cartesian coords of obs
    obs_cart = latlon_to_cartesian(obs_lat, obs_lon)
    _, neighbor_idxs = kdtree.query(obs_cart, k=k)
    if k == 1:
        neighbor_idxs = np.expand_dims(neighbor_idxs, axis=1)

    # Classify
    obs_super = obs_dict[["lat", "lon", "time", "obs_count"]].copy()
    obs_super["date"] = obs_super["time"].dt.floor("D")

    for i in range(len(obs_lat)):
        pt = Point(obs_lon[i], obs_lat[i])
        for idx in neighbor_idxs[i]:
            f, j, x = np.unravel_index(idx, shape)
            poly = polygons[f, j, x]
            if poly is not None and poly.contains(pt):
                obs_super["lat"][i] = lats[f,j,x]
                obs_super["lon"][i] = lons[f,j,x]

    return obs_super