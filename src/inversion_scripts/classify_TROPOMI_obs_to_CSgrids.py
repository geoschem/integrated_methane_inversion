#!/usr/bin/env python3
"""
Classify observations into a cubed-sphere grid based on geographic location.

"""

import numpy as np
import pandas as pd
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
    """
    Precompute simulation grid cell polygons from corner coordinates.

    Parameters:
    - corner_lats, corner_lons: numpy arrays of shape (nf, YC+1, XC+1)
    - project: optional projection function (lon, lat) -> (x, y)

    Returns:
    - polygons: numpy array of shape (nf, YC, XC) with shapely.Polygon objects
    """
    corner_lons = np.where(corner_lons > 180, corner_lons - 360, corner_lons)
    nf, YCdim, XCdim = corner_lats.shape
    Ydim, Xdim = YCdim - 1, XCdim - 1
    polygons = np.empty((nf, Ydim, Xdim), dtype=object)

    for f in range(nf):
        for j in range(Ydim):
            for i in range(Xdim):
                try:
                    lons = [
                        corner_lons[f, j, i],
                        corner_lons[f, j, i + 1],
                        corner_lons[f, j + 1, i + 1],
                        corner_lons[f, j + 1, i],
                    ]
                    lats = [
                        corner_lats[f, j, i],
                        corner_lats[f, j, i + 1],
                        corner_lats[f, j + 1, i + 1],
                        corner_lats[f, j + 1, i],
                    ]

                    # Handle antimeridian: wrap lons to avoid crossing discontinuity
                    lons = np.array(lons)
                    if np.ptp(lons) > 180:
                        lons = np.where(lons > 0, lons - 360, lons)

                    lonlat = list(zip(lons, lats))
                    poly = Polygon(lonlat)

                    # Optional: fix invalid polygons
                    if not poly.is_valid:
                        poly = poly.buffer(0)

                    # Sanity check: make sure we really have a shapely Polygon
                    assert isinstance(poly, Polygon), f"Expected Polygon but got {type(poly)} at f={f}, j={j}, i={i}"

                    polygons[f, j, i] = poly if poly.is_valid else None

                except Exception as e:
                    polygons[f, j, i] = None

    return polygons

def classify_obs_to_cs_grid(
    obs_df: pd.DataFrame,
    grid_path: str,
    k: int = 1
) -> pd.DataFrame:
    """
    Classify observation points into a cubed-sphere grid.

    Parameters:
    - obs_df: filtered TROPOMI observations as a pandas DataFrame
              with at least 'lat', 'lon', 'time', 'obs_count' columns
    - grid_path: Path to cubed-sphere grid NetCDF file
    - k: Number of KDTree neighbors to consider for polygon check

    Returns:
    - obs_super: DataFrame with observations mapped to CS grids
    """
    # Load datasets
    grid_ds = xr.open_dataset(grid_path)

    # Read in observation values
    obs_lat = obs_df['lat'].values
    obs_lon = obs_df['lon'].values
    obs_lon = np.where(obs_lon > 180, obs_lon - 360, obs_lon)

    # Read cubed-sphere grid
    nf = grid_ds["nf"].values
    Ydim = grid_ds["Ydim"].values
    Xdim = grid_ds["Xdim"].values
    lats = grid_ds["lats"].values
    lons = grid_ds["lons"].values
    corner_lats = grid_ds["corner_lats"].values
    corner_lons = grid_ds["corner_lons"].values

    # Build KDTree and polygons
    kdtree, shape = build_kdtree(lats, lons)
    polygons = precompute_polygons(corner_lats, corner_lons)

    # Cartesian coords of obs
    obs_cart = latlon_to_cartesian(obs_lat, obs_lon)
    _, neighbor_idxs = kdtree.query(obs_cart, k=k)
    if k == 1:
        neighbor_idxs = np.expand_dims(neighbor_idxs, axis=1)

    # Prepare output lat/lon arrays
    cs_lat = np.full_like(obs_lat, np.nan, dtype=float)
    cs_lon = np.full_like(obs_lon, np.nan, dtype=float)
    nfi = np.full_like(obs_lat, np.nan, dtype=float)
    yi = np.full_like(obs_lat, np.nan, dtype=float)
    xi = np.full_like(obs_lat, np.nan, dtype=float)

    # Classify each observation point
    for i in range(len(obs_lat)):
        pt = Point(obs_lon[i], obs_lat[i])
        for idx in neighbor_idxs[i]:
            f, j, x = np.unravel_index(idx, shape)
            poly = polygons[f, j, x]
            if poly is not None and poly.contains(pt):
                cs_lat[i] = lats[f, j, x]
                cs_lon[i] = lons[f, j, x]
                nfi[i] = nf[f]
                yi[i] = Ydim[j]
                xi[i] = Xdim[x]
                break

    # Create new DataFrame with modified lat/lon and additional 'date'
    obs_super = obs_df.copy()
    obs_super['lat'] = cs_lat
    obs_super['lon'] = cs_lon
    obs_super['nf'] = nfi
    obs_super['Ydim'] = yi
    obs_super['Xdim'] = xi
    obs_super['date'] = pd.to_datetime(obs_df['time']).dt.floor("D")

    return obs_super

def map_obs_to_CSgrid(obs_super, state_vector_labels, value_cols):
    """
    Map multiple observation columns from obs_super DataFrame onto
    an xarray Dataset grid based on matching (nf, Ydim, Xdim) indices and date.

    Parameters:
    - obs_super: pd.DataFrame with columns 'nf', 'Ydim', 'Xdim', 'date', plus value_cols
    - state_vector_labels: xarray.DataArray with dims ('nf', 'Ydim', 'Xdim')
    - value_cols: list of column names in obs_super to map, e.g. ['lat', 'lon', 'xch4']

    Returns:
    - obs_dataset: xarray.Dataset with DataArrays for each value_col + 'obs_count',
                   aligned on ('nf', 'Ydim', 'Xdim', 'date'), NaN elsewhere.
    """

    # First: aggregate obs_count + values
    # Exclude 'obs_count' if accidentally included in value_cols
    value_cols = [col for col in value_cols if col != 'obs_count']

    # Aggregate value columns
    agg_dict = {col: 'mean' for col in value_cols}
    grouped = obs_super.groupby(['nf', 'Ydim', 'Xdim', 'date']).agg(agg_dict).reset_index()

    # Add observation count explicitly
    grouped['obs_count'] = (
        obs_super.groupby(['nf', 'Ydim', 'Xdim', 'date']).size().values
    )

    # Build the full grid
    nf_vals = state_vector_labels['nf'].values
    Ydim_vals = state_vector_labels['Ydim'].values
    Xdim_vals = state_vector_labels['Xdim'].values
    date_vals = np.sort(obs_super['date'].unique())

    all_index = pd.MultiIndex.from_product(
        [nf_vals, Ydim_vals, Xdim_vals, date_vals],
        names=["nf", "Ydim", "Xdim", "date"]
    )
    full_df = pd.DataFrame(index=all_index).reset_index()

    # Merge full grid with grouped data
    merged = pd.merge(
        full_df,
        grouped,  # already aggregated
        on=["nf", "Ydim", "Xdim", "date"],
        how="left"
    )

    # Convert to xarray
    merged = merged.set_index(["nf", "Ydim", "Xdim", "date"])
    obs_dataset = merged.to_xarray()

    return obs_dataset