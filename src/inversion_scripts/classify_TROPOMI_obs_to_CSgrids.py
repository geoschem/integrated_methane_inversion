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
    lons[lons>180] -= 360
    corner_lats = grid_ds["corner_lats"].values
    corner_lons = grid_ds["corner_lons"].values

    # Build KDTree and polygons
    kdtree, shape = build_kdtree(lats, lons)

    # Cartesian coords of obs
    obs_cart = latlon_to_cartesian(obs_lat, obs_lon)
    _, neighbor_idxs = kdtree.query(obs_cart, k=k)

    # Create new DataFrame with modified lat/lon and additional 'date'
    obs_super = obs_df.copy()
    obs_super['sim_index'] = neighbor_idxs
    obs_super['date'] = pd.to_datetime(obs_df['time']).dt.floor("D")
    obs_super.drop(columns='time', inplace=True)

    return obs_super

def map_obs_to_CSgrid(obs_super, gridpath):
    """
    Map observation counts from obs_super DataFrame onto
    an xarray Dataset grid based on (nf, Ydim, Xdim) indices and date.

    Parameters:
    - obs_super: pd.DataFrame with columns 'date', 'sim_index', and 'obs_count'
    - gridpath: path to NetCDF grid file
    
    Returns:
    - obs_dataset: xarray.Dataset with 'obs_count' on ('date', 'nf', 'Ydim', 'Xdim'),
                   and associated coordinates: lats, lons, corner_lats, corner_lons.
    """

    grouped = obs_super.groupby(['date', 'sim_index'])['obs_count'].sum().reset_index()

    gridds = xr.open_dataset(gridpath)
    nf = gridds['nf'].values
    Ydim = gridds['Ydim'].values
    Xdim = gridds['Xdim'].values
    lats = gridds['lats']
    lons = gridds['lons']
    corner_lons = gridds['corner_lons']
    corner_lats = gridds['corner_lats']
    grid_shape = lats.shape
    dates = sorted(grouped['date'].unique())

    date_to_idx = {date: i for i, date in enumerate(dates)}
    obs_count_grid = np.full((len(dates), *grid_shape), np.nan, dtype=float)

    for _, row in grouped.iterrows():
        date_idx = date_to_idx[row['date']]
        f_idx, j_idx, x_idx = np.unravel_index(row['sim_index'], grid_shape)
        obs_count_grid[date_idx, f_idx, j_idx, x_idx] = row['obs_count']

    obs_count_da = xr.DataArray(
        obs_count_grid,
        dims=['date', 'nf', 'Ydim', 'Xdim'],
        coords=dict(
            date=(['date'], dates),
            nf=(['nf'], nf),
            Ydim=(['Ydim'], Ydim),
            Xdim=(['Xdim'], Xdim),
            lats=(['nf', 'Ydim', 'Xdim'], lats.values),
            lons=(['nf', 'Ydim', 'Xdim'], lons.values)
        ),
        attrs=dict(long_name='observation count', units='1')
    )

    # Set coordinate metadata
    obs_count_da.coords['lats'].attrs.update({"units": "degrees_north", "long_name": "Latitude"})
    obs_count_da.coords['lons'].attrs.update({"units": "degrees_east", "long_name": "Longitude"})

    obs_dataset = xr.Dataset({
        'obs_count': obs_count_da,
        'corner_lats': corner_lats,
        'corner_lons': corner_lons
    })

    return obs_dataset