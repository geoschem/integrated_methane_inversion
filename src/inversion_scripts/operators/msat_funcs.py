import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
from src.inversion_scripts.utils import get_strdate


def coordinate_edges(centers: np.ndarray) -> np.ndarray:
    """Infer cell edges from increasing, possibly nonuniform cell centers.

    Interior edges are the midpoint of each adjacent center pair. The two
    exterior edges are extrapolated using only their nearest center spacing,
    so no globally uniform resolution is assumed.

    Example: coordinate_edges([0, 1, 3, 6]) -> [-0.5, 0.5, 2, 4.5, 7.5]
    """
    centers = np.asarray(centers, dtype=np.float64)
    if centers.ndim != 1 or len(centers) < 2:
        raise ValueError("Coordinates must be a 1-D array with at least two cells")
    differences = np.diff(centers)
    if not np.all(differences > 0):
        raise ValueError("Coordinates must be strictly increasing")

    edges = np.empty(len(centers) + 1, dtype=np.float64)
    edges[1:-1] = centers[:-1] + differences / 2
    edges[0] = centers[0] - differences[0] / 2
    edges[-1] = centers[-1] + differences[-1] / 2
    return edges    


def values_with_nan(
    variable: nc.Variable,
    key: tuple[slice, slice],
) -> np.ndarray:
    """Read a 2-D NetCDF slice as float64 with masked values set to NaN.

    Args:
        variable: NetCDF variable containing a two-dimensional MSAT field.
        key: Latitude and longitude slices used to read the current chunk.

    Returns:
        A float64 NumPy array in which NetCDF masked/fill values are NaN.
    """
    values = variable[key]
    if np.ma.isMaskedArray(values):
        values = values.filled(np.nan)
    return np.asarray(values, dtype=np.float64)


def average_methanesat_observations(
    data_path: str,
    state_vector: str,
    species: str = "CH4",
    time_threshold: str = None,
    gc_startdate: str = None,
    gc_enddate: str = None,
    row_chunk_size: int = 256,
):
    """Average a rectilinear MethaneSAT L3 file onto the state-vector grid.

    The input can be much larger than memory. Only coordinate rows intersecting
    the state-vector domain are read, and the two-dimensional retrieval fields
    are accumulated in row chunks. A single joint validity mask is applied to
    XCH4, time, and surface pressure so all averaged fields and counts describe
    the same contributing pixels.
    """

    with xr.open_dataset(state_vector) as state_vector_ds:
        gc_lats = np.asarray(state_vector_ds["lat"].values, dtype=np.float64)
        gc_lons = np.asarray(state_vector_ds["lon"].values, dtype=np.float64)

    if gc_lats.ndim != 1 or gc_lons.ndim != 1:
        raise ValueError("MethaneSAT averaging requires 1-D state-vector lat/lon")
    if len(gc_lats) < 2 or len(gc_lons) < 2:
        raise ValueError("State-vector lat/lon must each contain at least two cells")

    lat_edges = coordinate_edges(gc_lats)
    lon_edges = coordinate_edges(gc_lons)
    n_lat = len(gc_lats)
    n_lon = len(gc_lons)
    n_cells = n_lat * n_lon

    def epoch_seconds(value):
        if value is None:
            return None
        return pd.Timestamp(value).timestamp()

    start_seconds = epoch_seconds(gc_startdate)
    end_timestamp = pd.Timestamp(gc_enddate) if gc_enddate is not None else None
    end_is_date = end_timestamp is not None and end_timestamp == end_timestamp.normalize()
    end_seconds = (
        (end_timestamp + pd.Timedelta(days=1)).timestamp()
        if end_is_date else epoch_seconds(gc_enddate)
    )

    sums_xch4 = np.zeros(n_cells, dtype=np.float64)
    sums_time = np.zeros(n_cells, dtype=np.float64)
    sums_surface_pressure = np.zeros(n_cells, dtype=np.float64)
    sums_latitude = np.zeros(n_cells, dtype=np.float64)
    sums_longitude = np.zeros(n_cells, dtype=np.float64)
    counts = np.zeros(n_cells, dtype=np.int64)



    with nc.Dataset(data_path, "r") as ncfile:
        required_root = {"lat", "lon", "time", "xch4"}
        missing_root = required_root.difference(ncfile.variables)
        if missing_root:
            raise ValueError(f"MSAT file is missing variables: {sorted(missing_root)}")
        if "apriori_data" not in ncfile.groups:
            raise ValueError("MSAT file is missing the apriori_data group")
        if "surface_pressure" not in ncfile.groups["apriori_data"].variables:
            raise ValueError("MSAT file is missing apriori_data/surface_pressure")

        source_lats = np.asarray(ncfile["lat"][:], dtype=np.float64)
        source_lons = np.asarray(ncfile["lon"][:], dtype=np.float64)
        lat_selection = np.flatnonzero(
            (source_lats >= lat_edges[0]) & (source_lats <= lat_edges[-1])
        )
        lon_selection = np.flatnonzero(
            (source_lons >= lon_edges[0]) & (source_lons <= lon_edges[-1])
        )
        if len(lat_selection) == 0 or len(lon_selection) == 0:
            return _empty_methanesat_observations(species)

        # Rectilinear L3 coordinates form contiguous domain subsets.
        lat_start, lat_stop = lat_selection[0], lat_selection[-1] + 1
        lon_start, lon_stop = lon_selection[0], lon_selection[-1] + 1
        selected_lons = source_lons[lon_start:lon_stop]
        lon_bins = np.searchsorted(lon_edges, selected_lons, side="right") - 1
        valid_lon = (lon_bins >= 0) & (lon_bins < n_lon)

        for row_start in range(lat_start, lat_stop, row_chunk_size):
            row_stop = min(row_start + row_chunk_size, lat_stop)
            selected_lats = source_lats[row_start:row_stop]
            lat_bins = np.searchsorted(lat_edges, selected_lats, side="right") - 1

            key = (slice(row_start, row_stop), slice(lon_start, lon_stop))
            xch4 = values_with_nan(ncfile["xch4"], key)
            times = values_with_nan(ncfile["time"], key)
            surface_pressure = values_with_nan(
                ncfile.groups["apriori_data"]["surface_pressure"], key
            )

            valid = (
                np.isfinite(xch4) & (xch4 > 0)
                & np.isfinite(times) & (times > 0)
                & np.isfinite(surface_pressure) & (surface_pressure > 0)
            )
            if start_seconds is not None:
                valid &= times >= start_seconds
            if end_seconds is not None:
                valid &= times < end_seconds if end_is_date else times <= end_seconds
            valid &= valid_lon[np.newaxis, :]
            valid &= (lat_bins[:, np.newaxis] >= 0)
            valid &= (lat_bins[:, np.newaxis] < n_lat)
            if not np.any(valid):
                continue

            cell_indices = (
                lat_bins[:, np.newaxis] * n_lon + lon_bins[np.newaxis, :]
            )
            flat_cells = np.broadcast_to(cell_indices, valid.shape)[valid]
            counts += np.bincount(flat_cells, minlength=n_cells)
            sums_xch4 += np.bincount(
                flat_cells, weights=xch4[valid], minlength=n_cells
            )
            sums_time += np.bincount(
                flat_cells, weights=times[valid], minlength=n_cells
            )
            sums_surface_pressure += np.bincount(
                flat_cells, weights=surface_pressure[valid], minlength=n_cells
            )
            chunk_latitudes = np.broadcast_to(
                selected_lats[:, np.newaxis], valid.shape
            )
            chunk_longitudes = np.broadcast_to(
                selected_lons[np.newaxis, :], valid.shape
            )
            sums_latitude += np.bincount(
                flat_cells, weights=chunk_latitudes[valid], minlength=n_cells
            )
            sums_longitude += np.bincount(
                flat_cells, weights=chunk_longitudes[valid], minlength=n_cells
            )

    populated = counts > 0
    if not np.any(populated):
        return _empty_methanesat_observations(species)

    mean_xch4 = sums_xch4[populated] / counts[populated]
    mean_time = sums_time[populated] / counts[populated]
    mean_surface_pressure = sums_surface_pressure[populated] / counts[populated]
    mean_latitude = sums_latitude[populated] / counts[populated]
    mean_longitude = sums_longitude[populated] / counts[populated]
    j_gc, i_gc = np.unravel_index(np.flatnonzero(populated), (n_lat, n_lon))

    # Fixed profiles retained from the original MethaneSAT implementation.
    fixed_ak = np.array([
        0.54948103, 0.5463017, 0.5425764, 0.5630824, 0.586732,
        0.61969376, 0.6751325, 0.7441332, 0.7991072, 0.8353182,
        0.87577033, 0.92045873, 0.9373637, 0.97244173, 0.99397856,
        1.0006262, 1.0219611, 1.0428349, 1.0573764,
    ], dtype=np.float64)[::-1]
    fixed_prior = np.array([
        224.05298, 685.95667, 1449.2792, 1765.1627, 1894.5641,
        1929.3362, 1993.0267, 2023.5267, 2023.5922, 2023.6508,
        2024.5969, 2024.6519, 2024.6962, 2024.9379, 2024.9979,
        2025.0558, 2038.752, 2040.5107, 2040.6727,
    ], dtype=np.float64)[::-1]

    tropopause_pressure = np.full_like(mean_surface_pressure, 100.0)
    pressure_edges = gosat_pressure_grid(
        mean_surface_pressure[:, np.newaxis],
        tropopause_pressure[:, np.newaxis],
    )[:, 0, :]
    dry_air = compute_profile_constgrav(
        pressure_edges[:, np.newaxis, :]
    )["dry_air_vcd"][:, 0, :]

    observations = _empty_methanesat_observations(species, len(i_gc))
    observations["iGC"] = i_gc
    observations["jGC"] = j_gc
    observations["lat"] = gc_lats[j_gc]
    observations["lon"] = gc_lons[i_gc]
    observations["lat_sat"] = mean_latitude
    observations["lon_sat"] = mean_longitude
    observations[species] = mean_xch4
    observations["p_sat"] = pressure_edges
    observations["surface_pressure"] = mean_surface_pressure
    observations["dry_air_subcolumns"] = dry_air
    observations["apriori"] = fixed_prior * 1e-9 * dry_air
    observations["avkern"] = fixed_ak
    observations["layer"] = np.arange(19, dtype=np.float32)
    observations["observation_count"] = counts[populated]
    observations["nir_albedo"] = np.nan
    observations["swir_albedo"] = np.nan

    for index, seconds in enumerate(mean_time):
        timestamp = pd.to_datetime(seconds, unit="s", utc=True).tz_localize(None)
        observations["time"][index] = (
            get_strdate(timestamp, time_threshold)
            if time_threshold is not None
            else timestamp.round("60min").strftime("%Y%m%d_%H")
        )
    return observations


def _empty_methanesat_observations(species, size=0):
    """Allocate the standard IMI representation for 19-layer MSAT data."""
    dtype = [
        ("iGC", "i4"), ("jGC", "i4"),
        ("lat_sat", "f4"), ("lon_sat", "f4"),
        (species, "f4"), ("time", "U13"),
        ("p_sat", "f4", (20,)), ("surface_pressure", "f4"),
        ("nir_albedo", "f4"), ("swir_albedo", "f4"),
        ("dry_air_subcolumns", "f4", (19,)),
        ("apriori", "f4", (19,)), ("avkern", "f4", (19,)),
        ("layer", "f4", (19,)), ("observation_count", "f4"),
        ("lat", "f4"), ("lon", "f4"),
    ]
    return np.zeros(size, dtype=dtype)


def gosat_pressure_grid(psurf,ptrop):

    ''' UoL GOSAT Pressure grid 

      ARGS:
        psurf[x,t]: Surface pressure [hPa]
        ptrop[x,t]: Tropopause pressure [hPa]

      RETURNS:
        pedge[x,t,z]: UoL Proxy GOSAT Pressure grid

    '''

    # Allocate output pressure edges
    lmx = 19 ; lmx_t = 13
    xmx = psurf.shape[0]
    tmx = psurf.shape[1]
    pedge = np.zeros((xmx, tmx, lmx+1))

    # Tropospheric sigma levels
    tsig = np.linspace(0,1,lmx_t+1)[::-1]
    for l in range(lmx_t+1):
        pedge[:,:,l] = ptrop + tsig[l]*(psurf-ptrop)

    # Constant stratospheric levels
    p_strat = [8.0000000e1, 5.0000000e1, 1.0000000e1, 1.0000000, 1.0000000e-1]
    for l,pl in zip(range(lmx_t+2,lmx+1),p_strat):
        pedge[:,:,l] = pl

    # Intermediate strat level
    pedge[:,:,lmx_t+1] = 0.5*(pedge[:,:,lmx_t]+pedge[:, :, lmx_t+2])

    return pedge

def compute_profile_constgrav(pedge, g0=9.80665):

    # molecular weights [kg/mol]
    mw_dry = 28.9647e-3  # dry air
    mw_h2o = 18.01528e-3 # water vapor

    outp = {}
    outp['pedge'] = pedge
    zmx = outp['pedge'].shape[-1]-1
    outp['pmid'] = 0.5*(outp['pedge'][:, :, 0:zmx]+outp['pedge'][:, :, 1:])

    outp['air_vcd'] = np.zeros_like(outp['pmid'])
    for iz in range(zmx):
#        outp['air_vcd'][:,:,iz] = (outp['pedge'][:,:,iz]-outp['pedge'][:,:,iz+1])/g0 * (1/(outp['mw'][:,:,iz]))
        outp['air_vcd'][:,:,iz] = (outp['pedge'][:,:,iz]-outp['pedge'][:,:,iz+1])*100/g0 * (1/mw_dry)
    outp['dry_air_vcd'] = outp['air_vcd']#*(1.0-xh2o)

    return outp
