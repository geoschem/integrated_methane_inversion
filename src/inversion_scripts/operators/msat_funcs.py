from GOOPy.parsers import FIXED_MSAT_AVERAGING_KERNEL, FIXED_MSAT_PRIOR_PROFILE
import netCDF4 as nc
import numpy as np
import pandas as pd
import xarray as xr
from src.inversion_scripts.utils import get_strdate
from src.inversion_scripts.operators.superobservation import (
    imi_superobservation_dtype,
)


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


def _read_state_vector_grid(
    state_vector: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read and validate rectilinear state-vector centers and derive edges."""
    with xr.open_dataset(state_vector) as state_vector_ds:
        gc_lats = np.asarray(state_vector_ds["lat"].values, dtype=np.float64)
        gc_lons = np.asarray(state_vector_ds["lon"].values, dtype=np.float64)

    if gc_lats.ndim != 1 or gc_lons.ndim != 1:
        raise ValueError("MethaneSAT averaging requires 1-D state-vector lat/lon")
    if len(gc_lats) < 2 or len(gc_lons) < 2:
        raise ValueError("State-vector lat/lon must each contain at least two cells")
    return gc_lats, gc_lons, coordinate_edges(gc_lats), coordinate_edges(gc_lons)


def _prepare_target_grid(
    state_vector: str | None,
    target_lats: np.ndarray | None,
    target_lons: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Validate target centers and derive edges for MethaneSAT binning."""
    if target_lats is None and target_lons is None:
        if state_vector is None:
            raise ValueError(
                "MethaneSAT averaging requires a state vector or target coordinates"
            )
        return _read_state_vector_grid(state_vector)
    if target_lats is None or target_lons is None:
        raise ValueError("Both target_lats and target_lons must be provided")

    lats = np.asarray(target_lats, dtype=np.float64)
    lons = np.asarray(target_lons, dtype=np.float64)
    if lats.ndim != 1 or lons.ndim != 1:
        raise ValueError("MethaneSAT averaging requires 1-D target lat/lon")
    if len(lats) < 2 or len(lons) < 2:
        raise ValueError("Target lat/lon must each contain at least two cells")
    return lats, lons, coordinate_edges(lats), coordinate_edges(lons)


def _methanesat_time_bounds(
    gc_startdate: str | None,
    gc_enddate: str | None,
) -> tuple[float | None, float | None, bool]:
    """Convert IMI date bounds to epoch seconds using existing end-date rules."""
    def epoch_seconds(value: str | None) -> float | None:
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
    return start_seconds, end_seconds, end_is_date


def _validate_methanesat_file(ncfile: nc.Dataset) -> None:
    """Validate fields required for chunked MethaneSAT L3 averaging."""
    required_root = {"lat", "lon", "time", "xch4"}
    missing_root = required_root.difference(ncfile.variables)
    if missing_root:
        raise ValueError(f"MSAT file is missing variables: {sorted(missing_root)}")
    if "apriori_data" not in ncfile.groups:
        raise ValueError("MSAT file is missing the apriori_data group")
    if "surface_pressure" not in ncfile.groups["apriori_data"].variables:
        raise ValueError("MSAT file is missing apriori_data/surface_pressure")


def _accumulate_methanesat_chunks(
    data_path: str,
    lat_edges: np.ndarray,
    lon_edges: np.ndarray,
    start_seconds: float | None,
    end_seconds: float | None,
    end_is_date: bool,
    row_chunk_size: int,
) -> tuple[np.ndarray, ...] | None:
    """Accumulate jointly valid MSAT fields into flattened state-grid cells."""
    n_lat = len(lat_edges) - 1
    n_lon = len(lon_edges) - 1
    n_cells = n_lat * n_lon
    sums_xch4 = np.zeros(n_cells, dtype=np.float64)
    sums_time = np.zeros(n_cells, dtype=np.float64)
    sums_surface_pressure = np.zeros(n_cells, dtype=np.float64)
    sums_latitude = np.zeros(n_cells, dtype=np.float64)
    sums_longitude = np.zeros(n_cells, dtype=np.float64)
    counts = np.zeros(n_cells, dtype=np.int64)

    with nc.Dataset(data_path, "r") as ncfile:
        _validate_methanesat_file(ncfile)
        source_lats = np.asarray(ncfile["lat"][:], dtype=np.float64)
        source_lons = np.asarray(ncfile["lon"][:], dtype=np.float64)
        lat_selection = np.flatnonzero(
            (source_lats >= lat_edges[0]) & (source_lats <= lat_edges[-1])
        )
        lon_selection = np.flatnonzero(
            (source_lons >= lon_edges[0]) & (source_lons <= lon_edges[-1])
        )
        if len(lat_selection) == 0 or len(lon_selection) == 0:
            return None

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
            valid &= lat_bins[:, np.newaxis] >= 0
            valid &= lat_bins[:, np.newaxis] < n_lat
            if not np.any(valid):
                continue

            cell_indices = lat_bins[:, np.newaxis] * n_lon + lon_bins[np.newaxis, :]
            flat_cells = np.broadcast_to(cell_indices, valid.shape)[valid]
            counts += np.bincount(flat_cells, minlength=n_cells)
            sums_xch4 += np.bincount(flat_cells, weights=xch4[valid], minlength=n_cells)
            sums_time += np.bincount(flat_cells, weights=times[valid], minlength=n_cells)
            sums_surface_pressure += np.bincount(
                flat_cells, weights=surface_pressure[valid], minlength=n_cells
            )
            chunk_latitudes = np.broadcast_to(selected_lats[:, np.newaxis], valid.shape)
            chunk_longitudes = np.broadcast_to(selected_lons[np.newaxis, :], valid.shape)
            sums_latitude += np.bincount(
                flat_cells, weights=chunk_latitudes[valid], minlength=n_cells
            )
            sums_longitude += np.bincount(
                flat_cells, weights=chunk_longitudes[valid], minlength=n_cells
            )

    return (
        sums_xch4, sums_time, sums_surface_pressure,
        sums_latitude, sums_longitude, counts,
    )


def _format_methanesat_observations(
    accumulators: tuple[np.ndarray, ...] | None,
    gc_lats: np.ndarray,
    gc_lons: np.ndarray,
    species: str,
    time_threshold: str | None,
) -> np.ndarray:
    """Convert flattened MSAT accumulators to the standard IMI structure."""
    if accumulators is None:
        return _empty_methanesat_observations(species)
    (
        sums_xch4, sums_time, sums_surface_pressure,
        sums_latitude, sums_longitude, counts,
    ) = accumulators
    populated = counts > 0
    if not np.any(populated):
        return _empty_methanesat_observations(species)

    mean_xch4 = sums_xch4[populated] / counts[populated]
    mean_time = sums_time[populated] / counts[populated]
    mean_surface_pressure = sums_surface_pressure[populated] / counts[populated]
    mean_latitude = sums_latitude[populated] / counts[populated]
    mean_longitude = sums_longitude[populated] / counts[populated]
    j_gc, i_gc = np.unravel_index(
        np.flatnonzero(populated), (len(gc_lats), len(gc_lons))
    )

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
    observations["apriori"] = FIXED_MSAT_PRIOR_PROFILE * 1e-9 * dry_air
    observations["avkern"] = FIXED_MSAT_AVERAGING_KERNEL
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


def average_methanesat_observations(
    data_path: str,
    state_vector: str | None,
    species: str = "CH4",
    time_threshold: str = None,
    gc_startdate: str = None,
    gc_enddate: str = None,
    row_chunk_size: int = 256,
    target_lats: np.ndarray | None = None,
    target_lons: np.ndarray | None = None,
) -> np.ndarray:
    """Average a rectilinear MethaneSAT L3 file onto a rectilinear target grid.

    The input can be much larger than memory. Only coordinate rows intersecting
    the target-grid domain are read, and the two-dimensional retrieval fields
    are accumulated in row chunks. A single joint validity mask is applied to
    XCH4, time, and surface pressure so all averaged fields and counts describe
    the same contributing pixels. If ``target_lats`` and ``target_lons`` are
    supplied, they define the output grid; otherwise the state-vector
    coordinates are retained as a compatibility fallback for preview workflows
    that run before GEOS-Chem output exists.
    """

    gc_lats, gc_lons, lat_edges, lon_edges = _prepare_target_grid(
        state_vector, target_lats, target_lons
    )
    start_seconds, end_seconds, end_is_date = _methanesat_time_bounds(
        gc_startdate, gc_enddate
    )
    accumulators = _accumulate_methanesat_chunks(
        data_path=data_path,
        lat_edges=lat_edges,
        lon_edges=lon_edges,
        start_seconds=start_seconds,
        end_seconds=end_seconds,
        end_is_date=end_is_date,
        row_chunk_size=row_chunk_size,
    )
    return _format_methanesat_observations(
        accumulators=accumulators,
        gc_lats=gc_lats,
        gc_lons=gc_lons,
        species=species,
        time_threshold=time_threshold,
    )


def _empty_methanesat_observations(species: str, size: int =0) -> np.ndarray:
    """Allocate the standard IMI representation for 19-layer MSAT data."""
    return np.zeros(
        size,
        dtype=imi_superobservation_dtype(
            species,
            n_pressure_edges=20,
            n_layers=19,
        ),
    )


def gosat_pressure_grid(
    psurf: np.ndarray,
    ptrop: np.ndarray,
) -> np.ndarray:
    """Construct the 19-layer University of Leicester GOSAT pressure grid.

    Args:
        psurf: Surface pressure in hPa. Its dimensions define the horizontal
            or observation dimensions of the returned array.
        ptrop: Tropopause pressure in hPa with the same shape as ``psurf``.

    Returns:
        Pressure edges in hPa with shape ``psurf.shape + (20,)``. Edges are
        ordered from the surface toward the top of the atmosphere.
    """

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

def compute_profile_constgrav(
    pedge: np.ndarray,
    g0: float = 9.80665,
) -> dict[str, np.ndarray]:
    """Derive layer pressures and dry-air columns using constant gravity.

    Args:
        pedge: Pressure edges in hPa. The final dimension contains vertical
            edges ordered from the surface toward the top of the atmosphere.
        g0: Constant gravitational acceleration in m s-2.

    Returns:
        Dictionary containing ``pedge`` and layer-midpoint pressure ``pmid``
        in hPa, plus ``air_vcd`` and ``dry_air_vcd`` in mol m-2. The layer
        arrays have one fewer element along the final dimension than ``pedge``.
        Water vapor is not currently removed, so ``dry_air_vcd`` equals
        ``air_vcd``.
    """

    # molecular weights [kg/mol]
    mw_dry = 28.9647e-3  # dry air
    mw_h2o = 18.01528e-3 # water vapor

    outp = {}
    outp['pedge'] = pedge
    zmx = outp['pedge'].shape[-1]-1
    outp['pmid'] = 0.5*(outp['pedge'][:, :, 0:zmx]+outp['pedge'][:, :, 1:])

    outp['air_vcd'] = np.zeros_like(outp['pmid'])
    for iz in range(zmx):
        outp['air_vcd'][:,:,iz] = (outp['pedge'][:,:,iz]-outp['pedge'][:,:,iz+1])*100/g0 * (1/mw_dry)
    outp['dry_air_vcd'] = outp['air_vcd']#*(1.0-xh2o)

    return outp
