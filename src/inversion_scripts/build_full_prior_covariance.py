from __future__ import annotations

import sys

import numpy as np
import xarray as xr
from utils import get_mean_emissions

"""
Build a full prior error covariance matrix with off-diagonal structure based 
on spatial proximity and emission sector similarity. The resulting covariance 
matrix is saved as a .npy file in order of the state vector IDs. Buffer elements 
are added as additional diagonal entries with no off-diagonal covariance.

Following Balasus et al. (2026), the covariance between two state vector elements 
is calculated as the product of:

1. A spatial decay term based on the haversine distance between the grid cells of 
   the two elements, with a specified length scale.
2. A sector similarity term based on the cosine similarity of the emission sector 
   proportions in the prior for the two grid cells. The sector proportions are 
   calculated from the prior emissions fields on the grid, excluding aggregate 
   fields and soil absorption. The cosine similarity is then calculated between 
   the sector proportion vectors for each grid cell to determine how similar the 
   emission profiles are between the two cells.
3. The diagonal elements of the covariance matrix are set to 1.0, representing the 
   variance of each state vector element. The matrix is then scaled in the 
   inversion step based on the specified prior error standard deviation.
   
Example usage:
python build_full_prior_covariance.py /path/to/StateVector.nc /prior/emissions/dir $length_scale_km $StartDate $EndDate $nbuffer_elements
"""

def select_state_vector_subset(
    state_vector: xr.Dataset, prior: xr.Dataset
) -> xr.DataArray:
    """Subset the state vector onto the prior-emissions grid."""
    subset = state_vector["StateVector"].sel(lat=prior.lat, lon=prior.lon)
    if subset.shape != (prior.sizes["lat"], prior.sizes["lon"]):
        raise ValueError("State vector subset shape does not match prior grid shape.")
    return subset


def get_sector_fields(prior: xr.Dataset) -> list[str]:
    """Return sector-resolved prior-emission fields used for similarity weighting."""
    # only consider emissions fields that are on the prior grid
    # and not already aggregated or excluded sectors
    AGGREGATE_FIELDS = {"EmisCH4_Total", "EmisCH4_Total_ExclSoilAbs"}
    EXCLUDED_SECTORS = {"EmisCH4_SoilAbsorb"}
    fields = []
    for name in prior.data_vars:
        if not name.startswith("EmisCH4_"):
            continue
        if name in AGGREGATE_FIELDS or name in EXCLUDED_SECTORS:
            continue
        if prior[name].dims != ("lat", "lon"):
            continue
        fields.append(name)
    if not fields:
        raise ValueError("No sector emission fields were found on the prior grid.")
    return fields


def build_sector_proportions(prior: xr.Dataset, sector_fields: list[str]) -> np.ndarray:
    """Convert sector emissions at each grid cell into normalized sector fractions."""
    stacked = np.stack(
        [prior[name].values.astype(np.float64) for name in sector_fields], axis=-1
    )
    flat = stacked.reshape(-1, stacked.shape[-1])
    totals = flat.sum(axis=1, keepdims=True)
    proportions = np.divide(flat, totals, out=np.zeros_like(flat), where=totals > 0.0)
    return proportions


def cosine_similarity_matrix(vectors: np.ndarray) -> np.ndarray:
    """Compute pairwise cosine similarity for a set of sector-fraction vectors."""
    norms = np.linalg.norm(vectors, axis=1)
    normalized = np.divide(
        vectors,
        norms[:, None],
        out=np.zeros_like(vectors),
        where=norms[:, None] > 0.0,
    )
    similarity = normalized @ normalized.T
    np.fill_diagonal(similarity, 1.0)
    return np.clip(similarity, 0.0, 1.0)


def haversine_distance_km(latitudes: np.ndarray, longitudes: np.ndarray) -> np.ndarray:
    """Compute the full pairwise great-circle distance matrix in kilometers."""
    lat_rad = np.deg2rad(latitudes)
    lon_rad = np.deg2rad(longitudes)
    dlat = lat_rad[:, None] - lat_rad[None, :]
    dlon = lon_rad[:, None] - lon_rad[None, :]

    sin_dlat = np.sin(dlat / 2.0)
    sin_dlon = np.sin(dlon / 2.0)
    a = sin_dlat**2 + np.cos(lat_rad)[:, None] * np.cos(lat_rad)[None, :] * sin_dlon**2
    a = np.clip(a, 0.0, 1.0)
    return 2.0 * 6371.0 * np.arcsin(np.sqrt(a))


def build_covariance(
    state_vector_subset: xr.DataArray,
    prior: xr.Dataset,
    length_scale_km: float,
    nbuffer_elements: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Build the normalized prior covariance ordered by state-vector ID."""
    state_vector_ids = state_vector_subset.values.reshape(-1)
    valid_mask = np.isfinite(state_vector_ids)
    if not np.any(valid_mask):
        raise ValueError("No ROI state-vector IDs were found on the prior grid.")
    state_vector_ids = state_vector_ids[valid_mask].astype(np.int32)
    if np.unique(state_vector_ids).size != state_vector_ids.size:
        raise ValueError("Duplicate state-vector IDs found on the prior grid.")

    lat_grid, lon_grid = np.meshgrid(prior.lat.values, prior.lon.values, indexing="ij")
    latitudes = lat_grid.reshape(-1)[valid_mask]
    longitudes = lon_grid.reshape(-1)[valid_mask]

    distances = haversine_distance_km(latitudes, longitudes)
    spatial_decay = np.exp(-distances / length_scale_km)

    sector_fields = get_sector_fields(prior)
    proportions = build_sector_proportions(prior, sector_fields)[valid_mask]
    similarity = cosine_similarity_matrix(proportions)

    covariance = spatial_decay * similarity
    np.fill_diagonal(covariance, 1.0)

    order = np.argsort(state_vector_ids)
    covariance = covariance[order][:, order]
    state_vector_ids = state_vector_ids[order]
    state_vector_ids, covariance = append_buffer_diagonal_elements(
        state_vector_ids=state_vector_ids,
        covariance=covariance,
        nbuffer_elements=nbuffer_elements,
    )
    return state_vector_ids, covariance.astype(np.float32)


def append_buffer_diagonal_elements(
    state_vector_ids: np.ndarray,
    covariance: np.ndarray,
    nbuffer_elements: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Append buffer elements as identity-only rows and columns."""
    if nbuffer_elements == 0:
        return state_vector_ids, covariance

    original_size = covariance.shape[0]
    expanded_size = original_size + nbuffer_elements
    expanded = np.zeros((expanded_size, expanded_size), dtype=np.float64)
    expanded[:original_size, :original_size] = covariance
    expanded[original_size:, original_size:] = np.diag(
        np.full(
            nbuffer_elements,
            1.0,
            dtype=np.float64,
        )
    )

    extra_ids = np.arange(
        original_size + 1, expanded_size + 1, dtype=state_vector_ids.dtype
    )
    expanded_ids = np.concatenate([state_vector_ids, extra_ids])
    return expanded_ids, expanded


def main(
    sv_path: str,
    prior_emis_dir: str,
    length_scale_km: float,
    start_date: str,
    end_date: str,
    nbuffer_elements: int,
) -> None:
    """Create and save the normalized prior covariance bundle for an inversion run."""
    # only have off diagonal covariance for elements in the ROI,
    # so set all buffer elements to False
    state_vector = xr.open_dataset(sv_path)
    last_ROI_element = int(state_vector["StateVector"].max().values) - nbuffer_elements
    sv_mask = state_vector["StateVector"] <= last_ROI_element
    state_vector["StateVector"] = state_vector["StateVector"].where(sv_mask, np.nan)

    print(f"Building prior covariance for matrix {start_date}-{end_date}")
    prior = get_mean_emissions(start_date, end_date, prior_emis_dir)
    state_vector_subset = select_state_vector_subset(state_vector, prior)
    state_vector_ids, covariance = build_covariance(
        state_vector_subset=state_vector_subset,
        prior=prior,
        length_scale_km=length_scale_km,
        nbuffer_elements=nbuffer_elements,
    )
    output_path = "prior_norm_error_covariance.npz"
    np.savez(
        output_path,
        covariance=covariance,
        state_vector_ids=state_vector_ids,
    )
    print(
        f"  wrote {output_path} with shape {covariance.shape} "
        f"and diagonal {float(covariance[0, 0]):.3f}"
    )


if __name__ == "__main__":
    sv_path = sys.argv[1]
    prior_emis_dir = sys.argv[2]
    length_scale_km = float(sys.argv[3])
    start_date = sys.argv[4]
    end_date = sys.argv[5]
    nbuffer_elements = int(sys.argv[6])
    main(
        sv_path, prior_emis_dir, length_scale_km, start_date, end_date, nbuffer_elements
    )
