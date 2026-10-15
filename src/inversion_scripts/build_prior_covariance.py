#!/usr/bin/env python3
"""Build the two-component prior error covariance (Sa) for the IMI.

Selected with PriorCovarianceMethod: aggregate_preserving. Each emission sector adds one block to Sa; the blocks
are summed (in absolute emission-squared units) and split into the two arrays the IMI already expects -- a
correlation matrix with ones on the diagonal, and a per-element relative uncertainty (sigma_scale) -- written
to prior_norm_error_covariance.npz and read back by build_prior_covariance in invert.py.

Anthropogenic sectors get a two-part block: one part where all cells in a region move together, one part where
each cell varies on its own. The parts are scaled so the region's total uncertainty equals a target value. A
region is the whole domain by default, or one country at a time when RespectCountryBorders is set.

Natural sectors get cells correlated by distance (exp(-distance/length)) plus a small fully-correlated floor
from the Saunois global methane budget. Wetlands take their per-cell uncertainty and length from a bundled
model-ensemble field; the other natural sectors use one flat value.

Usage: build_prior_covariance.py StateVectorFile PriorEmisDir config StartDate EndDate nBuffer
"""
import os
import sys

import numpy as np
import xarray as xr

from utils import get_mean_emissions
from config_utils import load_config
from build_full_prior_covariance import get_sector_fields, haversine_distance_km

SECONDS_PER_YEAR = 3.1536e7
GRID_RATIO_DEFAULT = 2.5     # how much a single native cell's uncertainty can exceed the region total
DEFAULT_UNCERTAINTY = 0.30   # target per-sector uncertainty when the config does not set one
WETLAND_SIGMA_FLOOR = 0.15   # smallest per-cell wetland uncertainty
WETLAND_SIGMA_CAP = 1.77     # largest per-cell wetland uncertainty (the South America 3-member ceiling)
GENERIC_LENGTH_KM = 161.0    # correlation length for the non-wetland natural sectors
GENERIC_SIGMA = 0.50         # per-cell uncertainty for the non-wetland natural sectors

ANTHRO_SECTORS = {
    "Livestock", "Rice", "Landfills", "Wastewater", "Coal", "Gas", "Oil", "OilGas",
    "Agriculture", "FossilFuels", "Industry", "Transport", "Other", "BiofuelBurn", "OtherAnth",
}

# For each natural sector: the global emission total (Tg CH4/yr, the denominator for domain_fraction) and the
# relative 1-sigma floor amplitude, both from the Saunois et al. (2020) Global Methane Budget.
SAUNOIS_GLOBAL_TOTALS = {
    "Wetlands": 149.0, "Livestock": 111.0, "Rice": 30.0, "Landfills": 45.0, "Wastewater": 23.0,
    "Coal": 42.0, "Gas": 40.0, "Oil": 40.0, "BiomassBurn": 17.0, "OtherAnth": 10.0,
    "Termites": 9.0, "Seeps": 45.0, "Lakes": 159.0, "Reservoirs": 30.0,
}
SAUNOIS_GLOBAL_SIGMA = {
    "Wetlands": 0.28, "Termites": 0.60, "Seeps": 0.60, "Reservoirs": 0.80,
    "Lakes": 0.80, "BiomassBurn": 0.41,
}


def add_anthropogenic_block(Sa, cells, emis, n_eff, target_uncertainty, grid_ratio):
    """Add one anthropogenic sector block for a region.

    One part has every cell move together; the other lets each cell vary on its own. Both are scaled by the
    amplitude so the region's total relative uncertainty equals target_uncertainty. grid_ratio sets how much
    a single native cell may exceed that total; the per-cell excess shrinks as a cell pools more native cells."""
    cell_excess = (grid_ratio ** 2 - 1.0) / n_eff
    region_excess = float((cell_excess * emis ** 2).sum()) / (emis.sum() ** 2)
    amplitude = target_uncertainty / np.sqrt(1.0 + region_excess)
    Sa[np.ix_(cells, cells)] += (amplitude ** 2) * np.outer(emis, emis)    # all cells move together
    Sa[cells, cells] += (amplitude ** 2) * cell_excess * emis ** 2         # plus an independent per-cell part


def add_natural_block(Sa, cells, emis, sigma, distance_km, length_km, floor_sigma=None, domain_fraction=1.0):
    """Add one natural sector block: per-cell uncertainty correlated by distance (exp(-distance/length)),
    plus an optional fully-correlated floor of size floor_sigma * domain_fraction."""
    cell_error = sigma * emis
    Sa[np.ix_(cells, cells)] += np.outer(cell_error, cell_error) * np.exp(-distance_km / length_km)
    if floor_sigma:
        Sa[np.ix_(cells, cells)] += (floor_sigma * domain_fraction) ** 2 * np.outer(emis, emis)


def decompose(Sa_abs, emis_total, fallback_sigma):
    """Split an absolute emission-squared covariance into a correlation matrix (ones on the diagonal) and a
    per-element relative uncertainty. Elements with no emissions fall back to fallback_sigma, uncorrelated."""
    emis = np.where(emis_total > 0, emis_total, 1.0)
    rel_cov = Sa_abs / np.outer(emis, emis)
    sigma = np.sqrt(np.clip(np.diag(rel_cov), 0.0, None))
    safe_sigma = np.where(sigma > 0, sigma, 1.0)
    corr = rel_cov / np.outer(safe_sigma, safe_sigma)
    np.fill_diagonal(corr, 1.0)
    return corr, np.where(sigma > 0, sigma, fallback_sigma)


def nearest_psd_correlation(corr, floor=1.0e-6):
    """Return the closest valid correlation matrix (positive semidefinite, ones on the diagonal)."""
    eigvals, eigvecs = np.linalg.eigh(0.5 * (corr + corr.T))
    fixed = (eigvecs * np.clip(eigvals, floor, None)) @ eigvecs.T
    scale = np.sqrt(np.clip(np.diag(fixed), floor, None))
    fixed = fixed / scale[:, None] / scale[None, :]
    np.fill_diagonal(fixed, 1.0)
    return fixed


def assign_country_cells(prior):
    """Label every grid cell (offshore included) with the id of its nearest country, matched against the
    Natural Earth layer geopandas ships with -- the same source the clustering uses for GroupByCountry."""
    import warnings
    import geopandas as gpd
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    countries = countries[["name", "geometry"]].dropna(subset=["name"]).dissolve(by="name", as_index=False)
    lon2d, lat2d = np.meshgrid(prior.lon.values, prior.lat.values)
    points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(lon2d.ravel(), lat2d.ravel()), crs="EPSG:4326")
    metric_crs = "ESRI:54030"   # Robinson, a distance-preserving projection, as in the clustering
    nearest = gpd.sjoin_nearest(points.to_crs(metric_crs), countries.to_crs(metric_crs))
    names = nearest[~nearest.index.duplicated(keep="first")].sort_index()["name"].to_numpy()
    name_to_id = {name: i for i, name in enumerate(sorted(dict.fromkeys(names)), start=1)}
    cell_ids = np.array([name_to_id[name] for name in names], int)
    return cell_ids.reshape(prior.lat.size, prior.lon.size)


def load_wetland_sigma(path, elem_lat, elem_lon, month):
    """Read per-cell wetland uncertainty and correlation length from the bundled field (monthly, global),
    sampled at the element centroids and clipped to [WETLAND_SIGMA_FLOOR, WETLAND_SIGMA_CAP]."""
    field = np.load(path)
    rel_sigma = field["rel"][month - 1]
    lat_idx = np.clip(np.searchsorted(field["lat"], elem_lat), 0, field["lat"].size - 1)
    lon_idx = np.clip(np.searchsorted(field["lon"], elem_lon), 0, field["lon"].size - 1)
    sigma = np.clip(np.nan_to_num(rel_sigma[lat_idx, lon_idx]), WETLAND_SIGMA_FLOOR, WETLAND_SIGMA_CAP)
    return sigma, float(field["L_km"])


def build_elements(sv, prior, sector_fields, roi_ids, country_cells):
    """For each state-vector element, gather: centroid (elem_lat, elem_lon), its most common country id, the
    emission per sector (emis[sector], kg/s), and the effective number of native cells it pools
    (n_eff[sector] = (sum emis)^2 / sum(emis^2); 1 for a native-resolution element)."""
    labels = sv.values
    lat2d, lon2d = np.meshgrid(prior.lat.values, prior.lon.values, indexing="ij")
    area = prior["AREA"].values
    row_of_label = {int(label): i for i, label in enumerate(roi_ids)}
    n = len(roi_ids)
    elem_lat, elem_lon, country = np.zeros(n), np.zeros(n), np.zeros(n, int)
    sectors = [name.replace("EmisCH4_", "") for name in sector_fields]
    emis = {s: np.zeros(n) for s in sectors}
    n_eff = {s: np.ones(n) for s in sectors}
    emis_field = {s: np.nan_to_num(np.asarray(prior[name].values, float)) * area
                  for s, name in zip(sectors, sector_fields)}
    for label in roi_ids:
        mask = labels == int(label)
        if not np.any(mask):
            continue
        i = row_of_label[int(label)]
        elem_lat[i], elem_lon[i] = float(lat2d[mask].mean()), float(lon2d[mask].mean())
        cell_countries = country_cells[mask][country_cells[mask] > 0]
        country[i] = int(np.bincount(cell_countries).argmax()) if cell_countries.size else 0
        for s in sectors:
            cell_emis = emis_field[s][mask]
            cell_emis = cell_emis[cell_emis > 0]
            emis[s][i] = float(cell_emis.sum())
            if cell_emis.size:
                n_eff[s][i] = float(cell_emis.sum() ** 2 / (cell_emis ** 2).sum())
    return elem_lat, elem_lon, country, emis, n_eff


def main(sv_path, prior_emis_dir, config_path, start_date, end_date, nbuffer_elements, install_path=None):
    config = load_config(config_path)
    nbuffer = int(nbuffer_elements)
    respect_borders = str(config.get("RespectCountryBorders", False)).strip().lower() in ("true", "1", "yes")
    prior_sigma = float((config.get("PriorError") or [0.5])[0])
    grid_ratio = float(config.get("PriorErrorGridRatio", GRID_RATIO_DEFAULT))
    target_uncertainty = float(config.get("DefaultDomainPriorUncertainty", DEFAULT_UNCERTAINTY))

    prior = get_mean_emissions(start_date, end_date, prior_emis_dir)
    sector_fields = get_sector_fields(prior)
    state_vector = xr.open_dataset(sv_path)
    if "time" in state_vector.dims:
        state_vector = state_vector.isel(time=0)
    sv = state_vector["StateVector"].sel(lat=prior.lat, lon=prior.lon)

    # region-of-interest elements: finite, positive, and excluding the trailing buffer elements
    in_roi = np.isfinite(sv.values) & (sv.values > 0)
    if nbuffer > 0:
        in_roi &= sv.values <= int(np.nanmax(sv.values)) - nbuffer
    roi_ids = np.unique(sv.values[in_roi].astype(np.int32))
    n = len(roi_ids)

    # a "region" is the whole domain by default, or one country at a time when respecting borders
    if respect_borders:
        country_cells = assign_country_cells(prior)
    else:
        country_cells = np.ones((prior.lat.size, prior.lon.size), int)
    elem_lat, elem_lon, country, emis, n_eff = build_elements(sv, prior, sector_fields, roi_ids, country_cells)
    sectors = list(emis.keys())
    emis_total = np.sum([emis[s] for s in sectors], axis=0)

    Sa = np.zeros((n, n))
    n_anthro_blocks, regions, n_natural_blocks = 0, set(), 0

    # one anthropogenic block per (region, sector). Respecting borders only stops the correlation from
    # crossing country lines; each sector's amplitude is raised by 1 / sqrt(sum share^2) (share = each
    # country's share of the sector's emissions) so the whole-domain uncertainty still equals the target.
    for s in [s for s in sectors if s in ANTHRO_SECTORS]:
        country_ids = np.unique(country[emis[s] > 0])
        country_emis = np.array([emis[s][(country == cid) & (emis[s] > 0)].sum() for cid in country_ids])
        if country_emis.sum() > 0:
            block_uncertainty = target_uncertainty / np.sqrt(((country_emis / country_emis.sum()) ** 2).sum())
        else:
            block_uncertainty = target_uncertainty
        for cid in country_ids:
            cells = np.where((country == cid) & (emis[s] > 0))[0]
            if cells.size:
                add_anthropogenic_block(Sa, cells, emis[s][cells], n_eff[s][cells], block_uncertainty, grid_ratio)
                n_anthro_blocks += 1
                regions.add(int(cid))

    # one natural block per sector (wetland uncertainty/length from the bundled field, else a flat value).
    # The floor is floor_sigma * domain_fraction, domain_fraction = the domain's clipped share of the
    # sector's Saunois global total.
    distance_km = haversine_distance_km(elem_lat, elem_lon)
    # the bundled wetland field lives in the install tree's resources/. At runtime run_inversion.sh runs a
    # copy of this script from the run dir and passes install_path; __file__ is only beside resources/ in-repo.
    resources_root = install_path or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
    wetland_path = os.path.join(resources_root, "resources", "wetlands", "wetland_ensemble.npz")
    month = int(str(start_date)[4:6])
    for s in [s for s in sectors if s not in ANTHRO_SECTORS]:
        cells = np.where(emis[s] > 0)[0]
        if not cells.size:
            continue
        floor_sigma = SAUNOIS_GLOBAL_SIGMA.get(s)
        global_total = SAUNOIS_GLOBAL_TOTALS.get(s)
        domain_tg = emis[s][cells].sum() * SECONDS_PER_YEAR * 1e-9
        domain_fraction = min(domain_tg / global_total, 1.0) if global_total else 1.0
        if s == "Wetlands":
            sigma, length_km = load_wetland_sigma(wetland_path, elem_lat[cells], elem_lon[cells], month)
        else:
            sigma, length_km = np.full(cells.size, GENERIC_SIGMA), GENERIC_LENGTH_KM
        add_natural_block(Sa, cells, emis[s][cells], sigma, distance_km[np.ix_(cells, cells)], length_km,
                          floor_sigma=floor_sigma, domain_fraction=domain_fraction)
        n_natural_blocks += 1

    corr, sigma = decompose(Sa, emis_total, prior_sigma)
    corr = nearest_psd_correlation(corr)

    # append the buffer elements as an uncorrelated block with unit correlation
    ids = roi_ids
    if nbuffer:
        padded = np.zeros((n + nbuffer, n + nbuffer))
        padded[:n, :n], padded[n:, n:] = corr, np.eye(nbuffer)
        corr = padded
        ids = np.concatenate([roi_ids, np.arange(n + 1, n + nbuffer + 1, dtype=roi_ids.dtype)])
    sigma_scale = np.concatenate([sigma / prior_sigma, np.ones(nbuffer)]).astype(np.float32)
    np.savez("prior_norm_error_covariance.npz",
             covariance=corr.astype(np.float32), state_vector_ids=ids.astype(np.int32), sigma_scale=sigma_scale)

    region_label = "per-country" if respect_borders else "domain-wide"
    print(f"aggregate_preserving prior covariance ({region_label}): {corr.shape[0]} elements; "
          f"{n_anthro_blocks} anthro blocks across {len(regions)} region(s); {n_natural_blocks} natural blocks.")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], int(sys.argv[6]),
         sys.argv[7] if len(sys.argv) > 7 else None)
