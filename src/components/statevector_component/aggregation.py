#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import yaml
import xarray as xr
import numpy as np
import warnings
import functools

from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.neighbors import NearestNeighbors

from src.inversion_scripts.point_sources import get_point_source_coordinates
from src.inversion_scripts.imi_preview import estimate_averaging_kernel, map_sensitivities_to_sv
from src.inversion_scripts.classify_TROPOMI_obs_to_CSgrids import latlon_to_cartesian, build_kdtree

# Always flush prints in batch jobs
print = functools.partial(print, flush=True)


# -----------------------------------------------------------------------------
# Grid helpers
# -----------------------------------------------------------------------------
def precompute_flat_lonlat(config, sv_ds):
    """
    Flatten grid center lon/lat once and reuse everywhere.

    Returns:
        lon_all     1D float : flattened lon (wrapped to [-180, 180] for GCHP)
        lat_all     1D float : flattened lat
        grid_shape  tuple    : original grid shape for reshape
    """
    if config["UseGCHP"]:
        lats = sv_ds["lats"].values
        lons = sv_ds["lons"].values
        grid_shape = lats.shape  # (nf, Ydim, Xdim)
        lat_all = lats.reshape(-1)
        lon_all = lons.reshape(-1).copy()
        lon_all[lon_all > 180] -= 360
    else:
        lat = sv_ds["lat"].values
        lon = sv_ds["lon"].values
        LON, LAT = np.meshgrid(lon, lat)
        grid_shape = LON.shape
        lon_all = LON.reshape(-1)
        lat_all = LAT.reshape(-1)
    return lon_all, lat_all, grid_shape


def estimate_gridstep_xyz_from_roi(lon_flat, lat_flat, roi_idx, sample_n=4000, random_state=0):
    """
    Estimate a typical nearest-neighbor spacing in xyz chord distance (unit sphere).

    Uses a random subsample of ROI cells and k=2 NN (self + nearest neighbor).
    The returned distance is used to scale xyz so that ~1 grid step ≈ 1 unit.
    """
    roi_idx = np.asarray(roi_idx, dtype=np.int64)
    if roi_idx.size < 2:
        return 1.0

    rng = np.random.default_rng(random_state)
    n = int(min(sample_n, roi_idx.size))
    samp = roi_idx if n == roi_idx.size else rng.choice(roi_idx, size=n, replace=False)

    lon = np.asarray(lon_flat[samp], dtype=float).copy()
    lat = np.asarray(lat_flat[samp], dtype=float)

    xyz = latlon_to_cartesian(lat, lon)  # (n, 3)

    nn = NearestNeighbors(n_neighbors=2, algorithm="auto")
    nn.fit(xyz)
    dists, _ = nn.kneighbors(xyz, return_distance=True)

    step = np.nanmedian(dists[:, 1])  # nearest neighbor (exclude self)
    if not np.isfinite(step) or step <= 0:
        step = 1.0
    return float(step)


# -----------------------------------------------------------------------------
# Country mask sampling
# -----------------------------------------------------------------------------
def sample_country_mask_nearest(lat, lon, mask, res_deg=0.1):
    """Nearest-neighbor sampling on a regular lat/lon mask."""
    lat = np.asarray(lat)
    lon = np.asarray(lon).copy()
    lon[lon > 180] -= 360

    j = np.rint((lat - (-90 + res_deg / 2)) / res_deg).astype(np.int64)
    i = np.rint((lon - (-180 + res_deg / 2)) / res_deg).astype(np.int64)

    j = np.clip(j, 0, mask.shape[0] - 1)
    i = np.clip(i, 0, mask.shape[1] - 1)
    return mask[j, i]


def sample_country_mask_majority(lat, lon, mask, ocean_id=-1, res_deg=0.1):
    """
    Majority vote among 5 samples (center + 4 nudges).
    If any land exists, ocean votes are ignored.
    """
    nudge = 0.8 * res_deg
    s0 = sample_country_mask_nearest(lat, lon, mask, res_deg)
    s1 = sample_country_mask_nearest(lat + nudge, lon, mask, res_deg)
    s2 = sample_country_mask_nearest(lat - nudge, lon, mask, res_deg)
    s3 = sample_country_mask_nearest(lat, lon + nudge, mask, res_deg)
    s4 = sample_country_mask_nearest(lat, lon - nudge, mask, res_deg)

    samples = np.stack([s0, s1, s2, s3, s4], axis=0)  # (5, N)
    is_land = samples != ocean_id
    has_land = is_land.any(axis=0)

    eq = (samples[:, None, :] == samples[None, :, :])
    eq &= is_land[:, None, :]
    eq &= is_land[None, :, :]
    counts = eq.sum(axis=1)

    best = counts.argmax(axis=0)
    out = samples[best, np.arange(samples.shape[1])]
    out[~has_land] = ocean_id
    return out


def assign_country_ids_valid(
    config,
    sv_ds,
    valid_indices,
    country_mask_ds,
    mask_var="country_id",
    ocean_id=-1,
    res_deg=0.1,
    majority_vote=True,
    lats_flat=None,
    lons_flat=None,
):
    """
    Assign country_id for selected grid cells.

    valid_indices can be a boolean mask on flattened grid or flattened indices.
    """
    valid_indices = np.asarray(valid_indices)
    idx = np.flatnonzero(valid_indices) if valid_indices.dtype == bool else valid_indices.astype(np.int64)

    mask = country_mask_ds[mask_var].values

    if (lats_flat is None) or (lons_flat is None):
        if config["UseGCHP"]:
            lats_flat = sv_ds["lats"].values.reshape(-1)
            lons_flat = sv_ds["lons"].values.reshape(-1)
        else:
            lon = sv_ds["lon"].values
            lat = sv_ds["lat"].values
            lons, lats = np.meshgrid(lon, lat)
            lats_flat = lats.reshape(-1)
            lons_flat = lons.reshape(-1)

    lat_v = lats_flat[idx]
    lon_v = lons_flat[idx]

    if majority_vote:
        return sample_country_mask_majority(lat_v, lon_v, mask, ocean_id=ocean_id, res_deg=res_deg)
    return sample_country_mask_nearest(lat_v, lon_v, mask, res_deg=res_deg)


# -----------------------------------------------------------------------------
# Grouped KMeans utilities (exact total cluster count across groups)
# -----------------------------------------------------------------------------
def allocate_k_per_group(group_ids, total_k, min_k=1):
    """
    Allocate exactly total_k clusters across groups.

    Guarantees:
      - sum(k) == total_k
      - 0 <= k_i <= count_i
      - min_k applied only when feasible (otherwise some groups get 0)
    """
    ids, counts = np.unique(group_ids, return_counts=True)
    caps = counts.astype(np.int64)
    total_k = int(total_k)

    if total_k <= 0:
        return ids, np.zeros_like(caps, dtype=int)

    n_points = int(caps.sum())
    if total_k > n_points:
        raise ValueError(f"total_k={total_k} exceeds number of points={n_points}")

    n_groups = ids.size
    floor_k = np.zeros(n_groups, dtype=int) if total_k < n_groups * min_k else np.full(n_groups, int(min_k), dtype=int)

    frac = caps / caps.sum()
    k = np.floor(frac * total_k).astype(int)

    k = np.maximum(k, floor_k)
    k = np.minimum(k, caps)

    diff = int(total_k - k.sum())
    rema = frac * total_k - np.floor(frac * total_k)

    if diff > 0:
        while diff > 0:
            can = np.where(k < caps)[0]
            if can.size == 0:
                break
            j = can[np.argmax(rema[can])]
            k[j] += 1
            diff -= 1
    elif diff < 0:
        order = np.argsort(-k)
        for j in order:
            if diff == 0:
                break
            if k[j] > floor_k[j]:
                k[j] -= 1
                diff += 1

    diff = int(total_k - k.sum())
    if diff > 0:
        can = np.where(k < caps)[0]
        for j in can[np.argsort(-(caps[can] - k[can]))]:
            if diff == 0:
                break
            add = min(int(caps[j] - k[j]), diff)
            k[j] += add
            diff -= add
    elif diff < 0:
        order = np.argsort(-k)
        for j in order:
            if diff == 0:
                break
            if k[j] > floor_k[j]:
                k[j] -= 1
                diff += 1

    assert int(k.sum()) == total_k
    assert np.all(k >= 0)
    assert np.all(k <= caps)
    return ids, k.astype(int)


def kmeans_by_group(features, group_ids, num_clusters, mini_batch=True, random_state=0, min_k=1):
    """Run KMeans per group and merge labels (0-based)."""
    features = np.asarray(features)
    group_ids = np.asarray(group_ids)

    ids, k_per = allocate_k_per_group(group_ids, num_clusters, min_k=min_k)

    labels = np.full(group_ids.shape[0], -1, dtype=int)
    next_label = 0

    for gid, k_i in zip(ids, k_per):
        idx = np.where(group_ids == gid)[0]
        n = idx.size
        if n == 0 or k_i == 0:
            continue

        k_i = int(min(k_i, n))
        if k_i == 1:
            labels[idx] = next_label
            next_label += 1
            continue

        km = MiniBatchKMeans(n_clusters=k_i, random_state=random_state) if mini_batch else KMeans(
            n_clusters=k_i, random_state=random_state
        )
        labels[idx] = km.fit_predict(features[idx]) + next_label
        next_label += k_i

    assert next_label == int(num_clusters)
    return labels


# -----------------------------------------------------------------------------
# Cluster-size cap: split oversized clusters (may increase label count)
# -----------------------------------------------------------------------------
def split_oversized_clusters(features, labels, max_cluster_size, mini_batch=True, random_state=0):
    """
    Split clusters with size > max_cluster_size by reclustering within that cluster.
    labels are 0-based; returns updated 0-based labels (may increase label count).
    """
    labels = np.asarray(labels, dtype=np.int64)
    if labels.size == 0:
        return labels

    cap = int(max_cluster_size)
    if cap <= 0:
        return labels

    # --- fast path: if nothing is oversized, do nothing ---
    # (one bincount; avoids entering the split loop and avoids KMeans calls)
    counts = np.bincount(labels, minlength=int(labels.max()) + 1)
    if counts.size == 0 or int(counts.max()) <= cap:
        return labels

    next_label = int(labels.max()) + 1

    while True:
        # reuse counts from above on first iteration; recompute after splits
        oversized = np.where(counts > cap)[0]
        if oversized.size == 0:
            break

        c = int(oversized[0])
        idx = np.where(labels == c)[0]
        n = int(idx.size)

        k_split = int((n + cap - 1) // cap)  # ceil(n/cap)
        if k_split <= 1:
            # shouldn't happen for oversized, but keep loop safe
            counts[c] = n
            continue

        km = (MiniBatchKMeans(n_clusters=k_split, random_state=random_state)
              if mini_batch else KMeans(n_clusters=k_split, random_state=random_state))
        sub = km.fit_predict(features[idx])

        # Keep sub==0 as label c; others become new labels
        counts[c] = int(np.sum(sub == 0))
        for s in range(1, k_split):
            sel = idx[sub == s]
            if sel.size == 0:
                continue
            labels[sel] = next_label
            # extend counts array if needed
            if next_label >= counts.size:
                counts = np.pad(counts, (0, next_label - counts.size + 1), constant_values=0)
            counts[next_label] = int(sel.size)
            next_label += 1

    return labels


# -----------------------------------------------------------------------------
# Light island-splitting pass (removes obvious disconnected islands)
# -----------------------------------------------------------------------------
class _UnionFind:
    def __init__(self, n):
        self.p = np.arange(n, dtype=np.int64)
        self.r = np.zeros(n, dtype=np.int8)

    def find(self, a):
        p = self.p
        while p[a] != a:
            p[a] = p[p[a]]
            a = p[a]
        return a

    def union(self, a, b):
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        r = self.r
        p = self.p
        if r[ra] < r[rb]:
            p[ra] = rb
        elif r[ra] > r[rb]:
            p[rb] = ra
        else:
            p[rb] = ra
            r[ra] += 1


def split_disconnected_components(lat, lon, labels):
    """
    Split each label into connected components using cheap adjacency candidates.

    This is a lightweight cleanup step. It does not build a full neighbor graph;
    it mainly removes very obvious disconnected islands created by KMeans.
    """
    labels = np.asarray(labels, dtype=np.int64)
    n = labels.size
    if n == 0:
        return labels

    xyz = latlon_to_cartesian(np.asarray(lat, float), np.asarray(lon, float))
    uf = _UnionFind(n)

    # Link adjacent pairs in sorted order along each axis (cheap candidate edges)
    for dim in range(3):
        order = np.argsort(xyz[:, dim], kind="mergesort")
        a = order[:-1]
        b = order[1:]
        same = labels[a] == labels[b]
        aa = a[same]
        bb = b[same]
        for i in range(aa.size):
            uf.union(int(aa[i]), int(bb[i]))

    roots = np.fromiter((uf.find(i) for i in range(n)), dtype=np.int64, count=n)
    pairs = np.stack([labels, roots], axis=1)

    uniq = np.unique(pairs, axis=0)
    uniq = uniq[np.lexsort((uniq[:, 1], uniq[:, 0]))]

    key_to_new = {}
    next_lab = 0
    for lab, root in uniq:
        key_to_new[(int(lab), int(root))] = next_lab
        next_lab += 1

    out = np.empty_like(labels)
    for i in range(n):
        out[i] = key_to_new[(int(labels[i]), int(roots[i]))]
    return out


# -----------------------------------------------------------------------------
# Subset clustering: spatial geometry + DOFS-scaled sensitivity
# -----------------------------------------------------------------------------
def cluster_data_kmeans(
    config,
    sensi_flat,
    sv_ds,
    num_clusters,
    mini_batch=False,
    cluster_by_country=False,
    _country_mask_cache=None,
    lon_all=None,
    lat_all=None,
    subset_idx=None,
    grid_shape=None,
    max_cluster_size=None,
    gridstep_xyz=None,
):
    """
    Cluster subset_idx cells into ~num_clusters clusters.

    Output:
        full-grid labels (0 outside subset_idx), labels are 1-based.

    Feature design:
      - geometry: xyz chord distance (unit sphere), scaled so ~1 grid step ≈ 1 unit
      - sensitivity: monotone compression of (Z / ClusteringThreshold), scaled to ~O(1)
    """
    if subset_idx is None:
        raise ValueError("subset_idx must be provided for clustering")
    if lon_all is None or lat_all is None or grid_shape is None:
        raise ValueError("lon_all, lat_all, and grid_shape must be provided")

    subset_idx = np.asarray(subset_idx, dtype=np.int64)
    out = np.zeros(sensi_flat.size, dtype=np.int32)

    n_sel = int(subset_idx.size)
    K_req = int(num_clusters)
    if n_sel == 0 or K_req <= 0:
        return out.reshape(grid_shape)

    # Drop NaNs inside subset
    Z0 = sensi_flat[subset_idx]
    finite = np.isfinite(Z0)
    if not np.any(finite):
        return out.reshape(grid_shape)

    idx = subset_idx[finite]
    Z = Z0[finite]

    lon = np.asarray(lon_all[idx], dtype=float)
    lat = np.asarray(lat_all[idx], dtype=float)

    xyz = latlon_to_cartesian(lat, lon)  # (N, 3)
    step = float(gridstep_xyz) if (gridstep_xyz is not None) else 1.0
    if (not np.isfinite(step)) or step <= 0:
        step = 1.0
    xyz_feat = xyz / step

    thr = float(config.get("ClusteringThreshold", 1.0))
    thr = max(thr, 1e-12)

    # Monotone compression keeps “low info merges easily” behavior stable
    z_raw = np.log1p(np.maximum(Z, 0.0) / thr)

    # Put sensitivity feature on ~O(1) scale (avoid huge numeric imbalance)
    z_scale = np.nanmedian(z_raw)
    if (not np.isfinite(z_scale)) or z_scale <= 0:
        z_scale = 1.0
    s_feat = z_raw / z_scale

    features = np.column_stack((xyz_feat, s_feat))

    # Ensure enough clusters to respect MaxClusterSize (ROI-only; buffer handled elsewhere)
    K = int(K_req)
    if max_cluster_size is not None and int(max_cluster_size) > 0:
        cap = int(max_cluster_size)
        K_min = int((idx.size + cap - 1) // cap)  # ceil(N/cap)
        K = max(K, K_min)

    K = min(K, int(idx.size))
    if K <= 0:
        return out.reshape(grid_shape)

    if cluster_by_country:
        # Cache country mask across calls
        if _country_mask_cache is not None and "ds" in _country_mask_cache:
            country_mask_ds = _country_mask_cache["ds"]
        else:
            country_mask_ds = xr.open_dataset(config["CountryMaskPath"])
            if _country_mask_cache is not None:
                _country_mask_cache["ds"] = country_mask_ds

        country_id = assign_country_ids_valid(
            config,
            sv_ds,
            idx,
            country_mask_ds,
            lats_flat=lat_all,
            lons_flat=lon_all,
        )
        labels0 = kmeans_by_group(features, country_id, K, mini_batch, random_state=0, min_k=1)
    else:
        km = MiniBatchKMeans(n_clusters=K, random_state=0) if mini_batch else KMeans(n_clusters=K, random_state=0)
        labels0 = km.fit_predict(features)

    # Cap cluster size (may increase label count)
    if max_cluster_size is not None and int(max_cluster_size) > 0:
        labels0 = split_oversized_clusters(features, labels0, int(max_cluster_size), mini_batch=mini_batch, random_state=0)

    # Light island splitting (may increase label count)
    labels0 = split_disconnected_components(lat, lon, labels0)

    out[idx] = labels0.astype(np.int32) + 1
    return out.reshape(grid_shape)


# -----------------------------------------------------------------------------
# Rank clusters by summed sensitivity and apply threshold
# -----------------------------------------------------------------------------
def get_highest_labels_threshold(labels, sensitivities, threshold):
    """Return cluster IDs whose total sensitivity meets threshold (descending)."""
    lab = np.asarray(labels).ravel()
    sen = np.asarray(sensitivities).ravel()

    m = (lab >= 1) & np.isfinite(sen)
    if not np.any(m):
        return (np.array([], dtype=int), 0, np.array([], dtype=int), np.array([], dtype=float))

    lab_m = lab[m].astype(np.int64)
    sen_m = sen[m].astype(np.float64)

    max_label = int(lab_m.max())
    total_sensi = np.bincount(lab_m, weights=sen_m, minlength=max_label + 1)  # label 0 unused
    counts = np.bincount(lab_m, minlength=max_label + 1)

    ids = np.arange(1, max_label + 1, dtype=np.int64)
    order = np.argsort(total_sensi[1:])[::-1]

    n_clusters = ids[order]
    n_sensis = total_sensi[1:][order]

    keep = n_sensis >= threshold
    n = int(np.count_nonzero(keep))
    return (n_clusters[:n], n, counts[1:][order][:n], np.round(n_sensis[:n], 2))


# -----------------------------------------------------------------------------
# Max cluster size presets (or config override)
# -----------------------------------------------------------------------------
def get_max_cluster_size(config, sensitivities, desired_element_num):
    """Choose MaxClusterSize from presets unless overridden by config."""
    if config["UseGCHP"]:
        preset = {"720": 128, "360": 64, "180": 32, "90": 16, "48": 8, "24": 4}.get(str(config.get("CS_RES")), 64)
    else:
        preset = {
            "0.125x0.15625": 128,
            "0.25x0.3125": 64,
            "0.5x0.625": 32,
            "2.0x2.5": 16,
            "4.0x5.0": 8,
        }.get(str(config.get("Res")), 64)

    max_cluster_size = int(config["MaxClusterSize"]) if "MaxClusterSize" in config else int(preset)

    # Same feasibility check as your original logic
    background_elements_needed = np.ceil(len(sensitivities) / max_cluster_size)
    if background_elements_needed > desired_element_num:
        raise Exception(
            "Error: too few clusters to have background state vector elements\n"
            + f"aggregating {max_cluster_size} native resolution elements.\n"
            + "More state vector elements recommended. Alternatively, raise the\n"
            + '"MaxClusterSize" allowed for the region of interest.'
        )

    print(f"MaxClusterSize set to: {max_cluster_size} elements in a cluster")
    return max_cluster_size


# -----------------------------------------------------------------------------
# Force native-resolution elements at point-source locations
# -----------------------------------------------------------------------------
def force_native_res_pixels(config, clusters_ds, sensitivities):
    """Snap point sources to grid cells and raise their sensitivities above threshold."""
    dofs_max = float(config["ClusteringThreshold"]) + 0.1 if "ClusteringThreshold" in config else 1.1

    raw_coords = get_point_source_coordinates(config)
    if len(raw_coords) == 0:
        print("No ForcedNativeResolutionElements or PointSourceDatasets specified in config file.")
        return sensitivities

    if config["UseGCHP"]:
        grid_lats = clusters_ds["lats"].values
        grid_lons = clusters_ds["lons"].values
        lon_min, lon_max = -180, 180
        lat_min, lat_max = -90, 90
    else:
        lat = clusters_ds["lat"].values
        lon = clusters_ds["lon"].values
        grid_lons, grid_lats = np.meshgrid(lon, lat)
        delta_lon = float(np.median(np.abs(np.diff(lon))))
        delta_lat = float(np.median(np.abs(np.diff(lat))))
        lon_min = max(float(np.nanmin(grid_lons) - delta_lon / 2), -180.0)
        lon_max = min(float(np.nanmax(grid_lons) + delta_lon / 2), 180.0)
        lat_min = max(float(np.nanmin(grid_lats) - delta_lat / 2), -90.0)
        lat_max = min(float(np.nanmax(grid_lats) + delta_lat / 2), 90.0)

    pts = np.asarray(raw_coords, dtype=float)
    pts_lat = pts[:, 0]
    pts_lon = pts[:, 1].copy()
    pts_lon[pts_lon > 180] -= 360

    roi_mask = (pts_lat >= lat_min) & (pts_lat <= lat_max) & (pts_lon >= lon_min) & (pts_lon <= lon_max)
    raw_coords = np.stack([pts_lat[roi_mask], pts_lon[roi_mask]], axis=1).tolist()

    if len(raw_coords) == 0:
        print("No forced point sources found in the region of interest.")
        return sensitivities

    print(f"Found {len(raw_coords)} point sources at the grid resolution in the region of interest.")

    kdtree, _ = build_kdtree(grid_lats, grid_lons)
    pts = np.asarray(raw_coords, dtype=float)
    q_cart = latlon_to_cartesian(pts[:, 0], pts[:, 1])

    _, idx_flat = kdtree.query(q_cart, k=1)
    idx_flat = np.unique(np.asarray(idx_flat).reshape(-1))

    print(f"{len(raw_coords)} sources → {len(idx_flat)} grid cells")

    if "NumberOfElements" in config:
        max_n = int(config["NumberOfElements"])
        if len(idx_flat) > max_n:
            idx_flat = idx_flat[:max_n]

    clusters_flat = clusters_ds["StateVector"].values.reshape(-1)
    cluster_ids = clusters_flat[idx_flat]

    valid = np.isfinite(cluster_ids) & (cluster_ids >= 1)
    forced_sv_idx = np.unique(cluster_ids[valid].astype(int) - 1)

    sensitivities[forced_sv_idx] = dofs_max
    return sensitivities


# -----------------------------------------------------------------------------
# Core aggregation driver (ROI clustering + buffer reattachment)
# -----------------------------------------------------------------------------
def update_sv_clusters(config, flat_sensi, orig_sv_ds):
    """
    Create a reduced state vector.

    ROI behavior:
      - Target ROI labels = NumberOfElements - nBufferClusters
      - MaxClusterSize is enforced inside clustering (K floor + split pass)
      - Spatial coherence improved by xyz geometry + island splitting pass
    """
    if config["ClusteringMethod"] == "kmeans":
        mini_batch = False
    elif config["ClusteringMethod"] == "mini-batch-kmeans":
        mini_batch = True
    else:
        raise Exception("Error: Invalid Clustering Method. Use 'kmeans' or 'mini-batch-kmeans'.")

    orig_sv = orig_sv_ds["StateVector"]
    orig_sv_np = orig_sv.values

    desired_num_labels = int(config["NumberOfElements"] - config["nBufferClusters"])
    last_ROI_element = int(np.nanmax(orig_sv_np) - config["nBufferClusters"])

    # DOFS threshold: user-provided or estimated
    if "ClusteringThreshold" in config:
        dofs_threshold = float(config["ClusteringThreshold"])
    else:
        dofs_threshold = float(np.sum(flat_sensi) / desired_num_labels)
        if dofs_threshold > 1:
            print(f"Estimated dofs per element too high ({dofs_threshold}), resetting ClusteringThreshold to 1")
            dofs_threshold = 1.0
    print(f"Target DOFS per cluster (ClusteringThreshold): {dofs_threshold}")

    max_cluster_size = get_max_cluster_size(config, flat_sensi, desired_num_labels)

    if "GroupByCountry" in config:
        cluster_by_country = bool(config["GroupByCountry"])
    else:
        warnings.warn('"GroupByCountry" not found in config file. Continuing without clustering by country.')
        cluster_by_country = False

    # ROI/buffer masks based on original labels
    finite = np.isfinite(orig_sv_np)
    buffer_threshold = int(np.nanmax(orig_sv_np)) - int(config["nBufferClusters"])

    is_buffer = finite & (orig_sv_np > buffer_threshold)
    is_roi = finite & (orig_sv_np <= last_ROI_element)

    buffer_labels_np = np.where(is_buffer, orig_sv_np, 0.0)  # keep original buffer IDs
    labels_np = np.where(is_roi, 0.0, np.nan)                # new ROI labels live here
    sv_np = np.where(is_roi, orig_sv_np, 0.0)                # used when agg_level==1

    # Map element sensitivities to the SV grid once
    sensi_da = map_sensitivities_to_sv(flat_sensi, orig_sv, last_ROI_element)
    sensi_np = sensi_da.values
    print(f"Reducing to {desired_num_labels} elements")

    lon_all, lat_all, grid_shape = precompute_flat_lonlat(config, orig_sv_ds)

    sensi_flat = sensi_np.reshape(-1)
    labels_flat = labels_np.reshape(-1)
    roi_flat_idx = np.flatnonzero(is_roi.reshape(-1))

    # Geometry scaling: typical neighbor spacing in xyz
    gridstep_xyz = estimate_gridstep_xyz_from_roi(lon_all, lat_all, roi_flat_idx, sample_n=4000, random_state=0)

    # Iterate from fine (agg_level=1) to coarser (agg_level=max_cluster_size)
    cluster_pairs = np.arange(1, max_cluster_size + 1)
    fill_grid = False
    country_cache = {}
    current_max_label = 0

    for agg_level in cluster_pairs:
        if agg_level == max_cluster_size:
            fill_grid = True

        subset_idx = roi_flat_idx[labels_flat[roi_flat_idx] == 0]
        assert np.all(labels_flat[subset_idx] == 0)

        elements_left = int(subset_idx.size)
        if elements_left == 0:
            break

        clusters_left = int(desired_num_labels - current_max_label)
        backfill_num = int(elements_left / max_cluster_size) if max_cluster_size > 0 else 0

        if fill_grid or (clusters_left <= backfill_num):
            print("Filling grid with remaining clusters.")
            out_labels = cluster_data_kmeans(
                config,
                sensi_flat,
                orig_sv_ds,
                int(clusters_left),
                mini_batch,
                cluster_by_country,
                _country_mask_cache=country_cache,
                lon_all=lon_all,
                lat_all=lat_all,
                subset_idx=subset_idx,
                grid_shape=grid_shape,
                max_cluster_size=max_cluster_size,
                gridstep_xyz=gridstep_xyz,
            )
            local_threshold = -1.0  # accept all remaining clusters

        elif agg_level == 1:
            out_labels = sv_np.astype(np.int32)
            local_threshold = dofs_threshold

        else:
            n_clusters = int(np.round(elements_left / agg_level))
            if n_clusters == 0:
                continue
            out_labels = cluster_data_kmeans(
                config,
                sensi_flat,
                orig_sv_ds,
                int(n_clusters),
                mini_batch,
                cluster_by_country,
                _country_mask_cache=country_cache,
                lon_all=lon_all,
                lat_all=lat_all,
                subset_idx=subset_idx,
                grid_shape=grid_shape,
                max_cluster_size=max_cluster_size,
                gridstep_xyz=gridstep_xyz,
            )
            local_threshold = dofs_threshold

        # Rank clusters by total sensitivity; keep those meeting threshold
        n_max_labels, n_highest, num_elements, _ = get_highest_labels_threshold(out_labels, sensi_np, local_threshold)

        # Never assign more than remaining label budget
        if len(n_max_labels) > clusters_left:
            n_max_labels = n_max_labels[:clusters_left]
            num_elements = num_elements[:clusters_left]

        # Reserve space for backfill
        if (clusters_left - backfill_num) < n_highest and not fill_grid:
            fill_grid = True
            n_ind = int(clusters_left - backfill_num)
            n_max_labels = n_max_labels[:n_ind]
            num_elements = num_elements[:n_ind]

        # Avoid corner-case over-assignment when remaining cells are limited
        if not fill_grid and len(num_elements) > 0:
            total_elements = 0
            keep_n = len(num_elements)
            for i, n in enumerate(num_elements):
                total_elements += int(n)
                if elements_left - total_elements < clusters_left:
                    keep_n = i
                    break
            n_max_labels = n_max_labels[:keep_n]

        if len(n_max_labels) == 0:
            continue

        print(f"assigning {len(n_max_labels)} labels with agg level: {agg_level}")

        label_start = current_max_label + 1

        out_flat_subset = np.asarray(out_labels).reshape(-1)[subset_idx].astype(np.int64, copy=False)
        max_out = int(out_flat_subset.max()) if out_flat_subset.size else 0
        if max_out <= 0:
            continue

        # LUT: old cluster id -> new sequential ROI id
        lut = np.full(max_out + 1, -1, dtype=np.int32)
        new_ids = np.arange(label_start, label_start + len(n_max_labels), dtype=np.int32)
        lut[np.asarray(n_max_labels, dtype=np.int64)] = new_ids

        valid_out = (out_flat_subset >= 1) & (out_flat_subset <= max_out)
        if not np.any(valid_out):
            continue

        mapped = lut[out_flat_subset[valid_out]]
        keep = mapped != -1
        if not np.any(keep):
            continue

        labels_flat[subset_idx[valid_out][keep]] = mapped[keep]
        current_max_label = int(label_start + len(n_max_labels) - 1)

        if current_max_label >= desired_num_labels:
            labels_np[labels_np == 0] = desired_num_labels
            current_max_label = desired_num_labels
            break

    # Final fill: any leftover ROI cells become the last ROI label
    if current_max_label < desired_num_labels:
        labels_np[labels_np == 0] = desired_num_labels
        current_max_label = desired_num_labels

    # Compress buffer labels to follow reduced ROI range
    cluster_number_diff = int(last_ROI_element - current_max_label)
    buf = buffer_labels_np.copy()
    buf_mask = np.isfinite(buf) & (buf > 0)
    buf[buf_mask] = buf[buf_mask] - cluster_number_diff

    # Merge: buffer where present, otherwise ROI labels
    statevector_np = np.where(buf_mask, buf, labels_np)

    # Write output dataset
    refyear = 2000
    fillvalue = -9999

    if config["UseGCHP"]:
        da_statevector = xr.DataArray(
            np.where(np.isfinite(statevector_np), statevector_np, fillvalue)[None, ...],
            dims=["time", "nf", "Ydim", "Xdim"],
            coords=dict(
                time=(["time"], [0.0]),
                lats=(["nf", "Ydim", "Xdim"], orig_sv_ds["lats"].values),
                lons=(["nf", "Ydim", "Xdim"], orig_sv_ds["lons"].values),
            ),
            attrs=dict(units="1", missing_value=fillvalue, _FillValue=fillvalue),
        )
        ds_statevector = xr.Dataset({"StateVector": da_statevector})

        ds_statevector.lats.attrs["units"] = "degrees_north"
        ds_statevector.lats.attrs["long_name"] = "Latitude"
        ds_statevector.lons.attrs["units"] = "degrees_east"
        ds_statevector.lons.attrs["long_name"] = "Longitude"
        ds_statevector["time"].attrs = dict(
            units="days since {}-01-01 00:00:00".format(refyear),
            delta_t="0000-01-00 00:00:00",
            axis="T",
            standard_name="Time",
            long_name="Time",
            calendar="standard",
        )

        if "corner_lats" in orig_sv_ds.variables:
            ds_statevector["corner_lats"] = orig_sv_ds["corner_lats"]
        if "corner_lons" in orig_sv_ds.variables:
            ds_statevector["corner_lons"] = orig_sv_ds["corner_lons"]

        if config.get("STRETCH_GRID", False):
            for k in ["STRETCH_FACTOR", "TARGET_LAT", "TARGET_LON"]:
                if k in config:
                    ds_statevector.attrs[k] = np.float32(config[k])

    else:
        da_statevector = xr.DataArray(
            np.where(np.isfinite(statevector_np), statevector_np, fillvalue),
            dims=orig_sv.dims,
            coords=orig_sv.coords,
            attrs=dict(units="1", missing_value=fillvalue, _FillValue=fillvalue),
        ).expand_dims(time=[0.0])

        ds_statevector = da_statevector.to_dataset(name="StateVector")
        ds_statevector["time"].attrs = dict(
            units="days since {}-01-01 00:00:00".format(refyear),
            delta_t="0000-01-00 00:00:00",
            axis="T",
            standard_name="Time",
            long_name="Time",
            calendar="standard",
        )
        ds_statevector.lat.attrs["units"] = "degrees_north"
        ds_statevector.lat.attrs["long_name"] = "Latitude"
        ds_statevector.lon.attrs["units"] = "degrees_east"
        ds_statevector.lon.attrs["long_name"] = "Longitude"

    return ds_statevector


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    try:
        config_path = sys.argv[1]
        native_state_vector_path = sys.argv[2]
        state_vector_path = sys.argv[3]
        preview_dir = sys.argv[4]
        tropomi_cache = sys.argv[5]
        kf_index = int(sys.argv[6]) if len(sys.argv) > 6 else None

        config = yaml.load(open(config_path), Loader=yaml.FullLoader)
        original_clusters_ds = xr.open_dataset(native_state_vector_path).squeeze()

        sensitivity_args = [config, native_state_vector_path, preview_dir, tropomi_cache, False]
        if kf_index is not None:
            print(f"Dynamically generating clusters for period: {kf_index}.")
            sensitivity_args.append(kf_index)

        sensitivities = estimate_averaging_kernel(*sensitivity_args)
        sensitivities = force_native_res_pixels(config, original_clusters_ds, sensitivities)

        print(
            "Creating clusters based on information content and spatial proximity.\n"
            "Using ClusteringMethod: 'mini-batch-kmeans' may be faster but less accurate."
        )

        new_sv = update_sv_clusters(config, sensitivities, original_clusters_ds)

        new_sv.to_netcdf(
            state_vector_path,
            encoding={v: {"zlib": True, "complevel": 1} for v in new_sv.data_vars},
        )

        original_clusters_ds.close()

    except Exception as err:
        with open(".aggregation_error.txt", "w") as f:
            f.write("This file is used to tell the controlling script that state vector clustering failed")
        print(err)
        sys.exit(1)
