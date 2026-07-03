"""Sentinel-2 L2A NDVI pipeline for the Mandalay region (Myanmar).

Queries Microsoft Planetary Computer for the most recent cloud-free
Sentinel-2 Level-2A scene that fully covers the area of interest,
downloads the Red (B04), NIR (B08) and Scene Classification (SCL)
bands clipped to the AOI, computes NDVI, and exports:

  - data/S2_<date>_B04.tif / _B08.tif / _SCL.tif  (clipped input bands)
  - data/NDVI_<date>.tif                           (NDVI GeoTIFF, EPSG:4326)
  - examples/NDVI_Mandalay_<date>.jpg              (publication-ready map)

Scene selection: among scenes whose data footprint fully contains the
AOI, keep those with the lowest scene cloud cover (within 0.5 % of the
minimum) and no meaningful radiometric degradation, then take the most
recent. Cloud cover inside the AOI is re-checked from the SCL band.

Radiometry: processing baseline >= 04.00 stores BOA reflectance with a
+1000 DN offset; it is removed before the NDVI ratio is formed.
"""

import os

import numpy as np
import planetary_computer
import pystac_client
import rasterio
import rasterio.warp
import rasterio.windows
from rasterio.enums import Resampling
from shapely.geometry import box, shape

BBOX = [95.964203, 21.725697, 96.330872, 22.104726]  # W, S, N, E (lon/lat)
SEARCH_WINDOW = "2025-07-01/2026-07-03"
MAX_SCENE_CLOUD = 10.0     # % — STAC pre-filter
CLOUD_TIE_BAND = 0.5       # % — scenes within this of the minimum tie on cloud
MAX_DEGRADED = 0.5         # % — s2:degraded_msi_data_percentage ceiling
SCL_CLOUD_CLASSES = (3, 8, 9, 10)  # shadow, cloud med/high prob, cirrus
STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"

DATA_DIR = "data"
EXAMPLES_DIR = "examples"


def select_scene():
    """Return the best STAC item: full AOI coverage, minimal cloud, most recent."""
    aoi = box(*BBOX)
    catalog = pystac_client.Client.open(
        STAC_URL, modifier=planetary_computer.sign_inplace
    )
    search = catalog.search(
        collections=["sentinel-2-l2a"],
        bbox=BBOX,
        datetime=SEARCH_WINDOW,
        query={"eo:cloud_cover": {"lt": MAX_SCENE_CLOUD}},
    )
    items = [
        it
        for it in search.items()
        if shape(it.geometry).contains(aoi)
        and it.properties.get("s2:degraded_msi_data_percentage", 0) < MAX_DEGRADED
    ]
    if not items:
        raise RuntimeError(
            "No scene fully covers the AOI within the search window; "
            "widen SEARCH_WINDOW or relax MAX_SCENE_CLOUD."
        )
    min_cloud = min(it.properties["eo:cloud_cover"] for it in items)
    candidates = [
        it for it in items if it.properties["eo:cloud_cover"] <= min_cloud + CLOUD_TIE_BAND
    ]
    best = max(candidates, key=lambda it: it.properties["datetime"])
    p = best.properties
    print(
        f"Selected {best.id}\n"
        f"  date {p['datetime'][:10]}  tile {p.get('s2:mgrs_tile')}  "
        f"cloud {p['eo:cloud_cover']:.2f}%  "
        f"degraded {p.get('s2:degraded_msi_data_percentage', 0):.3f}%  "
        f"baseline {p.get('s2:processing_baseline')}"
    )
    return best


def download_band(item, asset_key, out_path, out_shape=None):
    """Windowed read of one COG asset clipped to BBOX; saves a GeoTIFF clip."""
    href = item.assets[asset_key].href
    with rasterio.open(href) as src:
        aoi_bounds = rasterio.warp.transform_bounds("EPSG:4326", src.crs, *BBOX)
        window = rasterio.windows.from_bounds(*aoi_bounds, transform=src.transform)
        window = window.round_offsets().round_lengths()
        shape_ = out_shape or (int(window.height), int(window.width))
        data = src.read(
            1, window=window, out_shape=shape_, resampling=Resampling.nearest
        )
        transform = src.window_transform(window)
        if out_shape:
            transform = transform * transform.scale(
                window.width / shape_[1], window.height / shape_[0]
            )
        profile = {
            "driver": "GTiff",
            "height": shape_[0],
            "width": shape_[1],
            "count": 1,
            "dtype": data.dtype,
            "crs": src.crs,
            "transform": transform,
            "compress": "deflate",
            "nodata": 0,
        }
    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(data, 1)
    print(f"  {asset_key} -> {out_path}  {data.shape[1]}x{data.shape[0]} px")
    return data, profile


def compute_ndvi(red_dn, nir_dn, scl, baseline):
    """NDVI from L2A DNs; removes the +1000 BOA offset for baseline >= 04.00."""
    offset = 1000.0 if baseline and float(baseline) >= 4.0 else 0.0
    valid = (red_dn > 0) & (nir_dn > 0) & ~np.isin(scl, SCL_CLOUD_CLASSES)
    red = np.clip(red_dn.astype("float64") - offset, 0, None)
    nir = np.clip(nir_dn.astype("float64") - offset, 0, None)
    denom = nir + red
    ndvi = np.full(red.shape, np.nan)
    ok = valid & (denom > 0)
    ndvi[ok] = (nir[ok] - red[ok]) / denom[ok]

    aoi_px = (red_dn > 0) & (nir_dn > 0)
    cloudy = float(np.isin(scl, SCL_CLOUD_CLASSES)[aoi_px].mean()) * 100
    print(f"  BOA offset removed: {offset:.0f} DN | cloud/shadow in AOI: {cloudy:.2f}%")
    return ndvi


def reproject_to_wgs84(ndvi, profile):
    """Reproject the NDVI array from UTM to EPSG:4326 for a lon/lat map."""
    src_transform, src_crs = profile["transform"], profile["crs"]
    dst_transform, width, height = rasterio.warp.calculate_default_transform(
        src_crs, "EPSG:4326", profile["width"], profile["height"], *rasterio.transform.array_bounds(profile["height"], profile["width"], src_transform)
    )
    dst = np.full((height, width), np.nan)
    rasterio.warp.reproject(
        ndvi,
        dst,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=dst_transform,
        dst_crs="EPSG:4326",
        src_nodata=np.nan,
        dst_nodata=np.nan,
        resampling=Resampling.bilinear,
    )
    # Crop to the exact requested bbox: the reprojected UTM footprint is
    # slightly rotated in WGS84 and leaves NaN slivers at the frame edges.
    crop = rasterio.windows.from_bounds(*BBOX, transform=dst_transform)
    crop = crop.round_offsets().round_lengths()
    r0, c0 = max(crop.row_off, 0), max(crop.col_off, 0)
    r1, c1 = min(r0 + crop.height, height), min(c0 + crop.width, width)
    dst = dst[r0:r1, c0:c1]
    dst_transform = rasterio.windows.transform(
        rasterio.windows.Window(c0, r0, c1 - c0, r1 - r0), dst_transform
    )
    dst_profile = {
        "driver": "GTiff",
        "height": dst.shape[0],
        "width": dst.shape[1],
        "count": 1,
        "dtype": "float32",
        "crs": "EPSG:4326",
        "transform": dst_transform,
        "compress": "deflate",
        "nodata": np.nan,
    }
    return dst, dst_profile


def export_map(ndvi, dst_profile, scene_date, out_path):
    """Render the publication JPEG: diverging brown-green map centered at 0."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    t = dst_profile["transform"]
    h, w = dst_profile["height"], dst_profile["width"]
    extent = [t.c, t.c + t.a * w, t.f + t.e * h, t.f]  # W, E, S, N

    # NDVI is a polarity measure: diverging colormap, neutral midpoint at 0.
    vmax = float(np.ceil(np.nanpercentile(ndvi, 99) * 10) / 10)
    vmin = float(np.floor(np.nanpercentile(ndvi, 1) * 10) / 10)
    vmin, vmax = min(vmin, -0.1), max(vmax, 0.5)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(9, 9.8), dpi=300)
    im = ax.imshow(ndvi, cmap="BrBG", norm=norm, extent=extent, interpolation="nearest")

    ax.set_title(
        "Vegetation Index (NDVI) — Mandalay Region, Myanmar",
        fontsize=13, fontweight="bold", pad=14,
    )
    ax.text(
        0.5, 1.012, f"Sentinel-2 L2A  ·  {scene_date}  ·  10 m resolution",
        transform=ax.transAxes, ha="center", fontsize=9.5, color="0.35",
    )
    ax.set_xlabel("Longitude (°E)", fontsize=9)
    ax.set_ylabel("Latitude (°N)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.grid(True, linestyle=":", linewidth=0.4, color="0.75", alpha=0.7)

    cbar = fig.colorbar(im, ax=ax, shrink=0.72, pad=0.02, aspect=32)
    cbar.set_label("NDVI", fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    for v, lbl in [(-0.05, "water"), (0.12, "bare / built"), (0.65, "dense vegetation")]:
        if vmin < v < vmax:
            cbar.ax.axhline(v, color="0.25", linewidth=0.6)
            cbar.ax.text(1.45, v, lbl, transform=cbar.ax.get_yaxis_transform(),
                         fontsize=7, va="center", color="0.3")

    # Scale bar (~10 km) and north arrow.
    km_per_deg = 111.32 * np.cos(np.deg2rad((extent[2] + extent[3]) / 2))
    bar_deg = 10.0 / km_per_deg
    x0 = extent[0] + 0.05 * (extent[1] - extent[0])
    y0 = extent[2] + 0.045 * (extent[3] - extent[2])
    ax.plot([x0, x0 + bar_deg], [y0, y0], color="black", linewidth=2.5,
            solid_capstyle="butt")
    ax.text(x0 + bar_deg / 2, y0 + 0.008 * (extent[3] - extent[2]), "10 km",
            ha="center", va="bottom", fontsize=8)
    ax.annotate(
        "N", xy=(0.955, 0.965), xytext=(0.955, 0.895), xycoords="axes fraction",
        ha="center", va="center", fontsize=11, fontweight="bold",
        arrowprops=dict(arrowstyle="-|>", color="black", linewidth=1.4),
    )

    ax.text(
        0.0, -0.075,
        "Data: Copernicus Sentinel-2 L2A via Microsoft Planetary Computer  ·  "
        "NDVI = (B08 − B04) / (B08 + B04)  ·  clouds & shadows masked (SCL)",
        transform=ax.transAxes, fontsize=7, color="0.4",
    )

    fig.savefig(out_path, bbox_inches="tight", pil_kwargs={"quality": 95})
    plt.close(fig)
    print(f"  map -> {out_path}")


def main():
    os.makedirs(DATA_DIR, exist_ok=True)
    os.makedirs(EXAMPLES_DIR, exist_ok=True)

    item = select_scene()
    date = item.properties["datetime"][:10]
    tag = date.replace("-", "")

    print("Downloading bands (windowed COG reads, clipped to AOI):")
    red, profile = download_band(item, "B04", f"{DATA_DIR}/S2_{tag}_B04.tif")
    nir, _ = download_band(item, "B08", f"{DATA_DIR}/S2_{tag}_B08.tif")
    # SCL is 20 m; resample to the 10 m grid so masks align pixel-for-pixel.
    scl, _ = download_band(item, "SCL", f"{DATA_DIR}/S2_{tag}_SCL.tif",
                           out_shape=red.shape)

    print("Computing NDVI:")
    ndvi = compute_ndvi(red, nir, scl, item.properties.get("s2:processing_baseline"))

    ndvi_wgs84, dst_profile = reproject_to_wgs84(ndvi, profile)
    ndvi_path = f"{DATA_DIR}/NDVI_{tag}.tif"
    with rasterio.open(ndvi_path, "w", **dst_profile) as dst:
        dst.write(ndvi_wgs84.astype("float32"), 1)
    print(f"  NDVI GeoTIFF -> {ndvi_path}")

    stats = ndvi_wgs84[np.isfinite(ndvi_wgs84)]
    print(
        f"  NDVI stats: min {stats.min():.3f}  mean {stats.mean():.3f}  "
        f"max {stats.max():.3f}"
    )

    export_map(ndvi_wgs84, dst_profile, date, f"{EXAMPLES_DIR}/NDVI_Mandalay_{tag}.jpg")


if __name__ == "__main__":
    main()
