"""
Aggregate multiple GRACE mascon products to HUC2 basins
=======================================================
Processes JPL, CSR, and GSFC mascon solutions using identical
area-weighted (cos-lat) spatial averaging over HUC2 boundaries.
Outputs a single CSV per product for WY2003.
"""
import json
import numpy as np
import xarray as xr
import pandas as pd
from pathlib import Path
from rasterio.features import geometry_mask
from rasterio.transform import from_bounds

BASE = Path(__file__).parent
OUT_DIR = BASE / "data"
HUC2_GEOJSON = Path.home() / "Projects" / "WTD_viewer" / "web_app" / "huc2.geojson"

CONUS_HUC2 = [f"{i:02d}" for i in range(1, 19)]

# ── Product definitions ───────────────────────────────────────────────────────
PRODUCTS = {
    "JPL": {
        "path": BASE / "_old" / "GRCTellus.JPL.200204_202512.GLO.RL06.3M.MSCNv04CRI.nc",
        "tws_var": "lwe_thickness",
        "unc_var": "uncertainty",
        "label": "JPL RL06.3 CRI",
    },
    "CSR": {
        "path": BASE / "data" / "CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc",
        "tws_var": "lwe_thickness",
        "unc_var": None,
        "label": "CSR RL06.3",
    },
    "GSFC": {
        "path": BASE / "data" / "GSFC_mascon_halfdegree.nc",
        "tws_var": "lwe_thickness",
        "unc_var": None,
        "label": "GSFC RL06v2.0",
    },
}


def load_and_normalize(product_info):
    """Load a GRACE product and normalize lon to -180..180, time to datetime."""
    ds = xr.open_dataset(product_info["path"])

    # Fix CSR time (numeric days since 2002-01-01, stored with capital 'Units')
    if not np.issubdtype(ds.time.dtype, np.datetime64):
        units_attr = ds.time.attrs.get("Units", ds.time.attrs.get("units", ""))
        if "days since" in units_attr:
            ref = pd.Timestamp(units_attr.split("since")[-1].strip().replace("T", " ").rstrip("Z"))
            ds["time"] = pd.to_datetime(ref + pd.to_timedelta(ds.time.values, unit="D"))

    # Convert lon 0-360 -> -180..180
    if ds.lon.values.max() > 180:
        ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180))
        ds = ds.sortby("lon")

    return ds


def build_masks(lon, lat, huc2_features):
    """Build rasterized HUC2 masks on a given lon/lat grid."""
    nlon, nlat = len(lon), len(lat)
    dlon = abs(float(lon[1] - lon[0]))
    dlat = abs(float(lat[1] - lat[0]))

    transform = from_bounds(
        float(lon.min()) - dlon / 2, float(lat.min()) - dlat / 2,
        float(lon.max()) + dlon / 2, float(lat.max()) + dlat / 2,
        nlon, nlat,
    )

    masks = {}
    for code, geom in huc2_features.items():
        mask = ~geometry_mask([geom], out_shape=(nlat, nlon), transform=transform)
        if lat[0] < lat[-1]:
            mask = mask[::-1, :]
        masks[code] = mask
    return masks


def aggregate_to_huc2(da, masks, lat):
    """Area-weighted (cos-lat) spatial average over HUC2 masks."""
    nlat, nlon = len(lat), da.shape[-1]
    weights_1d = np.cos(np.deg2rad(lat))
    weights_2d = np.broadcast_to(weights_1d[:, None], (nlat, nlon))

    records = {}
    for code, mask in masks.items():
        w = weights_2d * mask.astype(float)
        w_sum = w.sum()
        if w_sum == 0:
            continue
        vals = (da.values * w[None, :, :]).sum(axis=(1, 2)) / w_sum
        records[code] = vals
    return records


# ── Load HUC2 boundaries ─────────────────────────────────────────────────────
print("Loading HUC2 boundaries...")
with open(HUC2_GEOJSON) as f:
    huc2_gj = json.load(f)

huc2_features = {
    feat["properties"]["huc2"]: feat["geometry"]
    for feat in huc2_gj["features"]
    if feat["properties"]["huc2"] in CONUS_HUC2
}
print(f"  {len(huc2_features)} CONUS HUC2 regions")

# ── Process each product ─────────────────────────────────────────────────────
for name, info in PRODUCTS.items():
    print(f"\n{'='*60}")
    print(f"Processing {name}: {info['label']}")
    print(f"  File: {info['path'].name}")

    ds = load_and_normalize(info)
    lon = ds.lon.values
    lat = ds.lat.values
    print(f"  Grid: {len(lat)} x {len(lon)}, {len(ds.time)} timesteps")
    print(f"  Time: {pd.Timestamp(ds.time.values[0]).strftime('%Y-%m')} to "
          f"{pd.Timestamp(ds.time.values[-1]).strftime('%Y-%m')}")

    # Build masks for this grid resolution
    masks = build_masks(lon, lat, huc2_features)

    # Aggregate TWS
    da = ds[info["tws_var"]]
    records = aggregate_to_huc2(da, masks, lat)

    df_all = pd.DataFrame(records, index=pd.to_datetime(ds.time.values))
    df_all = df_all[sorted(df_all.columns)]
    df_all.index.name = "time"

    # Save full record
    full_csv = OUT_DIR / f"grace_{name.lower()}_huc2_tws.csv"
    df_all.to_csv(full_csv)
    print(f"  Saved full: {full_csv.name} ({df_all.shape})")

    # Filter to WY2003
    wy = df_all[(df_all.index >= "2002-10-01") & (df_all.index < "2003-10-01")]
    wy_csv = OUT_DIR / f"grace_{name.lower()}_huc2_wy2003.csv"
    wy.to_csv(wy_csv)
    print(f"  Saved WY2003: {wy_csv.name} ({wy.shape})")

    # Aggregate uncertainty if available
    if info["unc_var"] and info["unc_var"] in ds.data_vars:
        da_unc = ds[info["unc_var"]]
        unc_records = aggregate_to_huc2(da_unc, masks, lat)
        df_unc = pd.DataFrame(unc_records, index=pd.to_datetime(ds.time.values))
        df_unc = df_unc[sorted(df_unc.columns)]
        df_unc.index.name = "time"
        unc_wy = df_unc[(df_unc.index >= "2002-10-01") & (df_unc.index < "2003-10-01")]
        unc_csv = OUT_DIR / f"grace_{name.lower()}_huc2_wy2003_uncertainty.csv"
        unc_wy.to_csv(unc_csv)
        print(f"  Saved uncertainty: {unc_csv.name}")

    ds.close()

print("\nDone! All CSVs saved to:", OUT_DIR)
