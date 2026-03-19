"""
GRACE TWS Anomalies by HUC2 — Full Pipeline
"""
import numpy as np
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
from rasterio.features import geometry_mask
from rasterio.transform import from_bounds
from shapely.geometry import mapping

# ── Config ───────────────────────────────────────────────────────────────────
GRACE_NC = "/mnt/user-data/uploads/conus_grace.nc"
HUC2_GEOJSON = "/home/claude/watersheds-biodiversity/data/us-huc2.geojson"
TWS_VAR = "lwe_thickness"

CONUS_HUC2 = [f"{i:02d}" for i in range(1, 19)]

HUC2_NAMES = {
    "01": "New England", "02": "Mid-Atlantic", "03": "South Atlantic-Gulf",
    "04": "Great Lakes", "05": "Ohio", "06": "Tennessee",
    "07": "Upper Mississippi", "08": "Lower Mississippi", "09": "Souris-Red-Rainy",
    "10": "Missouri", "11": "Arkansas-White-Red", "12": "Texas-Gulf",
    "13": "Rio Grande", "14": "Upper Colorado", "15": "Lower Colorado",
    "16": "Great Basin", "17": "Pacific Northwest", "18": "California",
}

# ── Load data ────────────────────────────────────────────────────────────────
print("Loading GRACE...")
ds = xr.open_dataset(GRACE_NC)

# Convert lon 0-360 → -180 to 180
if ds.lon.values.max() > 180:
    ds = ds.assign_coords(lon=(((ds.lon + 180) % 360) - 180))
    ds = ds.sortby("lon")

print(f"  Grid: {ds.sizes['lat']} x {ds.sizes['lon']}, {ds.sizes['time']} timesteps")
print(f"  Lon: {ds.lon.values.min():.2f} to {ds.lon.values.max():.2f}")
print(f"  Lat: {ds.lat.values.min():.2f} to {ds.lat.values.max():.2f}")

print("Loading HUC2 boundaries...")
gdf = gpd.read_file(HUC2_GEOJSON)
gdf = gdf[gdf["HUC2"].isin(CONUS_HUC2)].copy()
gdf = gdf.to_crs("EPSG:4326")
print(f"  {len(gdf)} CONUS HUC2 regions")

# ── Build masks ──────────────────────────────────────────────────────────────
print("Rasterizing HUC2 masks onto GRACE grid...")
lon = ds.lon.values
lat = ds.lat.values
nlon, nlat = len(lon), len(lat)
dlon = abs(lon[1] - lon[0])
dlat = abs(lat[1] - lat[0])

transform = from_bounds(
    lon.min() - dlon/2, lat.min() - dlat/2,
    lon.max() + dlon/2, lat.max() + dlat/2,
    nlon, nlat,
)

masks = {}
for _, row in gdf.iterrows():
    code = row["HUC2"]
    geom = [mapping(row.geometry)]
    mask = ~geometry_mask(geom, out_shape=(nlat, nlon), transform=transform)
    # rasterio rasterizes top-down; flip if lat is ascending
    if lat[0] < lat[-1]:
        mask = mask[::-1, :]
    masks[code] = mask
    ncells = mask.sum()
    print(f"  HUC2 {code} ({HUC2_NAMES[code]}): {ncells} grid cells")

# ── Aggregate ────────────────────────────────────────────────────────────────
print("Computing area-weighted spatial averages...")
da = ds[TWS_VAR]
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

df = pd.DataFrame(records, index=pd.to_datetime(ds.time.values))
df = df[sorted(df.columns)]

csv_path = "/home/claude/grace_huc2_tws.csv"
df.to_csv(csv_path)
print(f"Saved CSV: {csv_path}")
print(f"  Shape: {df.shape}")

# ── Plot: Multi-panel ────────────────────────────────────────────────────────
print("Plotting multi-panel figure...")
codes = sorted(df.columns)
ncols = 3
nrows = int(np.ceil(len(codes) / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(18, 2.8 * nrows), sharex=True)
axes = axes.flatten()

for i, code in enumerate(codes):
    ax = axes[i]
    ts = df[code]
    ax.plot(ts.index, ts.values, lw=0.7, color="steelblue")
    ax.axhline(0, color="k", lw=0.4, ls="--")
    ax.fill_between(ts.index, ts.values, 0,
                    where=ts.values >= 0, color="steelblue", alpha=0.2)
    ax.fill_between(ts.index, ts.values, 0,
                    where=ts.values < 0, color="tomato", alpha=0.2)
    label = f"HUC {code} – {HUC2_NAMES.get(code, '')}"
    ax.set_title(label, fontsize=9, fontweight="bold", loc="left")
    ax.set_ylabel("TWS [cm]", fontsize=7)
    ax.tick_params(labelsize=7)
    ax.xaxis.set_major_locator(mdates.YearLocator(5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))

for j in range(i + 1, len(axes)):
    axes[j].set_visible(False)

fig.suptitle(
    "GRACE/GRACE-FO TWS Anomalies by HUC2 Region (CONUS)\n"
    "JPL Mascon RL06.3 V04 CRI  •  2002–2025",
    fontsize=13, fontweight="bold", y=1.02,
)
fig.tight_layout()
out1 = "/home/claude/grace_huc2_all.png"
fig.savefig(out1, dpi=200, bbox_inches="tight", facecolor="white")
print(f"Saved: {out1}")

# ── Plot: Western basins overlay ─────────────────────────────────────────────
print("Plotting western basins overlay...")
west_codes = ["10", "14", "15", "16", "17", "18"]
fig2, ax2 = plt.subplots(figsize=(14, 5))
colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"]
for j, code in enumerate(west_codes):
    if code not in df.columns:
        continue
    label = f"HUC {code} – {HUC2_NAMES.get(code, '')}"
    ax2.plot(df.index, df[code], lw=1.2, label=label, color=colors[j])

ax2.axhline(0, color="k", lw=0.5, ls="--")
ax2.set_ylabel("TWS Anomaly [cm]", fontsize=11)
ax2.set_title("GRACE/GRACE-FO TWS Anomalies — Western HUC2 Basins", fontsize=12, fontweight="bold")
ax2.legend(fontsize=9, ncol=2, loc="lower left")
ax2.xaxis.set_major_locator(mdates.YearLocator(2))
ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax2.tick_params(labelsize=9)
fig2.tight_layout()
out2 = "/home/claude/grace_huc2_western.png"
fig2.savefig(out2, dpi=200, bbox_inches="tight", facecolor="white")
print(f"Saved: {out2}")

# ── Plot: Eastern basins overlay ─────────────────────────────────────────────
print("Plotting eastern basins overlay...")
east_codes = ["01", "02", "03", "04", "05", "06", "07", "08"]
fig3, ax3 = plt.subplots(figsize=(14, 5))
colors_e = ["#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928", "#a6cee3", "#b2df8a"]
for j, code in enumerate(east_codes):
    if code not in df.columns:
        continue
    label = f"HUC {code} – {HUC2_NAMES.get(code, '')}"
    ax3.plot(df.index, df[code], lw=1.2, label=label, color=colors_e[j])

ax3.axhline(0, color="k", lw=0.5, ls="--")
ax3.set_ylabel("TWS Anomaly [cm]", fontsize=11)
ax3.set_title("GRACE/GRACE-FO TWS Anomalies — Eastern HUC2 Basins", fontsize=12, fontweight="bold")
ax3.legend(fontsize=9, ncol=2, loc="lower left")
ax3.xaxis.set_major_locator(mdates.YearLocator(2))
ax3.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
ax3.tick_params(labelsize=9)
fig3.tight_layout()
out3 = "/home/claude/grace_huc2_eastern.png"
fig3.savefig(out3, dpi=200, bbox_inches="tight", facecolor="white")
print(f"Saved: {out3}")

print("\nDone!")
