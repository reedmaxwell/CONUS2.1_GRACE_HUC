"""
WY2003 GRACE (multi-product) vs ParFlow/CLM TWS Comparison
===========================================================
Loads three GRACE mascon products (JPL, CSR, GSFC) with per-product
uncertainty (JPL formal leakage, GSFC leakage+noise, CSR inter-product
spread), plus ParFlow/CLM model output — all pre-aggregated to HUC2
basins — and generates comparison plots for Water Year 2003.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import pandas as pd
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE = Path(__file__).parent
DATA = BASE / "data"
MODEL_NC = BASE / "pfclm_huc2_tws.nc"
OUT_DIR = BASE / "figures"
OUT_DIR.mkdir(exist_ok=True)

HUC2_NAMES = {
    "01": "New England",       "02": "Mid-Atlantic",
    "03": "South Atlantic-Gulf", "04": "Great Lakes",
    "05": "Ohio",              "06": "Tennessee",
    "07": "Upper Mississippi",  "08": "Lower Mississippi",
    "09": "Souris-Red-Rainy",  "10": "Missouri",
    "11": "Arkansas-White-Red", "12": "Texas-Gulf",
    "13": "Rio Grande",        "14": "Upper Colorado",
    "15": "Lower Colorado",    "16": "Great Basin",
    "17": "Pacific Northwest", "18": "California",
}
CODES = sorted(HUC2_NAMES.keys())
WEST_CODES = ["10", "13", "14", "15", "16", "17", "18"]
EAST_CODES = ["01", "02", "03", "04", "05", "06", "07", "08"]


def load_grace_csv(path):
    """Load a GRACE HUC2 WY2003 CSV and snap to first-of-month."""
    df = pd.read_csv(path, index_col=0, parse_dates=True)
    df.columns = [c.strip().zfill(2) for c in df.columns]
    df.index = df.index.to_period("M").to_timestamp()
    df.index.name = "time"
    return df


# ── Load GRACE products ──────────────────────────────────────────────────────
GRACE_PRODUCTS = {
    "JPL": {
        "csv": DATA / "grace_jpl_huc2_wy2003.csv",
        "color": "#b2182b",   # dark red
        "marker": "s",
        "label": "JPL RL06.3",
    },
    "CSR": {
        "csv": DATA / "grace_csr_huc2_wy2003.csv",
        "color": "#ef8a62",   # orange-red
        "marker": "^",
        "label": "CSR RL06.3",
    },
    "GSFC": {
        "csv": DATA / "grace_gsfc_huc2_wy2003.csv",
        "color": "#f4a582",   # salmon
        "marker": "D",
        "label": "GSFC RL06v2",
    },
}

grace_dfs = {}
grace_uncs = {}
print("Loading GRACE WY2003 products...")
for name, info in GRACE_PRODUCTS.items():
    df = load_grace_csv(info["csv"])
    grace_dfs[name] = df
    print(f"  {name}: {df.shape[0]} months ({df.index[0].strftime('%Y-%m')} to {df.index[-1].strftime('%Y-%m')})")

# Per-product uncertainty
#   JPL: formal leakage uncertainty from mascon product
#   GSFC: formal leakage + noise (from native HDF5 mascons)
#   CSR: no formal uncertainty — use inter-product spread as proxy
grace_uncs["JPL"] = load_grace_csv(DATA / "grace_jpl_huc2_wy2003_uncertainty.csv")
grace_uncs["GSFC"] = load_grace_csv(DATA / "grace_gsfc_huc2_wy2003_uncertainty.csv")

# CSR proxy: half-range of (max - min) across JPL/CSR/GSFC at each point
stacked = np.stack([grace_dfs[n].values for n in ["JPL", "CSR", "GSFC"]], axis=0)
spread = (stacked.max(axis=0) - stacked.min(axis=0)) / 2
grace_uncs["CSR"] = pd.DataFrame(spread, index=grace_dfs["CSR"].index,
                                  columns=grace_dfs["CSR"].columns)

for name in GRACE_PRODUCTS:
    src = {"JPL": "formal leakage", "GSFC": "leakage+noise",
           "CSR": "inter-product spread"}[name]
    print(f"  {name} uncertainty ({src}): mean {grace_uncs[name].values.mean():.2f} cm")

# ── Load ParFlow/CLM ──────────────────────────────────────────────────────────
import xarray as xr
print("\nLoading ParFlow/CLM model output...")
ds = xr.open_dataset(MODEL_NC)
model = ds["tws_cm"].to_dataframe().reset_index()
model["huc2_id"] = model["huc2_id"].apply(lambda x: f"{int(x):02d}")
model = model.pivot(index="time", columns="huc2_id", values="tws_cm")
model.index.name = "time"
print(f"  Model: {model.shape[0]} months ({model.index[0].strftime('%Y-%m')} to {model.index[-1].strftime('%Y-%m')})")
ds.close()

# ── Compute stats (model vs JPL on common months) ────────────────────────────
jpl = grace_dfs["JPL"]
common_months = jpl.index.intersection(model.index)
print(f"\n  Common months (model vs JPL): {len(common_months)}")

# ── Style setup ───────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "axes.spines.top": False,
    "axes.spines.right": False,
})

CLR_MODEL = "#2166ac"   # dark blue
CLR_UNC = "#d6604d"     # lighter red for uncertainty band

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 1: 18-panel multi-product comparison
# ══════════════════════════════════════════════════════════════════════════════
print("\nPlotting multi-panel comparison...")
ncols = 3
nrows = int(np.ceil(len(CODES) / ncols))

fig, axes = plt.subplots(
    nrows, ncols, figsize=(17, 3.2 * nrows),
    sharex=True, sharey=False,
    gridspec_kw={"hspace": 0.48, "wspace": 0.28},
)
axes_flat = axes.flatten()

UNC_ALPHAS = {"JPL": 0.15, "CSR": 0.10, "GSFC": 0.10}

for i, code in enumerate(CODES):
    ax = axes_flat[i]
    m = model[code]

    # Per-product uncertainty bands
    for name, info in GRACE_PRODUCTS.items():
        g = grace_dfs[name][code]
        u = grace_uncs[name][code]
        cidx = g.index.intersection(u.index)
        ax.fill_between(
            cidx, g.loc[cidx] - u.loc[cidx], g.loc[cidx] + u.loc[cidx],
            color=info["color"], alpha=UNC_ALPHAS[name], linewidth=0,
            zorder=1,
        )

    # Model line
    ax.plot(m.index, m.values, lw=2.0, color=CLR_MODEL,
            marker="o", ms=3.5, markeredgewidth=0, label="ParFlow/CLM", zorder=4)

    # GRACE products
    for name, info in GRACE_PRODUCTS.items():
        g = grace_dfs[name][code]
        ax.plot(g.index, g.values, lw=1.4, color=info["color"],
                marker=info["marker"], ms=3, markeredgewidth=0,
                label=info["label"], zorder=3)

    ax.axhline(0, color="0.5", lw=0.5, ls="--", zorder=0)

    # Stats: model vs JPL
    g_jpl = jpl[code]
    cidx = g_jpl.dropna().index.intersection(m.dropna().index)
    if len(cidx) > 2:
        r = np.corrcoef(m.loc[cidx], g_jpl.loc[cidx])[0, 1]
        rmse = np.sqrt(np.mean((m.loc[cidx] - g_jpl.loc[cidx]) ** 2))
        ax.text(
            0.98, 0.05, f"r = {r:.2f}\nRMSE = {rmse:.1f} cm",
            transform=ax.transAxes, fontsize=7,
            ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8", alpha=0.85),
        )

    ax.set_title(f"HUC {code} \u2013 {HUC2_NAMES[code]}",
                 fontsize=10, fontweight="bold", loc="left", pad=4)
    ax.set_ylabel("TWS Anomaly [cm]", fontsize=8)
    ax.tick_params(labelsize=8, length=3)
    ax.yaxis.set_major_locator(mticker.MaxNLocator(5, symmetric=True))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.grid(axis="y", color="0.92", lw=0.5, zorder=0)

# Legend in first panel
axes_flat[0].legend(
    fontsize=7, loc="upper left", frameon=True,
    fancybox=False, edgecolor="0.8", framealpha=0.9,
)

for j in range(i + 1, len(axes_flat)):
    axes_flat[j].set_visible(False)

for ax in axes_flat[len(CODES) - ncols : len(CODES)]:
    ax.set_xlabel("WY2003", fontsize=9)

fig.suptitle(
    "GRACE (JPL / CSR / GSFC) vs ParFlow/CLM Total Water Storage Anomalies \u2014 WY2003\n"
    "18 CONUS HUC2 Regions  |  Shading = per-product uncertainty  |  Stats vs JPL RL06.3",
    fontsize=12, fontweight="bold", y=1.01,
)
out = OUT_DIR / "wy2003_grace_vs_pfclm_panels.png"
fig.savefig(out, dpi=250, bbox_inches="tight", facecolor="white")
print(f"  Saved: {out.name}")
plt.close(fig)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 2: CONUS-mean time series (all products)
# ══════════════════════════════════════════════════════════════════════════════
print("Plotting CONUS-mean time series...")

model_mean = model.mean(axis=1)
grace_means = {name: df.mean(axis=1) for name, df in grace_dfs.items()}
unc_means = {name: df.mean(axis=1) for name, df in grace_uncs.items()}

fig5, ax5 = plt.subplots(figsize=(11, 5))

# Per-product uncertainty bands
for name, info in GRACE_PRODUCTS.items():
    gm = grace_means[name]
    um = unc_means[name]
    cidx = gm.index.intersection(um.index)
    ax5.fill_between(
        cidx, gm.loc[cidx] - um.loc[cidx], gm.loc[cidx] + um.loc[cidx],
        color=info["color"], alpha=UNC_ALPHAS[name], linewidth=0,
    )

# Model
ax5.plot(model_mean.index, model_mean.values, lw=2.5, color=CLR_MODEL,
         marker="o", ms=5, markeredgewidth=0, label="ParFlow/CLM", zorder=4)

# GRACE products
for name, info in GRACE_PRODUCTS.items():
    gm = grace_means[name]
    ax5.plot(gm.index, gm.values, lw=1.8, color=info["color"],
             marker=info["marker"], ms=4, markeredgewidth=0,
             label=f"GRACE {info['label']}", zorder=3)

ax5.axhline(0, color="0.5", lw=0.5, ls="--")
ax5.set_ylabel("TWS Anomaly [cm]", fontsize=11)
ax5.set_title(
    "CONUS-Mean TWS Anomaly \u2014 WY2003 (Oct 2002 \u2013 Sep 2003)",
    fontsize=12, fontweight="bold",
)
ax5.legend(fontsize=9, frameon=True, fancybox=False, edgecolor="0.8")
ax5.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
ax5.tick_params(labelsize=9)
ax5.grid(axis="y", color="0.92", lw=0.5)
fig5.tight_layout()
out5 = OUT_DIR / "wy2003_conus_mean_timeseries.png"
fig5.savefig(out5, dpi=250, bbox_inches="tight", facecolor="white")
print(f"  Saved: {out5.name}")
plt.close(fig5)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 3: Western basins overlay (model + 3 GRACE products)
# ══════════════════════════════════════════════════════════════════════════════
print("Plotting western basins overlay...")
cmap_w = plt.cm.Dark2(np.linspace(0, 0.85, len(WEST_CODES)))

fig2, (ax2a, ax2b, ax2c, ax2d) = plt.subplots(4, 1, figsize=(13, 14), sharex=True)
ax_panels = [("ParFlow/CLM", ax2a, model),
             ("GRACE JPL RL06.3", ax2b, grace_dfs["JPL"]),
             ("GRACE CSR RL06.3", ax2c, grace_dfs["CSR"]),
             ("GRACE GSFC RL06v2", ax2d, grace_dfs["GSFC"])]

for title, ax, df in ax_panels:
    for j, code in enumerate(WEST_CODES):
        ax.plot(df.index, df[code], lw=1.5, color=cmap_w[j],
                marker="o", ms=2.5, markeredgewidth=0,
                label=f"HUC {code} \u2013 {HUC2_NAMES[code]}")
    ax.axhline(0, color="0.5", lw=0.5, ls="--")
    ax.set_ylabel("TWS [cm]", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold", loc="left")
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7, ncol=4, loc="upper left", framealpha=0.9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.grid(axis="y", color="0.92", lw=0.5)

fig2.suptitle("WY2003 TWS Anomalies \u2014 Western Basins",
              fontsize=13, fontweight="bold", y=1.0)
fig2.tight_layout()
out2 = OUT_DIR / "wy2003_western_comparison.png"
fig2.savefig(out2, dpi=200, bbox_inches="tight", facecolor="white")
print(f"  Saved: {out2.name}")
plt.close(fig2)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 4: Eastern basins overlay
# ══════════════════════════════════════════════════════════════════════════════
print("Plotting eastern basins overlay...")
cmap_e = plt.cm.tab10(np.linspace(0, 0.85, len(EAST_CODES)))

fig3, (ax3a, ax3b, ax3c, ax3d) = plt.subplots(4, 1, figsize=(13, 14), sharex=True)
ax_panels_e = [("ParFlow/CLM", ax3a, model),
               ("GRACE JPL RL06.3", ax3b, grace_dfs["JPL"]),
               ("GRACE CSR RL06.3", ax3c, grace_dfs["CSR"]),
               ("GRACE GSFC RL06v2", ax3d, grace_dfs["GSFC"])]

for title, ax, df in ax_panels_e:
    for j, code in enumerate(EAST_CODES):
        ax.plot(df.index, df[code], lw=1.5, color=cmap_e[j],
                marker="o", ms=2.5, markeredgewidth=0,
                label=f"HUC {code} \u2013 {HUC2_NAMES[code]}")
    ax.axhline(0, color="0.5", lw=0.5, ls="--")
    ax.set_ylabel("TWS [cm]", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold", loc="left")
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7, ncol=4, loc="upper left", framealpha=0.9)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.grid(axis="y", color="0.92", lw=0.5)

fig3.suptitle("WY2003 TWS Anomalies \u2014 Eastern Basins",
              fontsize=13, fontweight="bold", y=1.0)
fig3.tight_layout()
out3 = OUT_DIR / "wy2003_eastern_comparison.png"
fig3.savefig(out3, dpi=200, bbox_inches="tight", facecolor="white")
print(f"  Saved: {out3.name}")
plt.close(fig3)

# ══════════════════════════════════════════════════════════════════════════════
# PLOT 5: Difference panels (model minus each GRACE product)
# ══════════════════════════════════════════════════════════════════════════════
print("Plotting difference panels...")
fig4, axes4 = plt.subplots(nrows, ncols, figsize=(17, 3.0 * nrows), sharex=True)
axes4_flat = axes4.flatten()

diff_colors = {"JPL": "#b2182b", "CSR": "#ef8a62", "GSFC": "#f4a582"}

for i, code in enumerate(CODES):
    ax = axes4_flat[i]
    for name, info in GRACE_PRODUCTS.items():
        g = grace_dfs[name][code]
        cidx = g.dropna().index.intersection(model[code].dropna().index)
        d = model[code].loc[cidx] - g.loc[cidx]
        ax.plot(d.index, d.values, lw=1.3, color=info["color"],
                marker=info["marker"], ms=2.5, markeredgewidth=0,
                label=f"vs {info['label']}")
    ax.axhline(0, color="0.5", lw=0.5, ls="--")
    ax.set_title(f"HUC {code} \u2013 {HUC2_NAMES[code]}",
                 fontsize=9, fontweight="bold", loc="left")
    ax.set_ylabel("Diff [cm]", fontsize=7)
    ax.tick_params(labelsize=7)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax.grid(axis="y", color="0.92", lw=0.5, zorder=0)

axes4_flat[0].legend(fontsize=7, loc="upper left", framealpha=0.9)
for j in range(i + 1, len(axes4_flat)):
    axes4_flat[j].set_visible(False)

fig4.suptitle(
    "ParFlow/CLM minus GRACE TWS \u2014 WY2003\n"
    "Difference vs each mascon product (common months)",
    fontsize=12, fontweight="bold", y=1.01,
)
fig4.tight_layout()
out4 = OUT_DIR / "wy2003_difference_panels.png"
fig4.savefig(out4, dpi=200, bbox_inches="tight", facecolor="white")
print(f"  Saved: {out4.name}")
plt.close(fig4)

# ── Summary stats (model vs each product) ────────────────────────────────────
print("\n── Per-HUC2 stats: model vs GRACE (common months) ──")
print(f"{'HUC':<6} {'Region':<22} {'JPL r':>7} {'CSR r':>7} {'GSFC r':>7}  {'JPL RMSE':>9} {'CSR RMSE':>9} {'GSFC RMSE':>10}")
for code in CODES:
    vals = {}
    for name in GRACE_PRODUCTS:
        g = grace_dfs[name][code]
        cidx = g.dropna().index.intersection(model[code].dropna().index)
        if len(cidx) > 2:
            r = np.corrcoef(model[code].loc[cidx], g.loc[cidx])[0, 1]
            rmse = np.sqrt(np.mean((model[code].loc[cidx] - g.loc[cidx]) ** 2))
            vals[name] = (r, rmse)
        else:
            vals[name] = (np.nan, np.nan)
    print(f"  {code}   {HUC2_NAMES[code]:<22} "
          f"{vals['JPL'][0]:+.3f} {vals['CSR'][0]:+.3f} {vals['GSFC'][0]:+.3f}  "
          f"{vals['JPL'][1]:.2f} cm  {vals['CSR'][1]:.2f} cm  {vals['GSFC'][1]:.2f} cm")

print("\nDone! All figures saved to:", OUT_DIR)
