# CONUS GRACE vs ParFlow/CLM TWS Comparison

Compare GRACE satellite-observed Total Water Storage anomalies with ParFlow/CLM model output across 18 CONUS HUC2 basins for Water Year 2003.

## Three GRACE Mascon Products

- **JPL RL06.3 V04 CRI** — with formal leakage uncertainty
- **CSR RL06.3** — from UT Center for Space Research
- **GSFC RL06v2.0** — from NASA Goddard, with leakage + noise uncertainty

All products aggregated to HUC2 basins using area-weighted (cos-latitude) spatial averaging.

## Quick Start

```bash
conda activate subsettools

# Step 1: Aggregate GRACE products to HUC2 (requires raw NCs in data/)
python aggregate_grace_products.py

# Step 2: Generate all comparison plots
python plot_wy2003_comparison.py
```

## Directory Layout

```
├── grace_pfclm_wy2003.new.ipynb  # Notebook: PF-CLM TWS build + GRACE comparison (Verde/HPC)
├── plot_wy2003_comparison.py      # Main comparison plots (5 figures)
├── aggregate_grace_products.py    # GRACE → HUC2 aggregation pipeline
├── pfclm_huc2_tws.nc             # ParFlow/CLM model output
├── grace_huc2_handoff.md          # Detailed documentation
├── data/                          # GRACE CSVs and raw NetCDFs
├── figures/                       # Output plots and CSVs
└── _old/                          # Archived previous scripts
```

## Figures

| Figure | Description |
|--------|-------------|
| `wy2003_grace_vs_pfclm_panels.png` | 18-panel comparison with per-product uncertainty bands |
| `wy2003_conus_mean_timeseries.png` | CONUS-mean time series (all products) |
| `wy2003_western_comparison.png` | Western basins overlay (4 panels: model + 3 GRACE) |
| `wy2003_eastern_comparison.png` | Eastern basins overlay |
| `wy2003_difference_panels.png` | Model minus GRACE difference by HUC2 |

## PF-CLM TWS Computation Notes

The notebook (`grace_pfclm_wy2003.new.ipynb`) computes daily TWS on the CONUS2.1 1km grid
from ParFlow/CLM daily-average PFBs, then aggregates to HUC2 monthly means.

**TWS formula** (per grid cell, per day):
```
TWS_mm = (SUBstorage + SURF_WATstorage) / CELL_AREA * 1000 + SWE
```
- `SUBstorage`: column-summed subsurface storage (assumed m³ per cell) — **verify units against the script that produced the daily PFBs** (`compute_daily_PF/CLM_averages.py`)
- `SURF_WATstorage`: ponded surface water (assumed m³ per cell)
- `SWE`: CLM snow water equivalent (mm)
- `CELL_AREA`: 1e6 m² (1km × 1km LCC grid)

**HUC2 aggregation**: arithmetic mean over all mask cells per basin.
Equal-area assumption is valid on the LCC projection (no cos-lat weighting needed).
Mask sourced from `hf_hydrodata` (`huc_mapping`, `conus2` grid, level 2).

**Known caveats**:
- SWE nodata (−9999) is set to 0 rather than NaN before averaging. Could slightly bias coastal/edge basins where HUC2 mask cells fall outside the ParFlow active domain.
- No explicit intersection with a ParFlow domain mask — relies on the hf_hydrodata HUC2 mask aligning with the active domain.
- Daily PFBs run to day 364 (Sep 29); day 365 (Sep 30) is missing, so the September monthly mean uses 29 days.
- Anomalies are computed by subtracting the WY2003 annual mean (per HUC2). This is a different reference than GRACE products (which use a ~2004–2009 baseline), so only temporal patterns (correlations) are directly comparable, not absolute offsets.

## Dependencies

```
xarray netcdf4 matplotlib numpy pandas rasterio shapely h5py
```
