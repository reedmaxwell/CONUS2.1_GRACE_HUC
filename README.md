# CONUS2.1 GRACE vs ParFlow/CLM TWS Comparison

Compare GRACE satellite-observed Total Water Storage (TWS) anomalies with ParFlow/CLM model output across 18 CONUS HUC2 basins for Water Year 2003.

## Three GRACE Mascon Products

- **JPL RL06.3 V04 CRI** — with formal leakage uncertainty
- **CSR RL06.3** — from UT Center for Space Research
- **GSFC RL06v2.0** — from NASA Goddard, with leakage + noise uncertainty

All products aggregated to HUC2 basins using area-weighted (cos-latitude) spatial averaging.

## Quick Start

### GRACE aggregation & comparison plots (local)

```bash
conda activate subsettools

# Step 1: Aggregate GRACE products to HUC2 (requires raw NCs in data/)
python aggregate_grace_products.py

# Step 2: Generate all comparison plots
python plot_wy2003_comparison.py
```

Raw GRACE NetCDFs are not included in this repo (too large). Download them to `data/` — see [grace_huc2_handoff.md](grace_huc2_handoff.md) for URLs.

### PF-CLM TWS build (Verde/HPC)

The notebook `build_pfclm_tws.ipynb` builds the model TWS output from CONUS2.1 daily PFBs on Verde with direct `/hydrodata` access. It produces `pfclm_huc2_tws.nc` (already included in this repo). Runtime is ~25 min.

## Directory Layout

```
├── build_pfclm_tws.ipynb          # Notebook: PF-CLM TWS build (runs on Verde/HPC)
├── aggregate_grace_products.py    # GRACE -> HUC2 aggregation pipeline
├── plot_wy2003_comparison.py      # Main comparison plots (5 figures)
├── pfclm_huc2_tws.nc             # ParFlow/CLM model output (pre-built)
├── grace_huc2_handoff.md          # Detailed documentation
├── data/                          # Processed HUC2 CSVs (raw NCs not tracked)
└── figures/                       # Output plots and CSVs
```

## Figures

| Figure | Description |
|--------|-------------|
| `wy2003_grace_vs_pfclm_panels.png` | 18-panel comparison with per-product uncertainty bands |
| `wy2003_conus_mean_timeseries.png` | CONUS-mean time series (all products) |
| `wy2003_western_comparison.png` | Western basins overlay (4 panels: model + 3 GRACE) |
| `wy2003_eastern_comparison.png` | Eastern basins overlay |
| `wy2003_difference_panels.png` | Model minus GRACE difference by HUC2 |

## PF-CLM TWS Computation

See `build_pfclm_tws.ipynb` for the full workflow. Summary:

**TWS formula** (per grid cell, per day):
```
TWS_mm = (SUBstorage + SURF_WATstorage) / CELL_AREA * 1000 + SWE
```
- `SUBstorage`: column-summed subsurface storage (m³ per cell)
- `SURF_WATstorage`: ponded surface water (m³ per cell)
- `SWE`: CLM snow water equivalent (mm)
- `CELL_AREA`: 1e6 m² (1km x 1km LCC grid)

**HUC2 aggregation**: arithmetic mean over all mask cells per basin.
Equal-area assumption is valid on the LCC projection (no cos-lat weighting needed).
Mask sourced from `hf_hydrodata` (`huc_mapping`, `conus2` grid, level 2).

**Known caveats**:
- SWE nodata (-9999) is set to 0 rather than NaN before averaging. Could slightly bias coastal/edge basins where HUC2 mask cells fall outside the ParFlow active domain.
- Daily PFBs run to day 364 (Sep 29); day 365 (Sep 30) is missing, so the September monthly mean uses 29 days.
- Anomalies are computed by subtracting the WY2003 annual mean (per HUC2). This is a different reference than GRACE products (which use a ~2004-2009 baseline), so only temporal patterns (correlations) are directly comparable, not absolute offsets.

## Dependencies

```
xarray netcdf4 matplotlib numpy pandas rasterio shapely h5py
```

The PF-CLM notebook additionally requires `hf_hydrodata` and `parflow` (available on Verde).
