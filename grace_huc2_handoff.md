# GRACE TWS vs ParFlow/CLM — HUC2 Comparison

## Objective

Compare GRACE TWS anomalies from three mascon products (JPL, CSR, GSFC) with ParFlow/CLM model output, both aggregated to 18 CONUS HUC2 basins, for Water Year 2003 (Oct 2002 – Sep 2003).

## GRACE Products

| Product | Version | Grid | Source | Uncertainty |
|---------|---------|------|--------|-------------|
| JPL Mascon | RL06.3 V04 CRI | 0.5° | PO.DAAC | Formal `uncertainty` field (leakage) |
| CSR Mascon | RL06.3 | 0.25° | UT CSR direct server | None distributed; inter-product spread used as proxy |
| GSFC Mascon | RL06v2.0 | 0.5° (gridded) / native | NASA Goddard | `leakage_2sigma` + `noise_2sigma` from HDF5 |

All products provide `lwe_thickness` in cm (liquid water equivalent), anomalies relative to 2004–2009 baseline.

### Download URLs (no Earthdata login required)

- **CSR**: `https://download.csr.utexas.edu/outgoing/grace/RL0603_mascons/CSR_GRACE_GRACE-FO_RL0603_Mascons_all-corrections.nc`
- **GSFC NC**: `https://earth.gsfc.nasa.gov/sites/default/files/geo/gsfc.glb_.200204_YYYYMM_rl06v2.0_obp-ice6gd_halfdegree.nc`
- **GSFC HDF5** (with uncertainty): `https://earth.gsfc.nasa.gov/sites/default/files/geo/gsfc.glb_.200204_YYYYMM_rl06v2.0_obp-ice6gd.h5`

## Model Output

- **File**: `pfclm_huc2_tws.nc`
- **Variables**: `tws_cm` (float32), `tws_mm` (float32), indexed by `huc2_id` (1–18) and `time`
- **Period**: WY2003 (12 monthly timesteps, Oct 2002 – Sep 2003)

## HUC2 Regions

18 CONUS HUC2 basins (codes 01–18). Boundaries from `~/Projects/WTD_viewer/web_app/huc2.geojson`.

```python
HUC2_NAMES = {
    "01": "New England", "02": "Mid-Atlantic", "03": "South Atlantic-Gulf",
    "04": "Great Lakes", "05": "Ohio", "06": "Tennessee",
    "07": "Upper Mississippi", "08": "Lower Mississippi", "09": "Souris-Red-Rainy",
    "10": "Missouri", "11": "Arkansas-White-Red", "12": "Texas-Gulf",
    "13": "Rio Grande", "14": "Upper Colorado", "15": "Lower Colorado",
    "16": "Great Basin", "17": "Pacific Northwest", "18": "California",
}
```

## Scripts

| Script | Purpose |
|--------|---------|
| `aggregate_grace_products.py` | Aggregate all 3 GRACE products from raw NC to HUC2 CSVs |
| `plot_wy2003_comparison.py` | Generate all comparison plots (panels, overlays, CONUS mean, differences) |

## Uncertainty Convention

| Product | Type | Mean (cm) | Method |
|---------|------|-----------|--------|
| JPL | Formal leakage | ~3.25 | Area-weighted mean of `uncertainty` field over HUC2 |
| GSFC | Leakage + noise | ~1.65 | 1-sigma = sqrt(leakage_2sig² + noise_2sig²)/2, area-weighted by mascon area |
| CSR | Inter-product spread | ~1.24 | Half-range of (max − min) across JPL/CSR/GSFC at each timestep |

## Water Year Convention

- WY N = Oct 1 of year N-1 through Sep 30 of year N
- WY2003 = Oct 2002 – Sep 2003
- GRACE missing June 2003 across all products (11 common months)

## Key Results (WY2003, model vs JPL/CSR/GSFC)

Strong correlations in western basins:
- Pacific Northwest: r = 0.996 / 0.988 / 0.987
- California: r = 0.979 / 0.980 / 0.978
- Great Basin: r = 0.931 / 0.890 / 0.863

Weaker in some eastern basins:
- Mid-Atlantic: r = 0.096 / 0.355 / 0.314
- Souris-Red-Rainy: r = 0.449 / 0.454 / 0.346

## Dependencies

```
xarray netcdf4 matplotlib numpy pandas rasterio shapely h5py
```

## Notes

- No geopandas/fiona in `subsettools` env — HUC2 boundaries loaded via `json` + `rasterio.features.geometry_mask`
- CSR time is numeric "days since 2002-01-01" (capital `Units` attr, needs manual decoding)
- GRACE timestamps are mid-month; snapped to first-of-month for alignment with ParFlow
- GSFC native mascons: 41168 irregular tiles, assigned to HUC2 via point-in-polygon (shapely)
