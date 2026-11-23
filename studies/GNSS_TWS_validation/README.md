# GNSS-TWS validation

Validate a monthly terrestrial water storage (TWS) product against GNSS: the vertical elastic loading predicted from
the product is compared with the monthly vertical displacements of 4,807 GNSS stations, 2003–2016. The study gives
per station and per region how much of the GNSS signal the product explains (WRMS reduction), how well it follows the
seasonal and non-seasonal variations, and its annual amplitude ratio and phase. It reproduces the GNSS validation of
our water storage reanalysis paper for the public products (CSR and JPL mascons, GLWS 2.0, GLWS 4.0, GLDAS-2.2), and
runs on any other monthly TWS grid.

It uses the public PyGLoad API only: `read_tenv3`, `read_ngl_steps` and `GNSSNetwork` for the GNSS cleaning,
`SphericalHarmonicBasis` and `LoveNumber` (PREM, CF frame) for the loading.

## Data

Everything lives under the PyGLoad data root (see [doc/installation.md](../../doc/installation.md)). The complete
input set of the study is distributed as a data package that *is* such a data root; with it nothing has to be
downloaded except the GRACE mascons.

| folder | content | from |
|---|---|---|
| `input/LLN/` | `PREM-LLNs-complete.dat` (Wang et al., 2012) | PyGLoad installation |
| `input/validation/` | `land_domain_0.5deg.nc` (0/1), `common_months.txt` (the 149 JPL months 2003–2016), `stations.csv`, `station_exclusions.csv`, `ice_exclusions.csv` | data package |
| `input/gnss/NGL/` | `DataHoldings.txt`, `steps.txt`, `tenv3/IGS20/`, `tenv3/IGS14/`, `cmcf/`, `ntal.cf/`, `ntol.cf/` (only to redo the cleaning) | data package, or `download.py --gnss` |
| `processed/GNSS_TWS_validation/<frame>/` | cleaned monthly series `gnss_vertical_monthly.nc` (and `.csv`) | data package, or `run.py --gnss raw` |
| `input/mascon/` | CSR RL06.3 and JPL RL06.3M CRI mascons | `download.py` |
| `input/tws_products/<name>/` | GLWS 2.0, GLWS 4.0, GLDAS-2.2 (optional) | `download.py --products ...` |
| `cache/GNSS_TWS_validation/`, `output/GNSS_TWS_validation/<frame>/` | anomaly grids, SH coefficients; scores | `run.py` |

`python download.py --list` prints every source with version, md5 of the files we used and citation;
`python download.py --check` checks the md5 of the product files present.

## Steps

```bash
cd studies/GNSS_TWS_validation
python download.py                        # CSR mascons; JPL needs a free NASA Earthdata login (instructions printed)
python run.py                             # CSR + JPL, both frames, from the processed GNSS series
python verify.py                          # PASS / FAIL against the published numbers
python run.py --tws MYTWS=/path/tws.nc --baseline CSR --products CSR    # your own product, compared with CSR
python run.py --gnss raw                  # also redo the GNSS cleaning from the NGL daily series
```

`run.py --help` lists all options (`--data-root`, `--frames`, `--products`, `--tws`, `--tws-var`, `--tws-units`,
`--product-file`, `--gnss`, `--baseline`, `--workers`).

## Method

**GNSS.** NGL daily final series, vertical, 2003.0 ≤ NGL decimal year < 2017.0, at least 100 daily solutions. Per
station a trajectory model (bias, rate, the offsets of NGL `steps.txt` (equipment and earthquakes), annual and
semi-annual terms) is fitted by weighted least squares with iterative outlier rejection (|residual| > 5 × 1.4826 MAD,
all components); bias, rate and offsets are removed, the seasonal terms and the residual are kept. This is
`GNSSNetwork.clean(remove='trend+steps', outlier=5.0)` with its additional data selection, offset detection and
quality control switched off. Calendar-month means need ≥ 10 days (`monthly(min_days=10, sigma='sem')`). Then, as
the mean over exactly the days of each monthly mean:

- frame `igs20cf` (main): IGS20 series minus the NGL `cmcf` translation (CM → CF, from the ITRF2020 geocentre motion
  model) minus ESMGFZ NTAL and NTOL in the CF frame;
- frame `igs14` (sensitivity): IGS14 series taken as CF, minus NTAL and NTOL (CF).

No common-mode filtering. The station selection (≥ 60 % of the 149 common months, ≥ 5 years from the first to the last
of them; exclusion within 300 km of Greenland, Antarctica, 18 glacier-dominated basin units and Natural Earth ice
fields ≥ 1000 km²) is given by `stations.csv` (`used`, `exclusion`); `ice_exclusions.csv` lists the ice features.

**TWS products.** Each product is put on the common months as monthly 0.5° anomalies about its own mean over its
months (a month counts if more than half of the land cells have data), on the common land domain
`land_domain_0.5deg.nc` (64,440 cells). The load is zero outside the domain (ocean, Greenland, Antarctica) and in land
cells without data. CSR is block-averaged from 0.25° to 0.5°; GLDAS-2.2 daily 0.25° fields become monthly means
(cells with ≥ 80 % of the days) and then cos-latitude weighted 0.5° means.

**Loading.** Spherical-harmonic analysis of each monthly grid to degree 360 (block-mean quadrature) and the
vertical displacement at the station coordinates of `stations.csv` (`lat`, `lon`, 5 significant digits) with
`SphericalHarmonicBasis(360).design_matrix(..., Displacement.Vertical)`, PREM load Love numbers (Wang et al., 2012)
converted to the CF frame (degree 1 after Blewitt, 2003: h₁ = −0.2595, l₁ = 0.1297), degree 0 included
(h₀ = −0.13205; the land-only load does not conserve mass).

**Scores.** Per station and product on the months both series have; both are demeaned there, then the mean and
trend of a mean + trend + annual + semi-annual least-squares fit are removed from each series (the cleaning removed
the GNSS trend):

| score | definition |
|---|---|
| `wrms_red` | 1 − WRMS(GNSS − model) / WRMS(GNSS), weights 1/σ² of the GNSS monthly means (`wrms_red_nodt`: model trend kept) |
| `r_seas` | correlation of the 12 calendar-month climatologies (each calendar month needs ≥ 2 values) |
| `r_nonseas_dt` | correlation of the series without their own climatology (`r_nonseas`: model trend kept) |
| `r_total` | correlation of the series |
| `amp_ratio` | annual amplitude model / GNSS from the fit |
| `phase_diff_d` | day of the annual maximum, model − GNSS, wrapped to ±182.6 d (+ = model later) |

Medians over stations with a station bootstrap (2000 draws, 5–95 %), globally, without Japan (24–46° N, 122–146° E;
GEONET dominates Asia and is tectonics-dominated), per WMO region and by distance to the nearest basin-unit border
(< 50 km against interior). Each (group, product, score) draws from its own stream seeded with 270 and the CRC-32 of
its name, so the ranges do not depend on which products are run together. `--baseline P` adds the share of stations
where each product beats P.

## Your own product

`--tws NAME=FILE` reads one NetCDF file with a variable of dimensions (time, lat, lon) (`--tws-var` if there are
several): CF time axis, one field per month (the first field of a calendar month is used); units mm, kg m⁻², cm or m
(attribute `units` or `--tws-units`); a regular grid of 0.5°, of a divisor of 0.5° (block means) or of a multiple of
0.5°; longitudes −180…180 or 0…360, latitudes in either order. Absolute storage or anomalies both work: the study takes
the anomaly about the product's own mean over its months within the 149 common months. Values outside the land domain
are ignored. Months missing in the product are skipped for that product.

## Expected numbers

Global medians over 4,807 stations (`verify.py` checks the WRMS reduction):

| product | WRMS red. IGS20-CF | WRMS red. IGS14 | r_seas IGS20-CF | r_nonseas_dt IGS20-CF | amp ratio IGS20-CF | phase d IGS20-CF |
|---|---|---|---|---|---|---|
| CSR RL06.3 | 0.15 | 0.06 | 0.83 | 0.28 | 0.63 | 14 |
| JPL RL06.3M | 0.13 | 0.06 | 0.83 | 0.26 | 0.57 | 16 |
| GLWS 2.0 | 0.13 | 0.04 | 0.88 | 0.20 | 0.66 | −5 |
| GLWS 4.0 | 0.17 | 0.07 | 0.87 | 0.24 | 0.67 | 5 |
| GLDAS-2.2 | 0.12 | 0.03 | 0.85 | 0.25 | 0.76 | 8 |

The data package holds the per-station reference results; `verify.py --reference` compares every station to 1e-6.

## Known limitations and notes

- The network is concentrated in North America, Europe and Japan; Japan is tectonics-dominated (report it apart).
- No common-mode filtering; barystatic ocean and ice-sheet loading are not modelled (load only on the land domain).
- NGL decimal years use 365.25-day years, while `GNSSNetwork.monthly()` assigns days to months with 365/366-day years
  (`_decyear_to_month` in `src/preprocessing/series.py`); a day near a month boundary can fall into the neighbouring
  month. The published series were made this way and the study keeps it, so the numbers reproduce.
- `LGF_Kummer.getVertical_exclude_d0` (`src/basis/green/point_load.py`) is not used by the study.
- CSR replaces its mascon file every month; the November 2026 file has the same 2003–2016 fields as the October 2024
  file used for the paper. NGL reprocesses its series and extends `steps.txt`; use the files of the data package to
  reproduce the published numbers.

## References

- Blewitt, G. (2003). Self-consistency in reference frames, geocenter definition, and surface loading of the solid
  Earth. *J. Geophys. Res.*, 108(B2), 2103. https://doi.org/10.1029/2002JB002082
- Blewitt, G., Hammond, W. C., & Kreemer, C. (2018). Harnessing the GPS data explosion for interdisciplinary science.
  *Eos*, 99. https://doi.org/10.1029/2018EO104623
- Altamimi, Z., Rebischung, P., Collilieux, X., Métivier, L., & Chanard, K. (2023). ITRF2020: an augmented reference
  frame refining the modeling of nonlinear station motions. *J. Geod.*, 97, 47. https://doi.org/10.1007/s00190-023-01738-w
- Dill, R., & Dobslaw, H. (2013). Numerical simulations of global-scale high-resolution hydrological crustal
  deformations. *J. Geophys. Res. Solid Earth*, 118, 5008–5017. https://doi.org/10.1002/jgrb.50353
- Wang, H., et al. (2012). Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and
  modified models with refined crustal structure from Crust 2.0. *Comput. Geosci.*, 49, 190–199.
  https://doi.org/10.1016/j.cageo.2012.06.022
- Product citations: `python download.py --list`.
