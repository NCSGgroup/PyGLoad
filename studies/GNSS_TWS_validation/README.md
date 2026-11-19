# GNSS-TWS validation

Monthly vertical GNSS displacements at 4,807 stations (2003–2016) compared with the elastic loading predicted from
global terrestrial water storage (TWS) products. It measures how much of the observed seasonal and non-seasonal
vertical motion a TWS product explains. The products are CSR and JPL GRACE mascons, GLWS 2.0, GLWS 4.0 and
GLDAS-2.2.

The study uses only the public `pygwater` API: the GNSS trajectory model and monthly means (`GNSSNetwork`), the
spherical-harmonic basis (`SphericalHarmonicBasis`) and the PREM load Love numbers in the CF frame (`LoveNumber`).

## Data

Everything lives under the PyGWater data root (see [doc/installation.md](../../doc/installation.md)):

| folder | content | source |
|---|---|---|
| `input/GNSS_TWS_validation/` | the study's reference files: station list and exclusions, cleaned monthly GNSS series, common land domain and months, reference scores | study data set (see its `README.md`) |
| `input/tws_products/<name>/` | the TWS products | `download.py` |
| `input/gnss/NGL/` | NGL daily series, steps file, ESMGFZ load predictions (only to rebuild the GNSS series) | `download.py --gnss` |
| `cache/GNSS_TWS_validation/` | product anomaly grids and spherical-harmonic coefficients | `run.py` |
| `output/GNSS_TWS_validation/` | scores and summaries | `run.py`, `verify.py` |

`python download.py --list` prints every source with its version, citation and DOI. CSR, GLWS 2.0 and GLWS 4.0 are
open downloads. JPL and GLDAS-2.2 need a free NASA Earthdata login. Fetch them by hand and put them in
`input/tws_products/JPL/` and `input/tws_products/GLDAS22/`. If a product is already stored somewhere else, point to
it with `GNSS_TWS_<NAME>_PATH` instead of copying it (see `config.py`).

## Steps

```bash
pip install pandas                                  # needed by the study, in addition to PyGWater's requirements
cd studies/GNSS_TWS_validation
python download.py                                  # products (about 3.5 GB, plus JPL / GLDAS-2.2 by hand)
python run.py                                       # all products; GNSS_TWS_PRODUCTS=CSR,JPL for a subset
python verify.py                                    # compare with the reference scores
```

`run.py` does the following:

1. **GNSS**: reads the cleaned monthly series of the reference files.
2. **Products**: puts each product on the common months (the 149 JPL months of 2003–2016) as 0.5° anomalies about
   its own mean, on the common land domain (WaterGAP land cells without Greenland and Antarctica, 64,440 cells). The
   load is zero outside the domain. Cells without data are set to 0. CSR is block-averaged from 0.25° to 0.5°.
   GLDAS-2.2 daily 0.25° fields become monthly means (≥ 80 % of the days) and then cos-latitude weighted 0.5° means.
   Each grid then goes through spherical-harmonic analysis to degree 360.
3. **Forward model**: computes the vertical displacement at the stations with the SH design matrix and PREM load Love
   numbers (Wang et al., 2012) in the CF frame, degree 0 included.
4. **Scores**: computes per station, on the months both series share, `r_total`, `r_seas`, `r_nonseas(_dt)`, the
   annual amplitude ratio and phase difference, and the WRMS reduction. It also gives medians with station-bootstrap
   5–95 % ranges (2000 draws, seed 270). These are reported globally, without Japan (GEONET dominates Asia), per WMO
   region, and by distance to the nearest basin-unit border (< 50 km against interior). The run writes
   `gnss_scores_per_station.csv`, `gnss_summary.csv` and `table_gnss.csv`.

Runtime on a desktop (4 threads): about 20 s per product for the SH analysis (GLDAS-2.2: about 15 min to read the
5,000 daily files), 1–2 min for the forward model and about 10 s for the scores.

## Expected numbers

The global medians over 4,807 stations (`results/table_gnss_public.csv` of the reference files):

| product | WRMS reduction % | r_total | amplitude error factor | WRMS red. % without Japan |
|---|---|---|---|---|
| CSR RL06.3 | 6.1 | 0.43 | 1.80 | 10.7 |
| JPL RL06.3M | 5.7 | 0.43 | 1.80 | 10.7 |
| GLWS 2.0 | 3.5 | 0.40 | 1.70 | 6.9 |
| GLWS 4.0 | 7.2 | 0.46 | 1.54 | 10.9 |
| GLDAS-2.2 | 2.9 | 0.41 | 2.02 | 6.4 |

The amplitude error factor is exp(median |ln(amplitude ratio)|). `verify.py` checks every median and station count
of the reference summary to rounding. The bootstrap bounds can move, because the random draws depend on which
products are scored together. The typical change is 0.001. It reaches 0.04 in small regions and 2 days for the
phase.

The CSR mascon file on the CSR server is replaced every month. The reference numbers were made with the
October 2024 file. The November 2026 file gives the same scores, because its 2003–2016 fields are unchanged.

## Rebuilding the GNSS series

```bash
python download.py --gnss                           # NGL IGS14 series, steps.txt, ESMGFZ NTAL / NTOL (CF), ~7 GB
GNSS_TWS_GNSS=rebuild python run.py                 # output files get the suffix _rebuild
python verify.py _rebuild
```

This rebuilds the cleaned series from the daily data with the settings of the reference series:

- daily epochs 2003-01-01 to 2016-12-31;
- `GNSSNetwork.clean(remove='trend+steps', outlier=5.0)` with the NGL offsets, and the extra data selection and
  quality control of `clean()` switched off;
- monthly means of months with ≥ 10 days (`sigma='sem'`);
- ESMGFZ NTAL + NTOL subtracted as the mean over the same days.

The station selection (≥ 60 % of the common months, ≥ 5 years) and the exclusions near Greenland, Antarctica,
glacier-dominated basins and large ice fields come from the station list. `run.py --default-cleaning` with
`GNSS_TWS_GNSS=rebuild` uses `clean()` with its default data selection, offset detection and quality control
instead. The output suffix is then `_default_cleaning`.

NGL keeps adding offsets to `steps.txt`. A newer steps file can therefore change a few series.

## References

- Blewitt, G., Hammond, W. C., & Kreemer, C. (2018). Harnessing the GPS data explosion for interdisciplinary science.
  *Eos*, 99. https://doi.org/10.1029/2018EO104623
- Dill, R., & Dobslaw, H. (2013). Numerical simulations of global-scale high-resolution hydrological crustal
  deformations. *J. Geophys. Res. Solid Earth*, 118, 5008–5017. https://doi.org/10.1002/jgrb.50353
- Wang, H., et al. (2012). Load Love numbers and Green's functions for elastic Earth models PREM, iasp91, ak135, and
  modified models with refined crustal structure from Crust 2.0. *Comput. Geosci.*, 49, 190–199.
  https://doi.org/10.1016/j.cageo.2012.06.022
- Product citations: `python download.py --list`.
