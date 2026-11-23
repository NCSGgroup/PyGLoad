# PyGWater

**Surface loading and GNSS inversion of terrestrial water storage, with a choice of basis functions.**

PyGWater computes how the solid Earth deforms under surface water loads and inverts GNSS displacements for the water
storage that caused them. Spherical harmonics, point-mass and disk-load Green's functions (D-LGF, CD-LGF) and Slepian
functions share one interface, so every computation can be repeated with another basis (Yang et al., 2025).

## What it does

- **Forward modelling**: GNSS displacement (up, north, east), geoid, gravity and InSAR line of sight from a load.
- **GNSS preprocessing**: download of NGL station series, offsets, outliers, data selection, monthly means and
  quality control, following Schumacher et al. (2018).
- **Inversion**: regularised monthly inversion of a GNSS network for water storage, with the regularisation chosen
  automatically.
- **Comparison with GRACE**: regional water storage and mass from GNSS against GRACE mascons.

## Getting started

1. Create the Python environment from the files in [installation/](installation) and install PyGWater.
2. Choose a data folder outside the code and register it in the paths file in [installation/](installation).
3. Put the load Love numbers in the input folder of that data folder.
4. Run the tests to check the setup.

Step-by-step instructions are in the [installation guide](doc/installation.md).

## Where to go next

| To ... | See |
|---|---|
| learn the concepts and all options | [user guide](doc/user_guide.md) |
| try the basis functions on small cases | [examples/](examples) |
| invert GNSS for water storage and compare with GRACE | [studies/GNSS_GRACE/](studies/GNSS_GRACE) |
| reproduce published results | [studies/](studies) |

## Repository layout

| Folder | Content |
|---|---|
| [src/](src) | the core code, imported as `PyGWater` |
| [examples/](examples) | short demonstrations |
| [studies/](studies) | complete, reproducible studies |
| [tests/](tests) | automated tests |
| [doc/](doc) | installation guide and user guide |
| [installation/](installation) | environment and data-path settings |

Data never live in the code folder: inputs, caches and results all go to the data folder chosen at installation.

## Citation

> Yang, F., Shen, Y., Zhang, W., & Forootan, E. (2025). On the choice of basis functions for modeling Earth's elastic
> deformation due to surface loading. *Journal of Geophysical Research: Solid Earth*, 130, e2025JB031662.
> https://doi.org/10.1029/2025JB031662

When you use GNSS series from the Nevada Geodetic Laboratory, please also cite Blewitt, Hammond & Kreemer (2018),
*Eos*, 99. See [CITATION.cff](CITATION.cff).

PyGWater belongs to the GravityPython toolboxes, with [PyHawk](https://github.com/NCSGgroup/PyHawk),
[SaGEA](https://github.com/NCSGgroup/SaGEA), [PyGLDA](https://github.com/AAUGeodesyGroup/PyGLDA) and
[SAGEA-fluid](https://github.com/NCSGgroup/SAGEA-fluid).

## License

[MIT](LICENSE). The load Love numbers (Wang et al., 2012; Bevis et al., 2016) are third-party data.

Contact: Fan Yang, Aalborg University — [fany@plan.aau.dk](mailto:fany@plan.aau.dk)
