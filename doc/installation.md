# PyGLoad — installation guide

This guide covers the Python environment, the external data PyGLoad needs, and how to run the code. Everything
needed for the environment lives in `installation/`:

| file | purpose |
|---|---|
| `setup_env.sh` | creates and verifies the environment, registers the data root, writes the lock file |
| `environment.yml` | conda/mamba environment `pygload` |
| `environment.lock.yml` | exact versions of a tested install (written by `setup_env.sh install`) |
| `requirements.txt` | the same dependencies for `pip` |
| `pygload_paths.cfg.example` | template of the local data-root setting |

As in PyGLDA v2, only the environment is installed, not the code: the code is used directly from the repository.

---

## 1. Environment

PyGLoad has no pinned binary stack (unlike PyGLDA v2, which is locked by WaterGAP 2.2e); any recent scientific
Python works. Tested with Python 3.13, numpy 2.5, scipy 1.18, h5py 3.16, netCDF4 1.7, geopandas 1.2 and shapely 2.2.

With conda / mamba (recommended), using the setup script:

```bash
git clone https://github.com/NCSGgroup/PyGLoad.git
cd PyGLoad
bash installation/setup_env.sh check                                            # shows what would happen; installs nothing
bash installation/setup_env.sh install --data-root /path/to/PyGLoad_external_data
conda activate pygload
```

`install` creates the environment `pygload` from `environment.yml` (an existing one is kept and only verified),
writes `installation/pygload_paths.cfg` with the data root (section 2), checks that every package and the PyGLoad
code import, runs the tests and writes `installation/environment.lock.yml`. Options: `--name NAME` for another
environment name, `--no-pygmt` without the plotting library. Without the script:

```bash
mamba env create -f installation/environment.yml
conda activate pygload
```

With pip, into an existing Python:

```bash
pip install -r installation/requirements.txt
```

| Dependency | Needed for |
|---|---|
| numpy, scipy, h5py, netCDF4 | everything |
| geopandas, shapely | regional masks from shapefiles (`src/mesh/shape_mask.py`) |
| pytest | the tests in `tests/` |
| pygmt | figures in `examples/` only (optional; can live in a separate environment) |

---

## 2. External data

All data PyGLoad reads or writes live outside the code, under one folder, the **data root**, e.g.
`/media/user/My Book/Fan/PyGLoad_external_data` (cf. `PyGLDA_v2_external_data` for PyGLDA v2). Nothing is written
inside the code folder. Set the data root once, in either of two ways:

1. a local config file (recommended; ignored by git): copy `installation/pygload_paths.cfg.example` to
   `installation/pygload_paths.cfg` and set

   ```
   data_root = /path/to/PyGLoad_external_data
   ```

2. or the environment variable, which takes precedence over the file:

   ```bash
   export PYGLOAD_DATA=/path/to/PyGLoad_external_data      # e.g. in ~/.bashrc or the IDE run configuration
   ```

Without either, PyGLoad stops with an error: there is no default data folder, and a data root inside the code
folder is refused. The code folder holds code only.

The data root is organised by role:

```
<data root>/
├── input/               external sources, only read by PyGLoad
│   ├── LLN/             load Love numbers and Green's functions
│   ├── gnss/NGL/        DataHoldings.txt, steps.txt, tenv3/ (studies/GNSS_GRACE/download.py)
│   ├── land_mask/       land_mask_res0.1.h5 (removing ocean cells in src.mesh)
│   ├── mascon/          GRACE mascon NetCDF (a symbolic link is fine) and its 2-degree extract mascon.txt
│   ├── shp/             region outlines: China/, NorthAmerica/, ...
│   ├── stations/        station lists and observations, e.g. fu2015_*_digitized.txt
│   └── reference/       solutions of other software for comparison, e.g. loadvcd.txt (LoadVCD)
├── cache/               built by PyGLoad, regenerable and reused
│   ├── LGD_template*/   CD-LGF templates
│   ├── mask/            regional masks (.h5) from src.mesh
│   └── shp/             generated grid shapefiles (EU/, exp/)
├── output/              results, one folder per study
│   ├── paper_JGR2025/   data of the figures of Yang et al. (2025): fig2/, fig8/, ..., BasisStudy/
│   ├── China_basis_comparison/
│   ├── inversion_China/
│   ├── inversion_PNW_Fu2015/
│   └── GNSS_GRACE/<region>/
└── temp/                scratch, safe to delete
```

In code: `FileTool.get_input_dir()`, `get_mascon_dir()`, `get_cache_dir()`, `get_output_dir()` and
`get_temp_dir()` (each takes an optional sub-path); the older demos use the variables `INPUT`, `CACHE`, `OUTPUT`
and `TEMP`. Never write a path to data directly in a script.

| File in `input/LLN/` | Content | Used by |
|---|---|---|
| `PREM-LLNs-complete.dat`, `iasp91-...`, `ak135-...`, `PREMhard-...`, `PREMsoft-...`, `iasp91hard-...`, `ak135hard-...` | load Love numbers, degrees 1 to 46 340, CE frame (Wang et al., 2012) | `LLN_Data.PREM`, ... |
| `Wang.dat` | load Love numbers to degree 13 000, CE frame | `LLN_Data.Wang` (default) |
| `REF_6371_loading_love_numbers.mat` | REF 6371 Love numbers, CM frame (Bevis et al., 2016) | `LLN_Data.REF` |
| `PREM-LGFs.dat` | tabulated PREM Green's functions (Wang et al., 2012) | one test, some demos |

These files are third-party data and are not covered by the PyGLoad license; please cite their sources.

`installation/pygload_paths.cfg` is ignored by git; it is the only local file in the code folder.

---

## 3. Running PyGLoad

The code is not installed. As in PyGLDA v2, it is used directly from the repository: the core code is the folder
`src/`, and every study, example and test puts the repository folder on the Python path in its first lines,

```python
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))   # repository root (parents[1] for examples/)
from src.basis.green.gf_basis import GreenFunctionBasis
import src as pgw                                               # main classes: pgw.LoveNumber, pgw.SlepianBasis, ...
```

so scripts and tests run from any folder, in the `pygload` environment:

```bash
python examples/slepian_truncation.py          # one feature per example
python studies/Fu2015_PNW/real_data.py         # a reproducible study
python -m pytest tests -q                      # all tests; those needing the Love numbers are skipped if input/LLN/ is missing
```

Your own scripts outside the repository add the same line with the path of the repository folder. In PyCharm, open
the repository folder as the project and select the `pygload` conda interpreter; `from src...` imports resolve
without any further setting.

### What is published

Everything in the repository except personal and local files: `src/`, the studies `GNSS_GRACE/`, `Fu2015_PNW/` and
`China/`, `examples/`, `tests/`, `doc/` and `installation/`. Not published (git-ignored):
`studies/JGR2025_basis_functions/`, `studies/GNSS_TWS_validation/`, `doc/cowork_sync.md`, `Claude outputs/` and
`installation/pygload_paths.cfg`.

See [user_guide.md](user_guide.md) for all options and usage.
