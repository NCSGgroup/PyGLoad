import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# The repository folder is the package PyGWater. Normally it is installed (`pip install -e .`); if not, register the
# folder under that name so that the tests also run from a plain checkout (e.g. on GitHub CI).
try:
    import PyGWater  # noqa: F401
except ImportError:
    spec = importlib.util.spec_from_file_location('PyGWater', ROOT / '__init__.py',
                                                  submodule_search_locations=[str(ROOT)])
    module = importlib.util.module_from_spec(spec)
    sys.modules['PyGWater'] = module
    spec.loader.exec_module(module)

_NEEDS_LLN = ('test_consistency.py', 'test_basis.py', 'test_workflow.py')


def pytest_collection_modifyitems(config, items):
    """Tests that read the Love numbers are skipped when <data root>/input/LLN is missing (e.g. on GitHub CI)."""
    from PyGWater.config import FileTool
    try:
        if (FileTool.get_input_dir('LLN') / 'PREM-LLNs-complete.dat').exists():
            return
    except (FileNotFoundError, ValueError):      # data root not set (e.g. GitHub CI)
        pass
    skip = pytest.mark.skip(reason='needs the load Love numbers in <data root>/input/LLN (see doc/installation.md)')
    for item in items:
        if item.fspath.basename in _NEEDS_LLN:
            item.add_marker(skip)
