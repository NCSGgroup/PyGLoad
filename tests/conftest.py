import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

# As in PyGLDA v2, the code is imported from the repository folder (`from src.basis ... import ...`); put that folder
# first on the path so that the tests run from any working directory without installation.
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

_NEEDS_LLN = ('test_consistency.py', 'test_basis.py', 'test_workflow.py')


def pytest_collection_modifyitems(config, items):
    """Tests that read the Love numbers are skipped when <data root>/input/LLN is missing (e.g. on GitHub CI)."""
    from src.config import FileTool
    try:
        if (FileTool.get_input_dir('LLN') / 'PREM-LLNs-complete.dat').exists():
            return
    except (FileNotFoundError, ValueError):      # data root not set (e.g. GitHub CI)
        pass
    skip = pytest.mark.skip(reason='needs the load Love numbers in <data root>/input/LLN (see doc/installation.md)')
    for item in items:
        if item.fspath.basename in _NEEDS_LLN:
            item.add_marker(skip)
