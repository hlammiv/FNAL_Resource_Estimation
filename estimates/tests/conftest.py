import sys
from pathlib import Path

import pytest

# Make `import estimates` work when pytest is run from the repository root or from the
# paper's doc/ or doc/scripts/ directories.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

# Some tests check the paper itself: that the chapter text quotes the model's numbers, or
# that the model agrees with cross-check data kept with the paper. They read those files from
# the paper tree (the directory above the package). In the standalone code repository that
# tree is absent; those tests skip instead of erroring.
PAPER_ROOT = Path(__file__).resolve().parents[3]
PAPER_PRESENT = (PAPER_ROOT / "applications").is_dir()


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item):
    outcome = yield
    if PAPER_PRESENT:
        return
    try:
        outcome.get_result()
    except FileNotFoundError as e:
        if str(e.filename or e).startswith(str(PAPER_ROOT)):
            outcome.force_exception(
                pytest.skip.Exception("paper source not present; paper cross-check skipped"))
