import sys

import pytest

from pytest_cython import PytestCython

sys.exit(pytest.main(['--benchmark-skip'], plugins = [PytestCython()]))
