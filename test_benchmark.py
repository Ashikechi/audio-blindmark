import sys

import pytest

from pytest_cython import PytestCython

sys.exit(pytest.main([
    '--benchmark-only',
    '--benchmark-timer=time.process_time',
    '--benchmark-sort=mean',
    '--benchmark-cprofile=tottime',
    '--benchmark-cprofile-dump=prof/benchmark'
], plugins = [PytestCython()]))
