import runpy


class PytestCython:
    def pytest_sessionstart(self):
        runpy.run_module('build_dev')
