**WIP: This is a repository of a research work in progress. All the code are not stable and can not be used in production.**
# Installation
You should install [Poetry](https://python-poetry.org/) first.
```sh
poetry install
```
For developers, you can also install with test tools.
```sh
poetry install --with dev
git config core.hooksPath .githooks
```
# Packaging
```sh
python -m build
```
Then the sdist and wheel will be in `dist/`
# Compiling Cython extensions
```sh
python build_dev.py
```
# Testing
You must compile all Cython extensions first
```sh
python test.py
```
or run benchmark test with profiling
```sh
python test_benchmark.py
```
# Clean
```sh
python clean.py
```
