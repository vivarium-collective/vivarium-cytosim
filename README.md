# vivarium-models

[![Build Status](https://github.com/allen-cell-animated/vivarium_models/workflows/Build%20Main/badge.svg)](https://github.com/allen-cell-animated/vivarium_models/actions)
[![Documentation](https://github.com/allen-cell-animated/vivarium_models/workflows/Documentation/badge.svg)](https://allen-cell-animated.github.io/vivarium_models/)
[![Code Coverage](https://codecov.io/gh/allen-cell-animated/vivarium_models/branch/main/graph/badge.svg)](https://codecov.io/gh/allen-cell-animated/vivarium_models)

Simularium prototypes of connecting models in Vivarium

---

## Installation

Install conda: https://docs.conda.io/en/latest/miniconda.html

Using conda, you can run `conda env create -f env.yml`, which will create a conda environment called `vivarium_models` with all the required dependencies (including ReaDDy) installed.

### Alternatively:

**Stable Release:** `pip install vivarium_models`<br>
**Development Head:** `pip install git+https://github.com/allen-cell-animated/vivarium-models.git`

ReaDDy models depend on ReaDDy, which requires conda. Install ReaDDy with `conda install -c readdy/label/dev readdy` after adding the conda-forge channel `conda config --add channels conda-forge`

## Documentation

For full package documentation please visit [allen-cell-animated.github.io/vivarium_models](https://allen-cell-animated.github.io/vivarium_models).

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for information related to developing the code.

## The Four Commands You Need To Know

1. `pip install -e .[dev]`

    This will install your package in editable mode with all the required development
    dependencies (i.e. `tox`).

2. `make build`

    This will run `tox` which will run all your tests in both Python 3.7
    and Python 3.8 as well as linting your code.

3. `make clean`

    This will clean up various Python and build generated files so that you can ensure
    that you are working in a clean environment.

4. `make docs`

    This will generate and launch a web browser to view the most up-to-date
    documentation for your Python package.


**Allen Institute Software License**

