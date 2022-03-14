# vivarium-cytosim

Provides a Vivarium Process wrapper for [Cytosim](https://gitlab.com/f-nedelec/cytosim)

---

## Installation with pyenv + conda

To see all pyenv versions:

```
pyenv install list
```

To install a particular version of python (or conda):

```
pyenv install anaconda3-5.3.1
```

Install dependencies using pyenv + conda:

```
pyenv local anaconda3-5.3.1 # or whatever version you have installed
pyenv virtualenv vivarium-models
pyenv local vivarium-models
conda env update -f env.yml
```

# Installation with conda alone

Install conda: https://docs.conda.io/en/latest/miniconda.html

Using conda, you can run

```
conda env create -f env.yml
```

which will create a conda environment called `vivarium-models` with all the required dependencies (including ReaDDy) installed.

To update:

```
conda env update -f env.yml
```

### CYTOSIM

First, clone the repo:

    git clone https://gitlab.com/f-nedelec/cytosim.git

Change the header to allow for 3d (in file src/math/dim.h)

    #define DIM 3 # instead of 2

Then, make the executable (avoid the GLEW functionality):

    make sim
    make report

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
