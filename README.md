# vivarium-cytosim

A Vivarium wrapper for [Cytosim](https://gitlab.com/f-nedelec/cytosim)

---

# Installation

**Stable Release:** `pip install vivarium_cytosim` (coming soon)<br>
**Development Head:** `pip install git+https://github.com/vivarium-collective/vivarium-cytosim.git`
**Local Editable Install** `pip install -e .[dev]` (or `pip install -e .\[dev\]` on mac) from repo root directory

Or use Conda with the `env.yml` file to create the environment: 
```
conda env create -f env.yml
conda activate vivarium-models
```

### Cytosim Installation

First, clone the repo:

    git clone https://gitlab.com/f-nedelec/cytosim.git

Change the header to allow for 3D in file src/math/dim.h:

    #define DIM 3 # instead of 2

Then, make the executable (avoid the GLEW functionality):

    make sim
    make report

## Development

See [CONTRIBUTING.md](CONTRIBUTING.md) for information related to developing the code.

## Commands You Need To Know

1. `black vivarium_cytosim`

    This will fix lint issues.

2. `make build`

    This will run `tox` which will run all your tests as well as lint your code.

3. `make clean`

    This will clean up various Python and build generated files so that you can ensure that you are working in a clean environment.
