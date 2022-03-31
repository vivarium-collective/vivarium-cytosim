# -*- coding: utf-8 -*-

"""Top-level package for vivarium-cytosim."""

__author__ = "Blair Lyons"
__email__ = "blair208@gmail.com"
# Do not edit this string manually, always use bumpversion
# Details in CONTRIBUTING.md
__version__ = "0.0.0"


def get_module_version():
    return __version__


from .processes.cytosim import CytosimProcess  # noqa: F401
