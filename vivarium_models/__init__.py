# -*- coding: utf-8 -*-

"""Top-level package for Vivarium Models."""

__author__ = "Blair Lyons"
__email__ = "blairl@alleninstitute.org"
# Do not edit this string manually, always use bumpversion
# Details in CONTRIBUTING.md
__version__ = "0.0.0"


def get_module_version():
    return __version__


from .processes import ReaddyActinProcess  # noqa: F401
