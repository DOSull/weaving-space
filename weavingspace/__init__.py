#!/usr/bin/env python
# coding: utf-8

"""
See [this notebook](https://github.com/DOSull/weaving-space/blob/main/using-the-library.ipynb) for some preliminary usage notes.
"""

# this is a workaround that allows Tileable to check if an instance
# is a WeaveUnit without circular import problems
def is_weaveunit(x):
  return isinstance(x, WeaveUnit)

from .tiling_utils import *
from ._loom import *
from ._weave_grid import *
from .tileable import *
from .tiling_geometries import *
from .tile_unit import *
from .weave_matrices import *
from .weave_unit import *
from .tile_map import *
from .symmetry import *
from .topology import *
