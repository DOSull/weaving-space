#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
from triaxial_weave_units import get_triaxial_weave_unit
from biaxial_weave_units import get_biaxial_weave_unit

import numpy as np
import geopandas 


@dataclass
class WeaveUnit:
    weave_unit: geopandas.GeoDataFrame
    transform: np.ndarray
    tile: geopandas.GeoDataFrame
    type: str
    
    def __init__(self, weave_unit, transform, tile, type):
        self.weave_unit = weave_unit
        self.transform = transform
        self.tile = tile
        self.type = type


def get_weave_unit(weave_type = "plain", spacing = 10000, aspect = 1,
                   margin = 0, n = (2, 2), strands = "a|b|c",
                   tie_up = None, tr = None, th = None, crs = 3857):
  
  if weave_type in ("hex", "cube"):
    unit = get_triaxial_weave_unit(spacing = spacing, aspect = aspect,
                                   margin = margin, weave_type = weave_type,
                                   strands = strands, crs = crs)
  else:
    unit = get_biaxial_weave_unit(spacing = spacing, aspect = aspect,
                                   margin = margin, weave_type = weave_type,
                                   strands = strands, crs = crs)

  return WeaveUnit(
      unit["weave_unit"],
      np.diag(np.ones(2)),
      unit["tile"],
      weave_type
  )

