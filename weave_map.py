#!/usr/bin/env python
# coding: utf-8

from itertools import chain
from dataclasses import dataclass
import logging

import numpy as np
import geopandas 
from shapely.affinity import translate
from shapely.geometry import Point

from triaxial_weave_units import get_triaxial_weave_unit
from biaxial_weave_units import get_biaxial_weave_unit


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
  parameter_info(margin, aspect)
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

def parameter_info(margin, aspect):
  max_margin = (1 - aspect) / 2
  if aspect == 0:
    logging.info("""Setting aspect to 0 is probably not a great plan.""")
  if aspect < 0 or aspect > 1:
    logging.warning("""Values of aspect outside the range 0 to 1 won't produce
                    tiles that will look like weaves, but they might be 
                    pretty anyway! Values less than -1 seem particularly
                    promising, especially with opacity set less than 1.""")
  if margin > max_margin:
    logging.warning(f"""With aspect set to {aspect} the largest margin that will 
                    work is {max_margin}. Lower values are required to produce 
                    proper tileable weaves. Higher values will make nice tilings 
                    that are not strictly weaves. An alternative is to make the 
                    tile with margin = 0, and then apply a negative buffer after 
                    you have tiled your map.""")
  


# returns width, height, x-origin and y-origin of GeoSeries
def _get_width_height(gs:geopandas.GeoSeries):
  extent = gs.total_bounds
  return extent[2] - extent[0], extent[3] - extent[1], extent[0], extent[1]


def _get_grid(ll, nums, tdim):
  return np.array(
    np.meshgrid(np.arange(nums[0]) * tdim[0] + ll[0],
                np.arange(nums[1]) * tdim[1] + ll[1])).reshape(2, nums[0] * nums[1]).transpose()
      

# takes tile and to_tile GeoSeries and makes a list of translation vectors
# both geoseries are expected to be orthogonally arranged so that the centres
# are a rectangular grid
def _get_rect_centres(tile_gs:geopandas.GeoSeries, to_tile_gs:geopandas.GeoSeries):
  to_tile_w, to_tile_h, to_tile_x0, to_tile_y0  = _get_width_height(to_tile_gs)
  tile_w, tile_h, tile_x0, tile_y0 = _get_width_height(tile_gs)
  nx = int(np.ceil(to_tile_w / tile_w))
  ny = int(np.ceil(to_tile_h / tile_h))
  x0 = ((nx * tile_w) - to_tile_w) / 2 + tile_x0 + to_tile_x0
  y0 = ((ny * tile_h) - to_tile_h) / 2 + tile_y0 + to_tile_y0
  return _get_grid((x0, y0), (nx, ny), (tile_w, tile_h))


def _get_hex_centres(tile_gs:geopandas.GeoSeries, to_tile_gs:geopandas.GeoSeries):
  to_tile_w, to_tile_h, to_tile_x0, to_tile_y0  = _get_width_height(to_tile_gs)
  tile_w, tile_h, tile_x0, tile_y0 = _get_width_height(tile_gs)
  nx = int(np.ceil(to_tile_w / (tile_w * 3 / 2))) + 1
  ny = int(np.ceil(to_tile_h / tile_h)) + 1
  x0 = ((nx * tile_w * 3 / 2) - to_tile_w) / 2 + tile_x0 + to_tile_x0
  y0 = ((ny * tile_h) - to_tile_h) / 2 + tile_y0 + to_tile_y0
  g1 = _get_grid((x0, y0 + tile_h / 4), (nx, ny), (tile_w * 3 / 2, tile_h))
  g2 = _get_grid((x0 + tile_w * 3 / 4, y0 - tile_h / 4), (nx, ny), (tile_w * 3 / 2, tile_h))
  return np.append(g1, g2).reshape((g1.shape[0] + g2.shape[0], 2))


# translates the geometries in a GeoSeries and returns a list of geometries
def translate_geoms(gs:geopandas.GeoSeries, dxdy:tuple = (0, 0)):
  return [translate(s, dxdy[0], dxdy[1]) for s in gs]


def rotate_gdf_to_geoseries(gdf:geopandas.GeoDataFrame, angle:float, about:tuple = None):
  centre = gdf.geometry.unary_union.centroid.coords[0] if about is None else about
  return gdf.geometry.rotate(angle, origin = centre), centre


# takes WeaveUnit and region GeoDataFrame and tiles the units over the region
def get_tiling(unit: WeaveUnit, region:geopandas.GeoDataFrame, rotation = 0):
  if region.geometry.name != "geometry":
    region.rename_geometry("geometry", inplace = True)

  tile_gs = unit.tile.geometry
  if rotation != 0:
    to_tile_gs, centre = rotate_gdf_to_geoseries(region, rotation)
  else:
    to_tile_gs = region.geometry # might need to think about geometry columns with other names...
    
  if unit.type in ("hex", "cube"):
    shifts = _get_hex_centres(tile_gs, to_tile_gs)
  else:
    shifts = _get_rect_centres(tile_gs, to_tile_gs)

  # TODO: filter the shifts by some buffered version of the region, so as to save 
  # time (possibly) on generating and then clipping the tiled weave units
  
  tiles = chain(*[translate_geoms(unit.weave_unit.geometry, dxdy) for dxdy in shifts])
  ids = list(unit.weave_unit.strand) * shifts.shape[0]
  
  tiles_gs = geopandas.GeoSeries(tiles)
  if rotation != 0:
    tiles_gs = tiles_gs.rotate(-rotation, origin = centre)

  return geopandas.GeoDataFrame(data = {"strand": ids}, crs = unit.weave_unit.crs,
                                geometry = tiles_gs)