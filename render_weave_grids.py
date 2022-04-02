#!/usr/bin/env python
# coding: utf-8

from itertools import chain
import logging

import numpy as np

from shapely.affinity import translate
from shapely.ops import unary_union
from shapely.wkt import dumps
from shapely.wkt import loads
from shapely.validation import make_valid
import geopandas

from weave_grids import _WeaveGrid

def centre_offset(shape, centre = (0, 0)):
  shape_c = shape.centroid.coords[0]
  return (centre[0] - shape_c[0], centre[1] - shape_c[1])


# builds the sf associate with a given weave supplied as 'loom' which is a list
# containing the coordinates in an appropriate grid (Cartesian or triangular)
# and the orderings of the strands at each coordinate location
def make_shapes_from_coded_weave_matrix(
      loom, spacing = 1, width = 1, margin = 0,
      axis1_threads = "a", axis2_threads = "b", axis3_threads = "c", crs = 3857):
  
  grid = _WeaveGrid(loom.n_axes, loom.orientations, spacing)
  ids1 = axis1_threads * (loom.dimensions[0] // len(axis1_threads))
  ids2 = axis2_threads * (loom.dimensions[1] // len(axis2_threads))
  if loom.n_axes == 3:
    ids3 = axis3_threads * (loom.dimensions[2] // len(axis3_threads))
  weave_polys = []
  bb_polys = []
  strands = []
  for coords, strand_order in zip(loom.indices, loom.orderings):
    ids = [ids1[coords[0]], ids2[coords[1]]]
    if loom.n_axes == 3:
      ids.append(ids3[coords[2]])
    cell = grid.get_grid_cell_at(coords)
    bb_polys.append(cell)
    # print(f"strand_order: {strand_order} ids: {ids}")
    if strand_order is None: continue
    if strand_order == "NA": continue
    n_slices = [len(id) for id in ids]
    next_polys = grid.get_visible_cell_strands(width, coords, strand_order, n_slices)
    weave_polys.extend(next_polys)
    labels = [list(ids[i]) for i in strand_order] # a list of lists
    labels = list(chain(*labels))                 # flatten list of lists
    # print(f"n: {len(next_polys)} labels: {labels}")
    strands.extend(labels)
  tile = unary_union(bb_polys)
  shift = centre_offset(tile)
  tile = translate(tile, shift[0], shift[1])
  return { 
    "weave_unit": make_weave_gdf(weave_polys, strands, tile, shift, margin, crs),
    "tile": geopandas.GeoDataFrame(geometry = geopandas.GeoSeries([tile])).set_crs(crs) 
    }


def make_weave_gdf(polys, strand_ids, bb, offset, margin, crs):
  weave = geopandas.GeoDataFrame(
    data = { "strand": strand_ids },
    geometry = geopandas.GeoSeries([translate(p, offset[0], offset[1]) for p in polys])
  )
  weave = weave[weave.strand != "-"]
  weave = weave.dissolve(by = "strand", as_index = False)
  weave.geometry = weave.buffer(-margin)
  return weave.clip(bb).set_crs(crs)


