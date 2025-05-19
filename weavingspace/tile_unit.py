#!/usr/bin/env python
# coding: utf-8

"""The `TileUnit` subclass of `weavingspace.tileable.Tileable` implements many
'conventional' tilings of the plane.

Examples:
  A `TileUnit` is initialised like this

    tile_unit = TileUnit(tiling_type = "cairo")

  The `tiling_type` may be one of the following

  + "cairo" the Cairo tiling more formally known as the Laves
  [3<sup>2</sup>.4.3.4] tiling. The author's favourite tiling, hence it has its
  own tiling_type.
  + "hex-slice" a range of dissections of the regular hexagon into, 2, 3, 4, 6,
  or 12 'pie slices'. The number of slices is set by specifying an additional
  argument `n`. Slices are cut either starting at the corners of  the hexagon
  or from the midpoints of hexagon edges, by specifying an additional argument
  `offset` set to either 0 or 1 respectively.
  + "hex-dissection" a range of 4, 7 or 9-fold dissections of the hexagon.  
  + "laves" a range of isohedral tilings. See [this article](https://en.wikipedia.org/wiki/List_of_Euclidean_uniform_tilings#Laves_tilings).
  The desired tiling is specified by the additional argument `code` which is a
  string like "3.3.4.3.4". Not all the possible Laves tilings are implemented.
  + "archimedean" a range of tilings by regular polygons. See [this
  article](https://en.wikipedia.org/wiki/Euclidean_tilings_by_convex_regular_polygons#Archimedean,_uniform_or_semiregular_tilings). 
  Many of these are the dual tilings of the Laves tilings. The desired tiling
  is specified by the additional argument `code` which is a string like
  "3.3.4.3.4". Not all the possible Archimedean tilings are implemented.
  + "hex-colouring" or "square-colouring" colourings of the regular hexagonal
  and square tilings of between 2 and 10 colours, as specified by the argument
  `n`.
  + "hex-slice" or "square-slice" dissections of the regular hexagonal and
  square tilings of between 2 and 12 colours, as specified by the arguments `n`
  and `offset`.
  + "crosses" colourings of cross-shaped pentominoes of between 

  Spacing and coordinate reference of the tile unit are specified by the
  `weavingspace.tileable.Tileable` superclass variables
  `weavingspace.tileable.Tileable.spacing` and
  `weavingspace.tileable.Tileable.crs`.

  Base tilings by squares, hexagons or triangles can also be requested
  using

    tile_unit = TileUnit()  # square tiling, the default
    tile_unit = TileUnit(tile_shape = TileShape.HEXAGON)
    tile_unit = TileUnit(tile_shape = TileShape.TRIANGLE)

  The first two of these have only one tile_id value, and so cannot be
  used for multivariate mapping. The triangle case has two tile_id
  values so may be useful in its base form.

  To create custom tilings start from one of the base tiles above, and
  explicitly set the `weavingspace.tileable.Tileable.tiles` variable
  by geometric construction of suitable shapely.geometry.Polygons. TODO: A detailed example of this usage can be found here ....
"""

from __future__ import annotations
import copy
from dataclasses import dataclass
from typing import Iterable

import geopandas as gpd
import numpy as np
import shapely.geometry as geom

from weavingspace import Tileable
from weavingspace import tiling_utils
from weavingspace import tiling_geometries


@dataclass
class TileUnit(Tileable):
  """Class to represent the tiles of a 'conventional' tiling. Most of the
  functionality of TileUnit is either in `the weavingspace.Tileable` superclass
  or, for setup, in the functions in `weavingspace.tiling_geometries`.
  """
  tiling_type:str = ""
  """tiling type as detailed in the class documentation preamble."""
  offset:float = 1
  """offset for 'dissection' and 'slice' tilings. Defaults to 1."""
  n:int = 3
  """number of dissections or colours in 'hex-dissection', 'hex-slice' and
  'hex-colouring' tilings. Defaults to 3."""
  code:str = "3.3.4.3.4"
  """the code for 'laves' or 'archimedean' tiling types."""

  def __init__(self, **kwargs) -> None:
    # pass the kwargs to the superclass constructor
    # it will delegate setting up the tiles back to `TileUnit._setup_tiles()`
    super(TileUnit, self).__init__(**kwargs)


  def _setup_tiles(self) -> None|str:
    """Delegates setup of unit to functions in `tiling_geometries` depending on
    `self.tiling_type`. If there is a problem a string is returned. If all is
    OK then None is returned (some care is required to do this!)

    The logical tests applied to check for tiling type are expansive to allow
    scope for user errors.

    Returns:
      str|None: if a problem occurs a message string, otherwise None.
    """
    if self.tiling_type == "cairo":
      return tiling_geometries.setup_cairo(self)
    elif "hex" in self.tiling_type:
      if "slice" in self.tiling_type:
        return tiling_geometries.setup_hex_slice(self)
      elif "dissect" in self.tiling_type:
        return tiling_geometries.setup_hex_dissection(self)
      elif "col" in self.tiling_type:
        return tiling_geometries.setup_hex_colouring(self)
    elif "square" in self.tiling_type:
      if "slice" in self.tiling_type:
        return tiling_geometries.setup_square_slice(self)
      elif "col" in self.tiling_type:
        return tiling_geometries.setup_square_colouring(self)
      elif "trisect" in self.tiling_type:
        return tiling_geometries.setup_square_trisection(self)
      elif "dissect" in self.tiling_type:
        return tiling_geometries.setup_square_dissection(self)
    elif "lave" in self.tiling_type:
      return tiling_geometries.setup_laves(self)
    elif "archi" in self.tiling_type:
      # easy to misspell this one!
      return tiling_geometries.setup_archimedean(self)
    elif "cross" in self.tiling_type:
      return tiling_geometries.setup_crosses(self)
    elif "star1" in self.tiling_type:
      return tiling_geometries.setup_star_polygon_1(self)
    elif "star2" in self.tiling_type:
      return tiling_geometries.setup_star_polygon_2(self)
    elif "chavey" in self.tiling_type:
      return tiling_geometries.setup_chavey(self)
    else:
      return tiling_geometries.setup_base_tiling(self)


  def _setup_regularised_prototile(self, override = False) -> None:
    """Sets up a 'regularised prototile' which fully contains all the tile
    elements, i.e. it does not cut them. In all TileUnit cases a suitable shape
    is the union of the elements.
    """
    # For whatever reasons in the web app version the unioning operation as
    # operated by tiling_utils.safe_union() is anything but and produces
    # TopologyException crashes... so here is a safe_union avoidant way...
    if self.regularised_prototile is None or override:
      tiles = copy.deepcopy(self.tiles.geometry)
      tiles = gpd.GeoSeries(
        [tiling_utils.gridify(p.buffer(
          tiling_utils.RESOLUTION * 10, join_style=2, cap_style=3))
          for p in tiles])
      self.regularised_prototile = gpd.GeoDataFrame(
        geometry = gpd.GeoSeries(
          [tiles.union_all().buffer(-tiling_utils.RESOLUTION * 10, 
                                    join_style=2, cap_style=3).simplify(
                                      tiling_utils.RESOLUTION * 10)]),
        crs = self.crs)
    return None


  def _get_legend_key_shapes(
      self, 
      polygon:geom.Polygon,
      counts:Iterable = [1] * 25, 
      angle:float = 0,
      radial:bool = False) -> list[geom.Polygon]:
    """Returns a set of shapes that can be used to make a legend key symbol for
    the supplied polygon. In TileUnit this is a set of 'nested' polygons formed
    by serially negative buffering the tile shapes.

    Args:
      polygon (geom.Polygon): the polygon to symbolise.
      count (Iterable, optional): iterable of the counts of each slice.
        Defaults to [1] * 25.
      rot (float, optional): rotation that may have to be applied. Not used in 
        the TileUnit case. Defaults to 0.
      radial (bool): if True then a pie slice dissection will be applied; if
        False then a set of 'nested' shapes will be applied. Defaults to False.

    Returns:
      list[geom.Polygon]: a list of nested polygons.
    """
    if not radial:
      n = sum(counts)
      # bandwidths = list(np.cumsum(counts))
      bandwidths = [c / n for c in counts]
      bandwidths = [bw if bw > 0.05 or bw == 0 else 0.05
              for bw in bandwidths]
      n = sum(bandwidths)
      bandwidths = [0] + [bw / n for bw in bandwidths]
      # # make buffer widths that will yield approx equal area 'annuli'
      # bandwidths = range(n_steps + 1)
      # sqrt exaggerates outermost annuli, which can otherwise disappear
      bandwidths = [np.sqrt(bw) for bw in bandwidths]
      distances = np.cumsum(bandwidths)
      # get the negative buffer distance that will 'collapse' the polygon
      # radius = tiling_utils.get_collapse_distance(polygon)
      radius = tiling_utils.get_apothem_length(polygon)
      distances = distances * radius / distances[-1]
      nested_polys = [polygon.buffer(-d, join_style = 2, cap_style = 3)
                      for d in distances]
      # DON'T CONVERT TO ANNULI - it washes the final colours out in rendering
      # nested_polys = [g1.difference(g2) for g1, g2 in
      #         zip(nested_polys[:-1], nested_polys[1:])]
      return [p for c, p in zip(counts, nested_polys) if c > 0]
    else:
      n = sum(counts)
      slice_posns = list(np.cumsum(counts))
      slice_posns = [0] + [p / n for p in slice_posns]
      return [tiling_utils.get_polygon_sector(polygon, i, j)
              for i, j in zip(slice_posns[:-1], slice_posns[1:])]


  def inset_prototile(self, d:float = 0) -> "TileUnit":
    """Returns a new TileUnit clipped by `self.regularised_prototile` after
    a negative buffer d has been applied. Note that geopandas clip is not order
    preserving hence we do this one polygon at a time.

    Args:
      d (float, optional): the inset distance. Defaults to 0.

    Returns:
      TileUnit: the new TileUnit with inset applied.
    """
    if d == 0: return self
    inset_tile = \
      self.regularised_prototile.loc[0, "geometry"] \
        .buffer(-d, join_style = 2, cap_style = 3)
    new_tiles = [tiling_utils.get_clean_polygon(inset_tile.intersection(e))
                 for e in self.tiles.geometry]
    result = copy.deepcopy(self)
    result.tiles.geometry = gpd.GeoSeries(new_tiles)
    return result
  
  
  def _as_circles(self) -> "TileUnit":
    """Experimental implementation of returning tiles as their incircles.

    Returns:
        TileUnit: a tiling which replaces each tile with its incircle.
    """
    result = copy.deepcopy(self)
    result.tiles.geometry = gpd.GeoSeries([tiling_utils.in_circle(p) 
                                           for p in self.tiles.geometry])
    return result


