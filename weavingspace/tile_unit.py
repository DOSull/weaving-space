#!/usr/bin/env python
# coding: utf-8

"""The `TileUnit` subclass of `weavingspace.tileable.Tileable` implements
many 'conventional' tilings of the plane.

Examples:
  A `TileUnit` is initialised like this

    tile_unit = TileUnit(tiling_type = "cairo")

  The `tiling_type` may be one of the following

  + "cairo" the Cairo tiling more formally known as the Laves
  [3<sup>2</sup>.4.3.4] tiling. The author's favourite tiling, hence it
  has its own tiling_type.
  + "hex-slice" a range of dissections of the regular hexagon into,
  2, 3, 4, 6, or 12 'pie slices'. The number of slices is set by
  specifying an additional argument `n`. Slices are cut either starting
  at the corners of  the hexagon or from the midpoints of hexagon edges,
  by specifying an additional argument `offset` set to either
  0 or 1 respectively.
  + "hex-dissection" a range of 4, 7 or 9-fold dissections of the hexagon.  
  + "laves" a range of isohedral tilings. See [this article](https://en.wikipedia.org/wiki/List_of_Euclidean_uniform_tilings#Laves_tilings).
  The desired tiling is specified by the additional argument `code` which
  is a string like "3.3.4.3.4".
  + "archimedean" a range of tilings by regular polygons. See [this
  article](https://en.wikipedia.org/wiki/Euclidean_tilings_by_convex_regular_polygons#Archimedean,_uniform_or_semiregular_tilings). Many of these are the dual tilings of
  the Laves tilings. The desired tiling is specified by the additional
  argument `code` which is a string like "3.3.4.3.4". Not all the
  possible Archimedean tilings are implemented.
  + "hex-colouring" three colourings of the regular hexagon tiling, of
  either 3, 4, or 7 colours, as specified by the argument `n`.
  + "square-colouring" one colouring of the regular square tiling, of 5
  colours as specified by the argument `n = 5`.

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

import copy
from dataclasses import dataclass
from typing import Iterable
import string

import geopandas as gpd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine

# from weavingspace.tileable import Tileable
# from weavingspace.tileable import TileShape

# import weavingspace.tiling_utils as tiling_utils
# import weavingspace.tiling_geometries as tiling_geometries

from weavingspace import Tileable
from weavingspace import TileShape

from weavingspace import tiling_utils
from weavingspace import tiling_geometries


@dataclass
class TileUnit(Tileable):
  """Class to represent the tiles of a 'conventional' tiling.
  """
  tiling_type:str = ""
  """tiling type as detailed in the class documentation preamble."""
  offset:int = 1
  """offset for 'hex-dissection' and 'hex-slice tilings. Defaults to 1."""
  n:int = 3
  """number of dissections or colours in 'hex-dissection', 'hex-slice' and
  'hex-colouring' tilings. Defaults to 3."""
  code:str = "3.3.4.3.4"
  """the code for 'laves' or 'archimedean' tiling types."""

  def __init__(self, **kwargs) -> None:
    super(TileUnit, self).__init__(**kwargs)
    # this next line makes all TileUnit geometries shapely 2.x precision-aware
    self.tiles.geometry = tiling_utils.gridify(self.tiles.geometry)
    if not self.tiling_type is None:
      self.tiling_type = self.tiling_type.lower()
    if self.base_shape == TileShape.TRIANGLE:
      self._modify_tile()
      self._modify_tiles()
      self.setup_vectors()
    self.setup_regularised_prototile_from_tiles()
    # if self.regularised_prototile is None:
    #   self.setup_regularised_prototile_from_tiles()


  def _setup_tiles(self) -> None:
    """Delegates setup of the unit to various functions depending
    on self.tiling_type.
    """
    if self.tiling_type == "cairo":
      tiling_geometries.setup_cairo(self)
    elif self.tiling_type == "hex-slice":
      tiling_geometries.setup_hex_slice(self)
    elif self.tiling_type == "hex-dissection":
      tiling_geometries.setup_hex_dissection(self)
    elif "cross" in self.tiling_type:
      tiling_geometries.setup_crosses(self)
    elif self.tiling_type == "square-slice":
      tiling_geometries.setup_square_slice(self)
    elif self.tiling_type == "laves":
      tiling_geometries.setup_laves(self)
    elif self.tiling_type == "archimedean":
      tiling_geometries.setup_archimedean(self)
    elif self.tiling_type in ("hex-colouring", "hex-coloring"):
      tiling_geometries.setup_hex_colouring(self)
    elif self.tiling_type in ("square-colouring", "square-coloring"):
      tiling_geometries.setup_square_colouring(self)
    else:
      tiling_geometries._setup_none_tile(self)
    return


  def _setup_regularised_prototile(self) -> None:
    self.regularise_tiles()
    self.regularised_prototile.geometry = tiling_utils.repair_polygon(
      self.regularised_prototile.geometry)
    

  def _modify_tiles(self) -> None:
    """It is not trivial to tile a triangle, so this function augments
    the tiles of a triangular tile to a diamond by 180 degree
    rotation. Operation is 'in place'.
    """
    tiles = self.tiles.geometry
    ids = list(self.tiles.tile_id)

    new_ids = list(string.ascii_letters[:(len(ids) * 2)])
    tiles = tiles.translate(0, -tiles.total_bounds[1])
    twins = [affine.rotate(tile, a, origin = (0, 0))
             for tile in tiles
             for a in range(0, 360, 180)]
    self.tiles = gpd.GeoDataFrame(
      data = {"tile_id": new_ids},
      geometry = tiling_utils.gridify(gpd.GeoSeries(twins)),
      crs = self.tiles.crs)
    return None


  def _modify_tile(self) -> None:
    """It is not trivial to tile a triangular tile so this function
    changes the tile to a diamond by manually altering the tile in place
    to be a diamond shape.
    """
    tile = self.prototile.loc[0, "geometry"]
    # translate to sit on x-axis
    tile = affine.translate(tile, 0, -tile.bounds[1])
    pts = [p for p in tile.exterior.coords]
    pts[-1] = (pts[1][0], -pts[1][1])
    self.prototile.geometry = tiling_utils.gridify(
      gpd.GeoSeries([geom.Polygon(pts)], crs = self.crs))
    self.base_shape = TileShape.DIAMOND
    return None


  def _get_legend_key_shapes(self, polygon:geom.Polygon,
                 counts:Iterable = [1] * 25, angle:float = 0,
                 radial:bool = False) -> list[geom.Polygon]:
    """Returns a set of shapes that can be used to make a legend key
    symbol for the supplied polygon. In TileUnit this is a set of 'nested'
    polygons.

    Args:
      polygon (geom.Polygon): the polygon to symbolise.
      count (Iterable, optional): iterable of the counts of each slice.
        Defaults to [1] * 25.
      rot (float, optional): rotation that may have to be applied.
        Not used in the TileUnit case. Defaults to 0.

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
      radius = tiling_utils.get_collapse_distance(polygon)
      distances = distances * radius / distances[-1]
      nested_polys = [polygon.buffer(-d, join_style = 2, cap_style = 3) 
                      for d in distances]
      # return converted to annuli (who knows someone might set alpha < 1)
      nested_polys = [g1.difference(g2) for g1, g2 in
              zip(nested_polys[:-1], nested_polys[1:])]
      return [p for c, p in zip(counts, nested_polys) if c > 0]
    else:
      n = sum(counts)
      slice_posns = list(np.cumsum(counts))
      slice_posns = [0] + [p / n for p in slice_posns]
      return [tiling_utils.get_polygon_sector(polygon, i, j)
          for i, j in zip(slice_posns[:-1], slice_posns[1:])]


  # Note that geopandas clip is not order preserving hence we do this
  # one polygon at a time...
  def inset_prototile(self, d:float = 0) -> "TileUnit":
    """Returns a new TileUnit clipped by `self.regularised_tile` after
    a negative buffer d has been applied.

    Args:
      d (float, optional): the inset distance. Defaults to 0.

    Returns:
      TileUnit: the new TileUnit with inset applied.
    """
    inset_tile = \
      self.regularised_prototile.loc[0, "geometry"].buffer(-d, join_style = 2, cap_style = 3)
    new_tiles = [inset_tile.intersection(e) for e in self.tiles.geometry]
    result = copy.deepcopy(self)
    result.tiles.geometry = gpd.GeoSeries(new_tiles)
    return result


  def scale_tiles(self, sf:float = 1) -> "TileUnit":
    """Scales the tiles by the specified factor, centred on (0, 0).

    Args:
      sf (float, optional): scale factor to apply. Defaults to 1.

    Returns:
      TileUnit: the scaled TileUnit.
    """
    result = copy.deepcopy(self)
    result.tiles.geometry = tiling_utils.gridify(
      self.tiles.geometry.scale(sf, sf, origin = (0, 0)))
    return result
  
  
  def _as_circles(self) -> "TileUnit":
    result = copy.deepcopy(self)
    result.tiles.geometry = gpd.GeoSeries([tiling_utils.in_circle(p) 
                                           for p in self.tiles.geometry])
    return result


