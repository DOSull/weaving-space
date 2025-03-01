#!/usr/bin/env python
# coding: utf-8

"""Implements `weavingspace.tileable.TileShape` and
`weavingspace.tileable.Tileable` the base classes for
`weavingspace.tile_unit.TileUnit` and `weavingspace.weave_unit.WeaveUnit`.

`Tileable` should not be called directly, but is instead accessed from the
`weavingspace.tile_unit.TileUnit` or `weavingspace.weave_unit.WeaveUnit`
constructor.

Several methods of `weavingspace.tileable.Tileable` are generally useful and
can be accessed through its subclasses.
"""

from enum import Enum
from typing import Union
from dataclasses import dataclass
import copy

import matplotlib.pyplot as pyplot

import geopandas as gpd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine

# import weavingspace.tiling_utils as tiling_utils
from weavingspace import tiling_utils


class TileShape(Enum):
  """The available base tile shapes.

  NOTE: the TRIANGLE type does not persist, but should be converted to a
  DIAMOND or HEXAGON type during `Tileable` construction.
  """

  RECTANGLE = "rectangle"
  HEXAGON = "hexagon"
  TRIANGLE = "triangle"
  DIAMOND = "diamond"


@dataclass
class Tileable:
  """Class to represent a tileable set of tile geometries.
  """

  tiles: gpd.GeoDataFrame = None
  """the geometries with associated `title_id` attribute encoding their 
  different colouring."""
  prototile: gpd.GeoDataFrame = None
  """the tileable polygon (rectangle, hexagon or diamond)"""
  spacing: float = 1000.0
  """the tile spacing effectively the resolution of the tiling. Defaults to
  1000"""
  base_shape: TileShape = TileShape.RECTANGLE
  """the tile shape. Defaults to 'RECTANGLE'"""
  vectors: dict[tuple[int], tuple[float]] = None
  """translation vector symmetries of the tiling"""
  regularised_prototile: gpd.GeoDataFrame = None
  """polygon containing the tiles of this tileable, usually a union of its
  tile polygons"""
  crs: int = 3857
  """coordinate reference system of the tile. Most often an ESPG code but
  any valid geopandas CRS specification is valid. Defaults to 3857 (i.e. Web
  Mercator)."""
  rotation: float = 0.0
  """cumulative rotation of the tileable."""
  debug: bool = False
  """if True prints debug messages. Defaults to False."""

  # Tileable constructor called by subclasses - should not be used directly
  def __init__(self, **kwargs):
    for k, v in kwargs.items():
      self.__dict__[k] = v
    if self.debug:
      print(
        f"""Debugging messages enabled for Tileable (but there aren't
          any at the moment...)"""
      )
    self._setup_tiles()
    self.setup_vectors()
    self._setup_regularised_prototile()
    return


  def setup_vectors(self) -> None:
    """Sets up the symmetry translation vectors as floating point pairs
    indexed by integer tuples with respect to either a rectangular or
    triangular grid location.

    Derived from the size and shape of the tile attribute. These are not
    the minimal translation vectors, but the 'face to face' vectors of the
    tile, such that a hexagonal tile will have 3 vectors, not the minimal
    parallelogram pair. Also supplies the inverse vectors.

    The vectors are stored in a dictionary indexed by their
    coordinates, e.g.

      {( 1,  0): ( 100, 0), ( 0,  1): (0,  100),
       (-1,  0): (-100, 0), ( 0, -1): (0, -100)}

    For a tileable of type `TileShape.HEXAGON`, the indexing tuples
    have three components. See https://www.redblobgames.com/grids/hexagons/
    """
    t = self.prototile.loc[0, "geometry"]
    pts = [p for p in t.exterior.coords][:-1]
    n_pts = len(pts)
    vec_dict = {}
    if n_pts == 4:
      vecs = [(q[0] - p[0], q[1] - p[1])
          for p, q in zip(pts, pts[1:] + pts[:1])]
      i = [1, 0, -1,  0]
      j = [0, 1,  0, -1]
      vec_dict = {(i, j): v for i, j, v in zip(i, j, vecs)}
    elif n_pts == 6:
      vecs = [(q[0] - p[0], q[1] - p[1])
          for p, q in zip(pts, pts[2:] + pts[:2])]
      # hex grid coordinates associated with each of the vectors
      i = [ 0,  1,  1,  0, -1, -1]
      j = [ 1,  0, -1, -1,  0,  1]
      k = [-1, -1,  0,  1,  1,  0]
      vec_dict = {(i, j, k): v for i, j, k, v in zip(i, j, k, vecs)}
    self.vectors = vec_dict


  def get_vectors(
      self, as_dict: bool = False
    ) -> Union[dict[tuple[int], tuple[float]], list[tuple[float]]]:
    """
    Returns symmetry translation vectors as floating point pairs.
    Optionally returns the vectors in a dictionary indexed by their
    coordinates, e.g.

      {( 1,  0): ( 100, 0), ( 0,  1): (0,  100),
       (-1,  0): (-100, 0), ( 0, -1): (0, -100)}

    Returns:
      Union[ dict[tuple[int],tuple[float]], list[tuple[float]] ]:
        either the vectors as a list of float tuples, or a dictionary
        of those vectors indexed by integer coordinate tuples.
    """
    if as_dict:
      return self.vectors
    else:
      return list(self.vectors.values())


  # Make up a regularised tile by carefully unioning the tiles
  def setup_regularised_prototile_from_tiles(self) -> None:
    """Sets the regularised tile to a union of the tiles."""
    self.regularised_prototile = copy.deepcopy(self.prototile)
    self.regularised_prototile.geometry = [tiling_utils.safe_union(
      self.tiles.geometry, as_polygon = True)]
    # This simplification seems very crude but fixes all kinds of issues...
    # particularly with the triaxial weave units... where intersection 
    # operations are prone to creating spurious vertices, etc.
    self.regularised_prototile.loc[0, 'geometry'] = \
      self.regularised_prototile.loc[0, 'geometry'].simplify(
        self.spacing * tiling_utils.RESOLUTION)
    # self.regularised_prototile.geometry[0] = \
    #   self.regularised_prototile.geometry[0].simplify(
    #     self.spacing * tiling_utils.RESOLUTION)
    return


  def merge_fragments(self, fragments:list[geom.Polygon]) -> list[geom.Polygon]:
    """
    Merges a set of polygons based on testing if they touch when subjected
    to the translation vectors provided by `get_vectors()`.

    Called by `regularise_tiles()` to combine tiles in a tile unit that
    may be fragmented as supplied but will combine after tiling into single
    tiles. This step makes for more efficient implementation of the
    tiling of map regions.

    Args:
      fragments (list[geom.Polygon]): A set of polygons to merge.

    Returns:
      list[geom.Polygon]: A minimal list of merged polygons.
    """
    if len(fragments) == 1:
      return [f for f in fragments if not f.is_empty]
    fragments = [f for f in fragments if not f.is_empty]
    prototile = self.prototile.loc[0, "geometry"]
    reg_prototile = copy.deepcopy(self.regularised_prototile.loc[0, "geometry"])
    changes_made = True
    while changes_made:
      changes_made = False
      for v in self.vectors.values():
        # empty list to collect the new fragments
        # assembled in this iteration
        next_frags = []
        t_frags = [affine.translate(f, v[0], v[1]) for f in fragments]
        # build a set of any near matching pairs of
        # fragments and their translated copies
        matches = set()
        for i, f1 in enumerate(fragments):
          for j, f2, in enumerate(t_frags):
            if i < j and tiling_utils.touch_along_an_edge(f1, f2):
              matches.add((i, j))
        # determine which of these when unioned has the larger area in common # with the prototile
        frags_to_remove = set()
        for i, j in matches:
          f1, f2 = fragments[i], t_frags[j]
          u1 = f1.buffer(tiling_utils.RESOLUTION, join_style = 2, cap_style = 3).union(
            f2.buffer(tiling_utils.RESOLUTION, join_style = 2, cap_style = 3))
          u2 = affine.translate(u1, -v[0], -v[1])
          if prototile.intersection(u1).area > prototile.intersection(u2).area:
            u1 = u1.buffer(-tiling_utils.RESOLUTION, join_style = 2, cap_style = 3)
            u2 = u2.buffer(-tiling_utils.RESOLUTION, join_style = 2, cap_style = 3)
            next_frags.append(u1)
            reg_prototile = reg_prototile.union(u1).difference(u2)
          else:
            u1 = u1.buffer(-tiling_utils.RESOLUTION, join_style = 2, cap_style = 3)
            u2 = u2.buffer(-tiling_utils.RESOLUTION, join_style = 2, cap_style = 3)
            next_frags.append(u2)
            reg_prototile = reg_prototile.union(u2).difference(u1)
          changes_made = True
          frags_to_remove.add(i)
          frags_to_remove.add(j)
        fragments = [f for i, f in enumerate(fragments)
                     if not (i in frags_to_remove)] + next_frags
    self.regularised_prototile.loc[0, "geometry"] = reg_prototile
    # self.regularised_prototile.geometry[0] = reg_prototile
    return [f for f in fragments if not f.is_empty] # don't return any duds


  def reattach_tiles(self) -> None:
    """Move tiles that are outside the regularised prototile main polygon
    back inside it adjusting regularised prototile if needed.
    """
    reg_prototile = self.regularised_prototile.loc[0, "geometry"]
    new_reg_prototile = copy.deepcopy(reg_prototile)
    new_tiles = list(self.tiles.geometry)
    for i, p in enumerate(self.tiles.geometry):
      if np.isclose(reg_prototile.intersection(p).area, p.area):
        new_tiles[i] = p
        continue
      for v in self.vectors.values():
        t_p = affine.translate(p, v[0], v[1])
        if reg_prototile.intersects(t_p):
          new_reg_prototile = new_reg_prototile.union(t_p)
          new_tiles[i] = t_p
    self.tiles.geometry = gpd.GeoSeries(new_tiles)
    self.regularised_prototile.loc[0, "geometry"] = new_reg_prototile
    # self.regularised_prototile.geometry[0] = new_reg_prototile
    return None


  def regularise_tiles(self) -> None:
    """Combines separate tiles that share a tile_id value into
    single tiles, if they would end up touching after tiling.

    Also adjusts the `Tileable.regularised_prototile`
    attribute accordingly.
    """
    self.regularised_prototile = copy.deepcopy(self.prototile)
    # This preserves order while finding uniques, unlike list(set()).
    # Reordering ids might cause confusion when colour palettes
    # are not assigned explicitly to each id, but in the order
    # encountered in the tile_id Series of the GeoDataFrame.
    tiles, tile_ids = [], []
    ids = list(self.tiles.tile_id.unique())
    for id in ids:
      fragment_set = list(
        self.tiles[self.tiles.tile_id == id].geometry)
      merge_result = self.merge_fragments(fragment_set)
      tiles.extend(merge_result)
      tile_ids.extend([id] * len(merge_result))

    self.tiles = gpd.GeoDataFrame(
      data = {"tile_id": tile_ids},
      crs = self.crs,
      geometry = gpd.GeoSeries([tiling_utils.get_clean_polygon(t) 
                                for t in tiles]))

    self.regularised_prototile = \
      self.regularised_prototile.explode(ignore_index = True)
    if self.regularised_prototile.shape[0] > 1:
      self.regularised_prototile.geometry = tiling_utils.get_largest_polygon(
        self.regularised_prototile.geometry)
    return None


  def get_local_patch(self, r: int = 1,
                      include_0: bool = False) -> gpd.GeoDataFrame:
    """Returns a GeoDataFrame with translated copies of the Tileable.

    The geodataframe takes the same form as the `Tileable.tile` attribute.

    Args:
      r (int, optional): the number of 'layers' out from the unit to
        which the translate copies will extendt. Defaults to `1`.
      include_0 (bool, optional): If True includes the Tileable itself at
        (0, 0). Defaults to `False`.

    Returns:
      gpd.GeoDataFrame: A GeoDataframe of the tiles extended by a number 
        of 'layers'.
    """
    # a dictionary of all the vectors we need, starting with (0, 0)
    vecs = (
      {(0, 0, 0): (0, 0)}
      if self.base_shape in (TileShape.HEXAGON,)
      else {(0, 0): (0, 0)}
    )
    steps = r if self.base_shape in (TileShape.HEXAGON,) else r * 2
    # a dictionary of the last 'layer' of added vectors
    last_vecs = copy.deepcopy(vecs)
    # get the translation vectors in a dictionary indexed by coordinates
    # we keep track of the sum of vectors using the (integer) coordinates
    # to avoid duplication of moves due to floating point inaccuracies
    vectors = self.get_vectors(as_dict = True)
    for i in range(steps):
      new_vecs = {}
      for k1, v1 in last_vecs.items():
        for k2, v2 in vectors.items():
          # add the coordinates to make a new key...
          new_key = tuple([k1[i] + k2[i] for i in range(len(k1))])
          # ... and the vector components to make a new value
          new_val = (v1[0] + v2[0], v1[1] + v2[1])
          # if we haven't reached here before store it
          if not new_key in vecs:
            new_vecs[new_key] = new_val
      # extend the vectors and set the last layer to the set just added
      vecs = vecs | new_vecs
      last_vecs = new_vecs
    if not include_0:  # throw away the identity vector
      vecs.pop((0, 0, 0) if self.base_shape in (TileShape.HEXAGON,) else (0, 0))
    ids, tiles = [], []
    # we need to add the translated prototiles in order of their distance from 
    # tile 0, esp. in the square case, i.e. something like this:
    #
    #      5 4 3 4 5
    #      4 2 1 2 4
    #      3 1 0 1 3
    #      4 2 1 2 4
    #      5 4 3 4 5
    #
    # this is important for topology detection, where filtering back to the
    # local patch of radius 1 is greatly eased if prototiles have been added in 
    # this order. We use the vector index tuples not the euclidean distances
    # because this may be more resistant to odd effects for non-convex tiles
    extent = self.prototile.geometry.scale(
      2 * r + tiling_utils.RESOLUTION, 2 * r + tiling_utils.RESOLUTION,
      origin = self.prototile.loc[0, "geometry"].centroid)[0]
    vector_lengths = {index: np.sqrt(sum([_ ** 2 for _ in index]))
                      for index in vecs.keys()}
    ordered_vector_keys = [k for k, v in sorted(vector_lengths.items(), 
                                                key = lambda item: item[1])]
    for k in ordered_vector_keys:
      v = vecs[k]
      if geom.Point(v[0], v[1]).within(extent):
        ids.extend(self.tiles.tile_id)
        tiles.extend(
          self.tiles.geometry.apply(affine.translate, xoff = v[0], yoff = v[1]))
    return gpd.GeoDataFrame(
      data = {"tile_id": ids}, crs=self.crs,
      # gridifying here causes an error in some cases
      geometry = gpd.GeoSeries(tiles) #tiling_utils.gridify(gpd.GeoSeries(tiles))
    )


  def fit_tiles_to_prototile(self, centre_tile: int = 0) -> None:
    """Fits the tiles so they sit inside the prototile boundary.

    If tiles project outside the boundaries of the prototile, this
    method will clip them so that they don't. This may result in
    'fragmented' tiles, i.e. pieces that would form a single tile
    after tiling which are separated into fragments.

    Args:
      centre_tile (int, optional): the index position of the central
        tile. Defaults to `0`.
    """
    dxy = self.tiles.geometry[centre_tile].centroid
    self.tiles.geometry = self.tiles.translate(-dxy.x, -dxy.y)
    # use r = 2 because rectangular tiles may need diagonal neighbours
    patch = (
      self.get_local_patch(r=2, include_0=True)
      if self.base_shape in (TileShape.RECTANGLE,)
      else self.get_local_patch(r=1, include_0=True)
    )
    self.tiles = patch.clip(self.prototile)
    # repair any weirdness...
    self.tiles.geometry = tiling_utils.repair_polygon(self.tiles.geometry)
    self.tiles = self.tiles[self.tiles.geometry.area > 0]
    self.regularised_prototile = copy.deepcopy(self.prototile)
    return None


  # applicable to both TileUnits and WeaveUnits
  def inset_tiles(self, inset: float = 0) -> "Tileable":
    """Returns a new Tileable with an inset applied around the tiles.

    Works by applying a negative buffer of specfied size to all tiles.
    Tiles that collapse to zero area are removed and the tile_id
    attribute updated accordingly.

    NOTE: this method is likely to not preserve the relative area of tiles.

    Args:
      inset (float, optional): The distance to inset. Defaults to `0`.

    Returns:
      "Tileable": the new inset Tileable.
    """
    inset_tiles, inset_ids = [], []
    for p, id in zip(self.tiles.geometry, self.tiles.tile_id):
      b = p.buffer(-inset, join_style = 2, cap_style = 3)
      if not b.area <= 0:
        inset_tiles.append(b)
        inset_ids.append(id)
    result = copy.deepcopy(self)
    result.tiles = gpd.GeoDataFrame(
      data={"tile_id": inset_ids},
      crs=self.crs,
      geometry=gpd.GeoSeries(inset_tiles),
    )
    return result


  def plot(self, ax = None, show_prototile: bool = True, 
    show_reg_prototile: bool = True, show_ids: str = "tile_id",
    show_vectors: bool = False, r: int = 0, prototile_edgecolour: str = "k", 
    reg_prototile_edgecolour: str = "r", r_alpha: float = 0.3, 
    cmap: list[str] = None, figsize: tuple[float] = (8, 8), **kwargs) -> pyplot.axes:
    """Plots a representation of the Tileable on the supplied axis. **kwargs
    are passed on to matplotlib.plot()

    Args:
      ax (_type_, optional): matplotlib axis to draw to. Defaults to None.
      show_prototile (bool, optional): if `True` show the tile outline.
        Defaults to `True`.
      show_reg_prototile (bool, optional): if `True` show the regularised tile
        outline. Defaults to `True`.
      show_ids (str, optional): if `tile_id` show the tile_ids. If
        `id` show index number. If None or `''` don't label tiles.
        Defaults to `tile_id`.
      show_vectors (bool, optional): if `True` show the translation
        vectors (not the minimal pair, but those used by
        `get_local_patch()`). Defaults to `False`.
      r (int, optional): passed to `get_local_patch()` to show context if
        greater than 0. Defaults to `0`.
      r_alpha (float, optional): alpha setting for units other than the
        central one. Defaults to 0.3.
      prototile_edgecolour (str, optional): outline colour for the tile.
        Defaults to `"k"`.
      reg_prototile_edgecolour (str, optional): outline colour for the
        regularised. Defaults to `"r"`.
      cmap (list[str], optional): colour map to apply to the central
        tiles. Defaults to `None`.
      figsize (tuple[float], optional): size of the figure.
        Defaults to `(8, 8)`.
    
    Returns:
      pyplot.axes: to which calling context may add things.
    """
    w = self.prototile.loc[0, "geometry"].bounds[2] - \
      self.prototile.loc[0, "geometry"].bounds[0]
    n_cols = len(set(self.tiles.tile_id))
    if cmap is None:
      cm = "Dark2" if n_cols <= 8 else "Paired"
    else:
      cm = cmap
    if ax is None:
      ax = self.tiles.plot(
        column="tile_id", cmap=cm, figsize=figsize, **kwargs)
    else:
      self.tiles.plot(
        ax=ax, column="tile_id", cmap=cm, figsize=figsize, **kwargs)
    if show_ids != None and show_ids != "":
      do_label = True
      if show_ids == "tile_id" or show_ids == True:
        labels = self.tiles.tile_id
      elif show_ids == "id":
        labels = [str(i) for i in range(self.tiles.shape[0])]
      else:
        do_label = False
      if do_label:
        for id, tile in zip(labels, self.tiles.geometry):
          ax.annotate(id, (tile.centroid.x, tile.centroid.y),
            ha = "center", va = "center", bbox = {"lw": 0, "fc": "#ffffff40"})
    if r > 0:
      self.get_local_patch(r=r).plot(
        ax = ax, column = "tile_id", alpha = r_alpha, cmap = cm, **kwargs)
    if show_prototile:
      self.prototile.plot(ax = ax, ec = prototile_edgecolour, lw = 0.5, 
                          fc = "#00000000", **kwargs)
    if show_vectors:  # note that arrows in mpl are dimensioned in plotspace
      vecs = self.get_vectors()
      for v in vecs[: len(vecs) // 2]:
        ax.arrow(0, 0, v[0], v[1], color = "k", width = w * 0.002,
          head_width = w * 0.05, length_includes_head = True, zorder = 3)
    if show_reg_prototile:
      self.regularised_prototile.plot(
        ax = ax, ec = reg_prototile_edgecolour, fc = "#00000000", 
        lw = 1.5, zorder = 2, **kwargs)
    return ax


  def _get_legend_tiles(self):
    """Returns the tiles augmented by a rotation column.

    This base implementation may be overridden by specific tile unit types.
    In particular see
    `weavingspace.weave_unit.WeaveUnit._get_legend_tiles()`.
    """
    tiles = copy.deepcopy(self.tiles)
    tiles["rotation"] = 0
    return tiles


  def transform_scale(self, xscale: float = 1.0, yscale: float = 1.0) -> "Tileable":
    """Transforms tileable by scaling.

    Args:
      xscale (float, optional): x scale factor. Defaults to 1.0.
      yscale (float, optional): y scale factor. Defaults to 1.0.

    Returns:
      Tileable: the transformed Tileable.
    """
    result = copy.deepcopy(self)
    result.tiles.geometry = tiling_utils.gridify(
      self.tiles.geometry.scale(xscale, yscale, origin=(0, 0)))
    result.prototile.geometry = tiling_utils.gridify(
      self.prototile.geometry.scale(xscale, yscale, origin=(0, 0)))
    result.regularised_prototile.geometry = tiling_utils.gridify(
      self.regularised_prototile.geometry.scale(xscale, yscale, origin=(0, 0)))
    result.setup_vectors()
    return result


  def transform_rotate(self, angle: float = 0.0) -> "Tileable":
    """Transforms tiling by rotation.

    Args:
      angle (float, optional): angle to rotate by. Defaults to 0.0.

    Returns:
      Tileable: the transformed Tileable.
    """
    result = copy.deepcopy(self)
    result.tiles.geometry = tiling_utils.gridify(
      self.tiles.geometry.rotate(angle, origin=(0, 0)))
    result.prototile.geometry = tiling_utils.gridify(
      self.prototile.geometry.rotate(angle, origin=(0, 0)))
    result.regularised_prototile.geometry = tiling_utils.gridify(
      self.regularised_prototile.geometry.rotate(angle, origin=(0, 0)))
    result.setup_vectors()
    result.rotation = result.rotation + angle
    return result


  def transform_skew(self, xa: float = 0.0, ya: float = 0.0) -> "Tileable":
    """Transforms tiling by skewing

    Args:
      xa (float, optional): x direction skew. Defaults to 0.0.
      ya (float, optional): y direction skew. Defaults to 0.0.

    Returns:
      Tileable: the transformed Tileable.
    """
    result = copy.deepcopy(self)
    result.tiles.geometry = tiling_utils.gridify(
      self.tiles.geometry.skew(xa, ya, origin=(0, 0)))
    result.prototile.geometry = tiling_utils.gridify(
      self.prototile.geometry.skew(xa, ya, origin=(0, 0)))
    result.regularised_prototile.geometry = tiling_utils.gridify(
      self.regularised_prototile.geometry.skew(xa, ya, origin=(0, 0)))
    result.setup_vectors()
    return result
