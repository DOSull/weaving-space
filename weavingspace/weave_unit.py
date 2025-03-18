#!/usr/bin/env python
# coding: utf-8

"""The `WeaveUnit` subclass of `weavingspace.tileable.Tileable` implements
tileable geometric patterns constructed by specifying 2- or 3-axial weaves.

Examples:
  Explain usage here...
"""

import itertools
from copy import deepcopy
import logging
from dataclasses import dataclass
from collections import defaultdict
from typing import Iterable, Union

import pandas as pd
import geopandas as gpd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine
import shapely.ops

from weavingspace import weave_matrices
from weavingspace import tiling_utils

from weavingspace import Loom
from weavingspace import WeaveGrid

from weavingspace import TileShape
from weavingspace import Tileable


@dataclass
class WeaveUnit(Tileable):
  """Extends Tileable to allow for tiles that appear like woven patterns.
  """

  weave_type:str = "plain"
  """type of weave pattern, one of `plain`, `twill`, `basket`, `cube`, `hex` or
  `this`. Defaults to `plain`."""
  aspect:float = 1.
  """width of strands relative to the `spacing`. Defaults to 1.0."""
  n:Union[int, tuple[int]] = (2, 2)
  """number of over-under strands in biaxial weaves. Only one item is 
  required in a plain weave. Twill and basket patterns expect an even number of 
  entries in the tuple."""
  strands:str = "a|b|c"
  """specification of the strand labels along each axis. Defaults to `a|b|c`."""
  _tie_up:np.ndarray = None
  """optional tie-up array to pass through for `this` weave type."""
  _tr:np.ndarray = None
  """optional treadling array to pass through for `this` weave type."""
  _th:np.ndarray = None
  """optional threading array to pass through for `this` weave type."""

  def __init__(self, **kwargs):
    super(WeaveUnit, self).__init__(**kwargs)
    self.weave_type = self.weave_type.lower()


  def _setup_tiles(self) -> None:
    """Returns dictionary with weave unit and tile GeoDataFrames based on
    parameters already supplied to the constructor.
    """
    self._parameter_info()

    if self.weave_type in ("hex", "cube"):
      self.base_shape = TileShape.HEXAGON
      self._setup_triaxial_weave_unit()
      # this is experimental for now:
      # self._minimise_triaxial_weave_unit()
    else:
      self.base_shape = TileShape.RECTANGLE
      self._setup_biaxial_weave_unit()
    return
  
  
  def _setup_regularised_prototile(self) -> None:
    self.regularise_tiles()
    if self.aspect < 1:
      self.reattach_tiles()
    self.regularised_prototile.geometry = tiling_utils.repair_polygon(
      self.regularised_prototile.geometry)
    

  def _parameter_info(self) -> None:
    """Outputs logging message concerning the supplied aspect settings.
    """
    if self.aspect == 0:
      logging.info("Setting aspect to 0 is probably not a great plan.")

    if self.aspect < 0 or self.aspect > 1:
      logging.warning(
        """Values of aspect outside the range 0 to 1 won't produce tiles 
        that will look like weaves, but they might be pretty anyway! Values
        less than -1 seem particularly promising, especially with opacity 
        set less than 1.""")

    return None


  def _setup_biaxial_weave_unit(self) -> None:
    """Returns weave tiles GeoDataFrame and tile GeoDataFrame in a 
    dictionary based on paramters already supplied to constructor.
    """
    warp_threads, weft_threads, _ = \
      tiling_utils.get_strand_ids(self.strands)

    if self.weave_type == "basket" and isinstance(self.n, (list, tuple)):
      self.n = self.n[0]

    p = weave_matrices.get_weave_pattern_matrix(
      weave_type = self.weave_type, n = self.n, warp = warp_threads,
      weft = weft_threads, tie_up = self._tie_up, tr = self._tr,
      th = self._th)

    self._make_shapes_from_coded_weave_matrix(
      Loom(p), strand_labels = [weft_threads, warp_threads, []])


  def _get_triaxial_weave_matrices(self,
    strands_1:Union[list[str],tuple[str]] = ["a"],
    strands_2:Union[list[str],tuple[str]] = ["b"],
    strands_3:Union[list[str],tuple[str]] = ["c"]) -> Loom:
    """Returns encoded weave pattern matrix as Loom of three biaxial matrices.

    Allowed weave_types: "cube" or "hex".

    "hex" is not flexible and will fail with any strand label lists that are 
    not length 3 or include more than one non-blank "-" item. You can 
    generate the "hex" weave with the default settings in any case!

    Strand lists should be length 3 or length 1. "cube" tolerates more 
    than "hex" for the items in the strand lists.

    Defaults will produce 'mad weave'.

    Args:
      strands_1 (Union[list[str],tuple[str]], optional): list of labels
      for warp strands. Defaults to ["a"].
      strands_2 (Union[list[str],tuple[str]], optional): list of labels
      for weft strands. Defaults to ["b"].
      strands_3 (Union[list[str],tuple[str]], optional): list of labels
      for weft strands. Defaults to ["c"].

    Returns:
      Loom: which combines the three biaxial weaves 12, 23 and 31 implied
      by the strand label lists.
    """
    if self.weave_type == "hex":
      loom = Loom(
        weave_matrices.get_weave_pattern_matrix(
          weave_type = "this", tie_up = np.ones((6, 6)),
          warp = strands_1, weft = strands_2),
        weave_matrices.get_weave_pattern_matrix(
          weave_type = "this", tie_up = np.ones((6, 6)),
          warp = strands_2, weft = strands_3),
        weave_matrices.get_weave_pattern_matrix(
          weave_type = "this", tie_up = np.ones((6, 6)),
          warp = strands_3, weft = strands_1))
    else: # "cube"
      loom = Loom(
      # Note n = (1,2,1,2) is required here to force 6x6 twill
        weave_matrices.get_weave_pattern_matrix(
          weave_type = "twill", n = (1, 2, 1, 2),
          warp = strands_1, weft = strands_2),
        weave_matrices.get_weave_pattern_matrix(
          weave_type = "twill", n = (1, 2, 1, 2),
          warp = strands_2, weft = strands_3),
        weave_matrices.get_weave_pattern_matrix(
          weave_type = "twill", n = (1, 2, 1, 2),
          warp = strands_3, weft = strands_1))
    return loom


  def _setup_triaxial_weave_unit(self) -> None:
    """Returns weave tiles GeoDataFrame and tile GeoDataFrame in a
    dictionary based on parameters already supplied to constructor.
    """
    strands_1, strands_2, strands_3 = \
      tiling_utils.get_strand_ids(self.strands)

    loom = self._get_triaxial_weave_matrices(
      strands_1 = strands_1, strands_2 = strands_2, strands_3 = strands_3)

    self._make_shapes_from_coded_weave_matrix(
      loom, strand_labels = [strands_1, strands_2, strands_3])


  def _minimise_triaxial_weave_unit(self) -> None:
    n = self.tiles.shape[0]
    tiles_ids = [x for x in zip(range(n), self.tiles.tile_id, 
                                list(self.tiles.geometry))]
    v_ij = defaultdict(list)
    for i, id_i, tile_i in tiles_ids:
      for j, id_j, tile_j in tiles_ids:
        if i != j and id_i == id_j:
          v = (tile_j.centroid.x - tile_i.centroid.x, 
              tile_j.centroid.y - tile_i.centroid.y)
          if tiling_utils.geometry_matches(
            tile_j, affine.translate(tile_i, v[0], v[1])):
              v_ij[v].append((i, j))
    v_ij = dict(sorted(v_ij.items(), key = lambda x: len(x[1]), reverse = True))
    v_ij = [k for k, v in v_ij.items()][:12]
    combos = itertools.combinations(v_ij, 3)
    result = [c for c in combos 
              if not (tiling_utils.are_parallel(c[0], c[1]) or
                      tiling_utils.are_parallel(c[1], c[2]) or
                      tiling_utils.are_parallel(c[2], c[0]))] # no two parallel
            #   if np.isclose(sum([x for x, y in c]), 0)
            #  and np.isclose(sum([y for x, y in c]), 0)] # resultant (0, 0)
    result = sorted(result, 
                    key = lambda vecs: sum([dx**2 + dy**2 for dx, dy in vecs]))
    poly = tiling_utils.get_prototile_from_vectors(result[0][:3])
    self.prototile.geometry = [poly]
    self.setup_vectors()
    new_tiles = [(id, t.intersection(poly)) for id, t in 
                 zip(self.tiles.tile_id, self.tiles.geometry) 
                 if t.intersection(poly).area > 0]
    self.tiles = gpd.GeoDataFrame({
      'geometry': gpd.GeoSeries([x[1] for x in new_tiles]),
      'tile_id': [x[0] for x in new_tiles]})
    # self.regularise_tiles()
    return None


  def _make_shapes_from_coded_weave_matrix(
    self, loom:Loom, strand_labels:list[list[str]] = [["a"], ["b"], ["c"]]
    ) -> None:
    """Returns weave tiles and prototile GeoDataFrames in a dictionary

    Builds the geometries associated with a given weave supplied as
    'loom' containing the coordinates in an appropriate grid (Cartesian or
    triangular) and the orderings of the strands at each coordinate location

    Args:
      loom (Loom): matrix or stack of matrices representing the weave
      pattern.
      strand_labels (list[list[str]], optional): list of lists of labels
      for strands in each direction. Defaults to [["a"], ["b"], ["c"]]
    """
    grid = WeaveGrid(loom.n_axes, loom.orientations, self.spacing)
    # expand the list of strand labels if needed in each direction
    # labels = []
    labels = [thread * int(np.ceil(dim // len(thread)))
              for dim, thread in zip(loom.dimensions, strand_labels)]
    weave_polys = [] 
    strand_ids = [] 
    cells = []
    for k, strand_order in zip(loom.indices, loom.orderings):
      ids = [thread[coord] for coord, thread in zip(k, labels)]
      cells.append(grid.get_grid_cell_at(k))
      if strand_order is None:
        continue  # No strands present
      if strand_order == "NA":
        continue  # Inconsistency in layer order
      n_slices = [len(id) for id in ids]
      next_polys = grid.get_visible_cell_strands(
        width = self.aspect, coords = k,
        strand_order = strand_order, n_slices = n_slices)
      weave_polys.extend(next_polys)
      next_labels = [list(ids[i]) for i in strand_order]  # list of lists
      next_labels = list(itertools.chain(*next_labels))   # flatten
      strand_ids.extend(next_labels)
    # sometimes empty polygons make it to here, so
    # filter those out along with the associated IDs
    real_polys = [not p.is_empty for p in weave_polys]
    weave_polys = [p for p, b in zip(weave_polys, real_polys) if b]
    strand_ids = [id for id, b in zip(strand_ids, real_polys) if b]
    # note that the tile is important for the biaxial case, which may
    # also a little hard to understand why the behaviour is so differe in 
    # the biaxial and triaxial cases; however below seems to work...
    if loom.n_axes == 3:
      tile = grid.get_tile_from_cells(cells)
      shift = (0, 0)
    else:
      tile = tiling_utils.safe_union(gpd.GeoSeries(cells), as_polygon = True)
      shift = (-tile.centroid.x, -tile.centroid.y)
      tile = grid.get_tile_from_cells(tile)
    self.tiles = self._get_weave_tiles_gdf(weave_polys, strand_ids, shift)
    self.prototile = gpd.GeoDataFrame(
      geometry = gpd.GeoSeries([tile]), crs = self.crs)
    # self.setup_vectors()
    return None


  def _get_weave_tiles_gdf(
      self, polys:list[geom.Polygon], strand_ids:list[str], 
      offset:tuple[float]) -> gpd.GeoDataFrame:
    """Makes a GeoDataFrame from weave tile polygons, labels, etc.

    Args:
      polys (list[Polygon | MultiPolygon]): list of weave tile
        polygons.
      strand_ids (list[str]): list of strand labels.
      offset (tuple[float]): offset to centre the weave tiles on the
        tile.

    Returns:
      geopandas.GeoDataFrame: GeoDataFrame clipped to the tile, with
        margin applied.
    """
    weave = gpd.GeoDataFrame(
      data = {"tile_id": strand_ids},
      geometry = gpd.GeoSeries([affine.translate(p, offset[0], offset[1])
                                for p in polys]))
    weave = weave[weave.tile_id != "-"]
    weave.geometry = gpd.GeoSeries(
      [tiling_utils.get_clean_polygon(p) for p in weave.geometry])
    
    # some buffering is required if aspect is 1 to safely dissolve and
    # explode weave unit tiles that meet at corners
    if self.aspect == 1:
      # grow for dissolve
      weave.geometry = weave.geometry.buffer(
        # self.spacing * tiling_utils.RESOLUTION, 
        tiling_utils.RESOLUTION, 
        join_style = 2, cap_style = 3)
      weave = weave.dissolve(by = "tile_id", as_index = False)
      # shrink by more to explode into separate polygons
      weave.geometry = weave.geometry.buffer(
        # -2 * self.spacing * tiling_utils.RESOLUTION, 
        -2 * tiling_utils.RESOLUTION, 
        join_style = 2, cap_style = 3)
      weave = weave.explode(ignore_index = True)
      weave.geometry = weave.geometry.buffer(
        # self.spacing * tiling_utils.RESOLUTION, 
        tiling_utils.RESOLUTION, 
        join_style = 2, cap_style = 3)
    else: # aspect < 1 is fine without buffering
      weave = weave.dissolve(by = "tile_id", as_index = False)
      weave = weave.explode(ignore_index = True)

    weave.geometry = gpd.GeoSeries(
      [tiling_utils.get_clean_polygon(p) for p in weave.geometry])
    return weave.set_crs(self.crs)


  def _get_axis_from_label(self, label:str = "a", strands:str = None):
    """Determines the axis of a tile_id from the strands spec string.

    Args:
      label (str, optional): the tile_id. Defaults to "a".
      strands (str, optional): the strand spec. Defaults to the WeaveUnit
      strands attribute.

    Returns:
      _type_: the axis in which the supplied tile is found.
    """
    if strands == None:
      strands = self.strands
    index = strands.index(label)
    return strands[:index].count("|")


  def _get_legend_tiles(self) -> gpd.GeoDataFrame:
    """Returns tiles suitable for use in a legend representation.

    One tile for each tile_id value will be chosen, close to the
    centre of the prototile extent, and not among the smallest tiles present
    (for example not a short length of strand mostly hidden by other
    strands)

    Returns:
      gpd.GeoDataFrame: the chosen tiles.
    """
    angles = ((0, 240, 120)
              if self.weave_type in ("hex", "cube")
              else (90, 0))
    tile_ids = pd.Series.unique(self.tiles.tile_id)
    groups = self.tiles.groupby("tile_id")
    tiles, rotations = [], []
    for id in tile_ids:
      candidates = groups.get_group(id)
      axis = self._get_axis_from_label(id, self.strands)
      tiles.append(self._get_most_central_large_tile(
      candidates, tiles))
      rotations.append(-angles[axis] + self.rotation)
    return gpd.GeoDataFrame(
      data = {"tile_id": tile_ids, "rotation": rotations},
      crs = self.crs,
      geometry = gpd.GeoSeries(tiles))


  def _get_most_central_large_tile(self, tiles:gpd.GeoDataFrame,
          other_tiles:list[geom.Polygon]) -> geom.Polygon:
    """Gets a large tile close to the centre of the WeaveUnit.

    Args:
      tiles (gpd.GeoDataFrame): the set of tiles to choose from.

    Returns:
      geom.Polygon: the chosen, large central tile.
    """
    areas = [g.area for g in tiles.geometry]
    min_area, max_area = min(areas), max(areas)
    if min_area / max_area > 0.5:
      geoms = list(tiles.geometry)
    else:
      mean_log_a = np.mean(np.log(areas))
      geoms = [g for g, a in zip(tiles.geometry, areas)
              if np.log(a) > mean_log_a]
    if len(other_tiles) == 0 or self.weave_type in ("cube", "hex"):
      d = [g.centroid.distance(geom.Point(0, 0)) for g in geoms]
    else:
      c = geom.MultiPolygon(other_tiles).centroid
      d = [geom.MultiPolygon([g] + other_tiles).centroid.distance(c)
          for g in geoms]
    return geoms[d.index(min(d))]


  def _get_legend_key_shapes(self, polygon:geom.Polygon,
                             counts:Iterable = [1] * 25, angle:float = 0,
                             radial:bool = False) -> list[geom.Polygon]:
    """Returns a list of polygons obtained by slicing the supplied polygon
    across its length into n slices. Orientation of the polygon is
    indicated by the angle.

    The returned list of polygons can be used to form a colour ramp in a
    legend.

    Args:
      polygon (geom.Polygon): the weave strand polygon to slice.
      counts (Iterable, optional): an iterable list of the numbers of
        slices in each category. Defaults to [1] * 25.
      angle (float, optional): orientation of the polygon. Defaults to 0.
      categorical (bool, optional): ignored by WeaveUnit.

    Returns:
      list[geom.Polygon]: a list of polygons.
    """
    c = polygon.centroid
    g = affine.rotate(polygon, -angle, origin = c)
    width, height, left, bottom = \
      tiling_utils.get_width_height_left_bottom(gpd.GeoSeries([g]))
    total = sum(counts)
    cuts = list(np.cumsum(counts))
    cuts = [0] + [c / total for c in cuts]
    cuts = [left + c * width for c in cuts]
    # add margin to avoid weird effects intersecting almost parallel lines.
    cuts[0] = cuts[0] - 1
    cuts[-1] = cuts[-1] + 1
    bottom = bottom - 1
    top = bottom + height + 1
    slices = []
    for l, r in zip(cuts[:-1], cuts[1:]):
      slice = geom.Polygon([(l, bottom), (r, bottom), (r, top), (l, top)])
      slices.append(slice.intersection(g))
    return [affine.rotate(s, angle, origin = c) for s in slices]


