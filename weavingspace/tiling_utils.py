#!/usr/bin/env python
# coding: utf-8

# don't understand it but this allows import of TileUnit for type hinting only
from __future__ import annotations
from copy import deepcopy
from typing import TYPE_CHECKING
if TYPE_CHECKING:
  from weavingspace import TileUnit
from typing import Iterable, Union

import re
import string
from collections import namedtuple
from collections import defaultdict
from itertools import combinations
from itertools import chain

import numpy as np
from scipy import interpolate

import geopandas as gpd
import pandas as pd
import shapely.algorithms.polylabel as polylabel
import shapely.geometry as geom
import shapely.affinity as affine
import shapely.ops


"""Utility functions for use by various classes in the `weavingspace`
package. Most are geometric convenience functions for commonly applied
operations.
"""

PRECISION = 6
RESOLUTION = 1e-6

def _parse_strand_label(s:str) -> list[str]:
  """Breaks a strand label specification in to a list of labels

  Args:
    s (str): see get_strand_ids() for details
  Returns:
    list[str]: list of strand labels.
  """
  clean_s = re.sub("[(]+", "(", re.sub("[)]+", ")", s))
  is_combo = False
  output = []
  current = ""
  for c in clean_s:
    if is_combo:
      if c == ")":
        output.append(current)
        current = ""
        is_combo = False
      else:
        current = current + c
    else:
      if c == "(":
        is_combo = True
      else:
        output.append(c)
  return output


def get_strand_ids(strands_spec: str) -> tuple[list[str]]:
  """Conversts a strands specification string to a list of lists of strand
  labels.

  Args:
    strands_spec (str): string format "a|bc|(de)f" | separates strands in
      each direction and () designates combining labels into a single
      strand that will be sliced lengthwise. Example output:

        "a|bc|(de)f" -> (["a"], ["b", "c"], ["de", "f"])

    Superflous parentheses are removed, but no other error-checks are
      applied.

  Returns:
    tuple[str]: tuple of lists of labels for each set of strands.
  """
  strand_ids = [_parse_strand_label(s) for s in strands_spec.split("|")]
  strand_ids = strand_ids if len(strand_ids) == 3 else strand_ids + [[""]]
  return tuple(strand_ids)


def centre_offset(shape: geom.Polygon,
                  target:tuple[float] = (0, 0)) -> tuple[float]:
  """Returns vector required to move centroid of polygon to target.

  Args:
    shape (Polygon): polygon to move.
    target (tuple[float], optional): target to move to.
      Defaults to (0, 0).

  Returns:
    tuple[float]: tuple of x, y movement required.
  """
  shape_c = shape.centroid
  return (target[0] - shape_c.x, target[1] - shape_c.y)


def get_regular_polygon(spacing, n:int) -> geom.Polygon:
  """Returns regular polygon with n sides centered on (0, 0) with a horizontal base, and height given by spacing.

  Args:
    spacing (_type_): required height.
    n (int): number of sides.

  Returns:
    geom.Polygon: required geom.Polygon.
  """
  R = spacing / np.cos(np.radians(180 / n)) / 2
  a0 = -90 + 180 / n
  a_diff = 360 / n
  angles = [a0 + a * a_diff for a in range(n)]
  corners = [(R * np.cos(np.radians(a)),
              R * np.sin(np.radians(a))) for a in angles]
  return gridify(geom.Polygon(corners))


def offset_polygon_corners(polygon:geom.Polygon, 
                          offset:int) -> geom.Polygon:
  """Returns this polygon but with its first corner offset from its
  original position in the coordinate sequence. The returned polygon will
  be identical but stored differently internally.

  Args:
    polygon (geom.Polygon): the polygon to reorder.
    offset (int): the number of corner positions by which to shift the
      sequence.

  Returns:
      geom.Polygon: the reordered polygon.
  """
  corners = get_corners(polygon, repeat_first = False)
  return geom.Polygon([c for c in corners[offset:] + corners[:offset]]) 


def rotate_preserving_order(polygon:geom.Polygon, angle:float,
                            centre:geom.Point) -> geom.Polygon:
  """Returns the supplied polygon rotated with the order of its corner points
  preserved (not guaranteed by shapely.affinity.rotate).

  Args:
      polygon (geom.Polygon): polygon to rotate.
      angle (float): desired angle of rotation (in degrees).
      centre (geom.Point): the rotation centre (passed on to shapely.affinity.
        rotate).

  Returns:
      geom.Polygon: rotated polygon.
  """
  corners = get_corners(polygon, repeat_first = False)
  return geom.Polygon([affine.rotate(c, angle, centre) for c in corners])


def geometry_matches(geom1:geom.Polygon, geom2:geom.Polygon):
  a, b = geom1.area, geom2.area
  return bool(
    np.isclose(a, b, rtol = RESOLUTION * 100, atol = RESOLUTION * 100) and 
    np.isclose(a, geom1.intersection(geom2).area, rtol = RESOLUTION * 100, 
               atol=RESOLUTION * 100))


def are_parallel(v1:tuple[float], v2:tuple[float]) -> bool:
  a1 = np.arctan2(v1[1], v1[0])
  a2 = np.arctan2(v2[1], v2[0])
  return bool(np.isclose(a1, a2) or np.isclose(np.abs(a1 - a2), np.pi))


def is_convex(shape:geom.Polygon) -> bool:
  """Tests for shape convexity. There are better ways to do this, like
  e.g. using something like the get_interior_angles() function, but simply
  checking if the area is close to that of its convex hull works too!

  Args:
    shape (geom.Polygon): polygon to check

  Returns:
    bool: True if the polygon is convex, False otherwise.
  """
  return np.isclose(shape.area, shape.convex_hull.area, rtol = RESOLUTION)


def get_corners(shape:geom.Polygon,
                repeat_first:bool = True) -> list[geom.Point]:
  """Returns a list of geom.Points around the boundary of a polygon, optionally
  repeating the first. Does no simplification (e.g. if a line segment has a 
  'corner' along its length, it is NOT removed; see get_clean_polygon for 
  that). Points have precision set to the package default tiling_utils.
  RESOLUTION.

  Args:
    shape (geom.Polygon): polygon whose corners are required.
    repeat_first (bool, optional): if True the first corner is repeated in the 
      returned list, if False it is omitted. Defaults to True.

  Returns:
    list[geom.Point]: list of geom.Point vertices of the polygon.
  """
  corners = [gridify(geom.Point(pt)) for pt in shape.exterior.coords]
  if repeat_first:
    return corners
  else:
    return corners[:-1]


def get_sides(shape:geom.Polygon) -> list[geom.LineString]:
  """Returns polygon sides as a list of geom.LineStrings, with resolution set
  to the package default tiling_utils.RESOLUTION. No simplification for 
  successive colinear sides is carried out.

  Args:
    shape (geom.Polygon): polygon whose edges are required.

  Returns:
    list[geom.LineString]: list of geom.LineString sides of the polygon.
  """
  corners = get_corners(shape)
  return [gridify(geom.LineString([p1, p2]))
          for p1, p2 in zip(corners[:-1], corners[1:])]


def get_side_lengths(shape:geom.Polygon) -> list[float]:
  """Returns list of lengths of polygon sides. No simplification for corners
  along sides is carried out.

  Args:
    shape (geom.Polygon): polygon whose edge lengths are required.

  Returns:
    list[float]: list of side lengths.
  """
  corners = get_corners(shape)
  return [p1.distance(p2) for p1, p2 in zip(corners[:-1], corners[1:])]


def get_side_bearings(shape:geom.Polygon) -> tuple[tuple[float]]:
  """Returns a list of angles (in degrees) between the sides of a polygon and
  the positive x-axis, when proceeding from the first point in each side to its
  end point. This should usually be CW around the polygon.

  Args:
    shape (geom.Polygon): polygon whose side bearings are required.

  Returns:
    tuple[tuple[float]]: tuple of bearings of each edge.
  """
  return tuple([np.degrees(np.arctan2(
    e.coords[1][1] - e.coords[0][1], e.coords[1][0] - e.coords[0][0])) 
    for e in get_sides(shape)])


def get_interior_angles(shape:geom.Polygon) -> list[float]:
  """Returns angles (in degrees) between successive edges of a polygon. No 
  polygon simplification is carried out so some angles may be 180 (i.e. a 
  'corner' along a side, such that successive sides are colinear). 

  Args:
    shape (geom.Polygon): polygon whose angles are required.

  Returns:
    list[float]: list of angles.
  """
  corners = get_corners(shape, repeat_first = False)
  wrap_corners = corners[-1:] + corners + corners[:1]
  triples = zip(wrap_corners[:-2], wrap_corners[1:-1], wrap_corners[2:])
  return [get_inner_angle(p1, p2, p3) for p1, p2, p3 in triples]


def get_clean_polygon(
    shape:Union[geom.MultiPolygon,geom.Polygon]) -> geom.Polygon:
  """Returns polygon with any successive corners that are co-linear along a 
  side or very close to one another removed. 
  
  Particularly useful for tidying polygons weave tilings that have been 
  assembled from multiple 'cells' in the weave grid.

  Args:
    shape (Union[geom.MultiPolygon,geom.Polygon]): polygon to clean.

  Returns:
    geom.Polygon: cleaned polygon.
  """
  if isinstance(shape, geom.MultiPolygon):
    return geom.MultiPolygon([get_clean_polygon(p) for p in shape.geoms])
  else:
    corners = get_corners(shape, repeat_first = False)
    # first remove any 'near neighbour' corners
    distances = get_side_lengths(shape)
    to_remove = [np.isclose(d, 0, rtol = RESOLUTION, atol = 10 * RESOLUTION) 
                 for d in distances]
    corners = [c for c, r in zip(corners, to_remove) if not r]
    # next remove any that are colinear
    p = geom.Polygon(corners)
    corners = get_corners(p, repeat_first = False)
    angles = get_interior_angles(p)
    to_remove = [np.isclose(a, 0, rtol = RESOLUTION, atol = RESOLUTION) or
                 np.isclose(a, 180, rtol = RESOLUTION, atol = RESOLUTION)
                 for a in angles]
    corners = [c for c, r in zip(corners, to_remove) if not r]
    return gridify(geom.Polygon(corners))


def get_inner_angle(p1:geom.Point, p2:geom.Point, p3:geom.Point) -> float:
  r"""Returns the angle (in degrees) between line p1-p2 and p2-p3, i.e., the 
  angle A below
  
            p2
           / \ 
          / A \ 
        p1     p3
  
  Args:
    p1 (geom.Point): first point.
    p2 (geom.Point): second 'corner' point.
    p3 (geom.Point): third point.

  Returns:
      float: angle in degrees.
  """
  return np.degrees(np.arctan2(p3.y - p2.y, p3.x - p2.x) -
                    np.arctan2(p1.y - p2.y, p1.x - p2.x)) % 360


def get_outer_angle(p1:geom.Point, p2:geom.Point, p3:geom.Point) -> float:
  r"""Returns outer angle (in degrees) between lines p1-p2 and p2-p3, i.e., the 
  angle A below
                
              /
             /
            p2 A
           / \ 
          /   \ 
        p1     p3
  
  Args:
    p1 (geom.Point): first point.
    p2 (geom.Point): second 'corner' point.
    p3 (geom.Point): third point.

  Returns:
    float: angle in degrees.
  """
  return 180 - get_inner_angle(p1, p2, p3)


def is_regular_polygon(shape:geom.Polygon) -> bool:
  """Tests if supplied polygon is regular (i.e. equal sides and angles).

  Args:
      shape (geom.Polygon): polygon to test.

  Returns:
      bool: True if polygon is regular, False if not.
  """
  side_lengths = get_side_lengths(shape)
  angles = get_interior_angles(shape)
  return all(np.isclose(side_lengths, side_lengths[0])) \
     and all(np.isclose(angles, angles[0]))


def is_tangential(shape:geom.Polygon) -> bool:
  """Determines if the supplied polygon is tangential i.e., it can have
  circle inscribed tangential to all its sides. 
  
  Note that this will fail for polygons with successive colinear sides,
  meaning that polygons should be fully simplified...
  """
  if is_regular_polygon(shape):
    return True
  side_lengths = get_side_lengths(shape)
  n = len(side_lengths)
  if n % 2 == 1:
    if n == 3: #triangles are easy!
      return True
    # odd number of sides there is a nice solvable system  of equations
    # see https://math.stackexchange.com/questions/4065370/tangential-polygons-conditions-on-edge-lengths
    #
    # TODO: this is not reliable for reasons I don't fully understand
    # there is a corresponding crude catch exception in incentre... but
    # this needs a more complete investigation
    mat = np.identity(n) + np.roll(np.identity(n), 1, axis = 1)
    fractions = np.linalg.inv(mat) @ side_lengths / side_lengths
    if not(np.isclose(np.mean(fractions), 
                      0.5, rtol= RESOLUTION, atol = RESOLUTION)):
      return False
    ones = [np.isclose(f, 1, rtol = RESOLUTION, atol = RESOLUTION) 
            for f in fractions]
    zeros = [np.isclose(f, 0, rtol = RESOLUTION, atol = RESOLUTION) 
             for f in fractions]
    negatives = [f < 0 for f in fractions]
    return not (any(ones) or any(zeros) or any(negatives))
  elif n == 4:
    # for quadrilaterals odd and even side lengths are equal
    return np.isclose(sum(side_lengths[0::2]), 
                      sum(side_lengths[1::2]), rtol = RESOLUTION)
  else:
    # other even numbers of sides... hard to avoid brute force
    bisectors = [get_angle_bisector(shape, i) for i in range(n)]
    intersects = [b1.intersection(b2) for b1, b2 in
                  zip(bisectors, bisectors[1:] + bisectors[:1])
                  if b1.intersects(b2)]
    distances = [i1.distance(i2) for i1, i2 in
                 zip(intersects, intersects[1:] + intersects[:1])]
    return all([d <= 2 * RESOLUTION for d in distances])


def incentre(shape:geom.Polygon) -> geom.Point:
  """A different polygon centre, which produces better results for some
  dual tilings where tiles are not regular polygons... see
  https://en.wikipedia.org/wiki/Incenter

  This method relies on the polygon being tangential, i.e. there is an
  inscribed circle to which all sides of the polygon are tangent. It will
  work on all the polygons encountered in the Laves tilings, but is not
  guaranteed to work on all polygons.

  Given that the polygon is tangential, the radius of the inscribed circle is
  the [apothem of the polygon](https://en.wikipedia.org/wiki/Apothem) given
  by 2 x Area / Perimeter. We apply a parallel offset of this size to two
  sides of the polygon and find their intersection to determine the centre of
  the circle, givng the incentre of the polygon.

    Args:
    shape (geom.Polygon): the polygon.

  Returns:
    geom.Point: the incentre of the polygon.
  """
  if is_regular_polygon(shape):
    return shape.centroid
  return polylabel.polylabel(shape)
  # shape = ensure_cw(shape)
  # if is_convex(shape):
  #   if is_tangential(shape):  # find the incentre
  #     corners = get_corners(shape, repeat_first = False)
  #     r = get_apothem_length(shape)
  #     e1 = geom.LineString(corners[:2]).parallel_offset(r, side = "right")
  #     e2 = geom.LineString(corners[1:3]).parallel_offset(r, side = "right")
  #     c = e1.intersection(e2)
  #     # is_tangential is unreliable, and sometimes we get to here and do not
  #     # get a Point, but a LineString, or even no intersection, so...
  #     return c if isinstance(c, geom.Point) else shape.centroid
  # return shape.centroid


def in_circle(shape:geom.Polygon) -> geom.Polygon:
  c = incentre(shape)
  r = min([c.distance(e) for e in get_sides(shape)])
  return c.buffer(r)


def get_apothem_length(shape:geom.Polygon) -> float:
  return 2 * shape.area / geom.LineString(shape.exterior.coords).length


def ensure_cw(shape:geom.Polygon) -> geom.Polygon:
  """Returns the polygon with its outer boundary vertices in clockwise order.

  It is important to note that shapely.set_precision() imposes clockwise order
  on polygons, and since it is used widely throughout theses modules, it makes
  sense to impose this order.

    Args:
    shape (geom.Polygon): the polygon.

  Returns:
    geom.Polygon: the polygon in clockwise order.
  """
  if geom.LinearRing(shape.exterior.coords).is_ccw:
    return shape.reverse()
  else:
    return shape


def get_angle_bisector(shape:geom.Polygon, v = 0) -> geom.LineString:
  """Returns a line which is the angle bisector of the specified corner of the
  supplied polygon.

  Args:
      shape (geom.Polygon): the polygon
      v (int, optional): index of the corner whose bisector is required. 
      Defaults to 0.

  Returns:
      geom.LineString: line which bisects the specified corner.
  """
  pts = get_corners(shape, repeat_first = False)
  n = len(pts)
  p0 = pts[v % n]
  p1 = pts[(v + 1) % n]
  p2 = pts[(v - 1) % n]
  d01 = p0.distance(p1)
  d02 = p0.distance(p2)
  if d01 < d02:
    p2 = geom.LineString([p0, p2]).interpolate(d01)
  else:
    p1 = geom.LineString([p0, p1]).interpolate(d02)
  m = geom.Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2)
  # we have to scale this to much larger in case of an angle near 180
  return affine.scale(geom.LineString([p0, m]), 10, 10, origin = p0)


def get_side_bisector(shape:geom.Polygon, i:int = 0) -> geom.LineString:
  return affine.scale(affine.rotate(get_sides(shape)[i], 90), 100, 100)


def _get_interior_vertices(tiles:gpd.GeoDataFrame) -> gpd.GeoSeries:
  """Returns points not on the outer boundary of the supplied set
  of tiles.

  Args:
    shapes (gpd.GeoDataFrame): a set of polygons.

  Returns:
    gpd.GeoSeries: interior vertices of the set of polygons.
  """
  tiles = gridify(tiles.geometry)
  uu = safe_union(tiles, as_polygon = True)
  interior_pts = set()
  for tile in tiles:
    for corner in tile.exterior.coords:
      if uu.contains(geom.Point(corner).buffer(1e-3, cap_style = 3)):
        interior_pts.add(corner)
  return gridify(gpd.GeoSeries([geom.Point(p) for p in interior_pts]))


def gridify(
      gs:Union[gpd.GeoSeries, gpd.GeoDataFrame, geom.Polygon]
    ) -> Union[gpd.GeoSeries, gpd.GeoDataFrame, geom.Polygon]:
  """Returns the supplied GeoSeries rounded to the specified precision.

  Args:
    gs (gpd.GeoSeries): geometries to gridify.
    precision (int, optional): digits of precision. Defaults to 6.

  Returns:
    gpd.GeoSeries: the rounded geometries.
  """
  # return gpd.GeoSeries(
  #   list(gs.apply(
  #     wkt.dumps, rounding_precision = precision).apply(wkt.loads)))
  if isinstance(gs, (geom.Polygon, geom.Point, geom.MultiPoint,
                     geom.MultiPolygon, geom.LineString)) :
    return shapely.set_precision(gs, grid_size = RESOLUTION)
  if isinstance(gs, gpd.GeoSeries):
    return gpd.GeoSeries(
      [shapely.set_precision(g, grid_size = RESOLUTION) for g in list(gs)],
      crs = gs.crs)
  if isinstance(gs, gpd.GeoDataFrame):
    gs.geometry = gridify(gs.geometry)
    return gs


def get_dual_tile_unit(unit: TileUnit) -> gpd.GeoDataFrame:
  """Converts supplied TileUnit to a candidate GeoDataFrame of its dual
  TileUnit.

  NOTE: this is complicated and not remotely guaranteed to work!
  a particular issue is that where to place the vertices of the faces
  of the dual with respect to the tiles in the original is ill-defined.
  This is because the dual process is topologically not metrically defined,
  so that exact vertex locations are ambiguous. Tiling duality is defined in
  Section 4.2 of Grunbaum B, Shephard G C, 1987 _Tilings and Patterns_ (W. H.
  Freeman and Company, New York)

  NOTE: In general, this method will work only if all supplied tiles are
  regular polygons. A known exception is if the only non-regular polygons are
  triangles.

  NOTE: 'clean' polygons are required. If supplied polygons have messy
  vertices with multiple points where there is only one proper point, bad
  things are likely to happen! Consider using `clean_polygon()` on the
  tile geometries.

  Because of the above limitations, we only return a GeoDataFrame
  for inspection. However some `weavingspace.tile_unit.TileUnit` setup
  methods in `weavingspace.tiling_geometries` use this method, where we are
  confident the returned dual is valid.

  Args:
    unit (TileUnit): the tiling for which the dual is required.

  Returns:
    gpd.GeoDataFrame: GeoDataFrame that could be the tiles attribute for
      a TileUnit of the dual tiling.
  """
  # get a local patch of this Tiling
  local_patch = unit.get_local_patch(r = 1, include_0 = True)
  # Find the interior points of these tiles - these will be guaranteed
  # to have a sequence of surrounding tiles incident on them
  interior_pts = _get_interior_vertices(local_patch)
  # Compile a list of the polygons incident on the interior points
  cycles = []
  for pt in interior_pts:
    cycles.append(
      set([poly_id for poly_id, p in enumerate(local_patch.geometry)
           if pt.distance(p) < RESOLUTION * 2]))
  # convert the polygon ID sequences to (centroid.x, centroid.y, ID) tuples
  dual_faces = []
  for cycle in cycles:
    ids, pts = [], []
    for poly_id in cycle:
      ids.append(local_patch.tile_id[poly_id])
      poly = local_patch.geometry[poly_id]
      pts.append(incentre(poly))
    # sort them into CW order so they are well formed
    sorted_coords = sort_cw([(pt.x, pt.y, id)
                              for pt, id in zip(pts, ids)])
    dual_faces.append(
      (geom.Polygon([(pt_id[0], pt_id[1]) for pt_id in sorted_coords]),
       "".join([pt_id[2] for pt_id in sorted_coords])))
  # ensure the resulting face centroids are inside the original tile
  # displaced a little to avoid uncertainties at corners/edges
  # TODO: Check  the logic of this - it seems like dumb luck that it works...
  dual_faces = [(f, id) for f, id in dual_faces
                if affine.translate(
                  unit.prototile.loc[0, "geometry"], RESOLUTION * 10, 
                  RESOLUTION * 10).contains(f.centroid)]
  gdf = gpd.GeoDataFrame(
    data = {"tile_id": [f[1] for f in dual_faces]}, crs = unit.crs,
    geometry = gridify(gpd.GeoSeries([f[0] for f in dual_faces])))
  # ensure no duplicates
  gdf = gdf.dissolve(by = "tile_id", as_index = False).explode(
    index_parts = False, ignore_index = True)

  gdf.tile_id = relabel(gdf.tile_id)
  return gdf


def relabel(data:Iterable) -> list:
  """Returns supplied data reassigned with unique values from
  string.ascii_letters.

  Args:
    data (Iterable): the data to relabel

  Returns:
    list: the reassigned data
  """
  new_data = {}
  d_count = 0
  for d in data:
    if not d in new_data:
      new_data[d] = string.ascii_letters[d_count]
      d_count = d_count + 1
  return [new_data[d] for d in data]


def sort_cw(pts_ids:list[tuple[float, float, str]]):
  """Sorts supplied tuple of x, y, ID into clockwise order relative to their
  mean centre.

  Args:
    pts_ids (list[tuple[float, float, str]]): A tuple of a pair of
      floats and a string.

  Returns:
    list: a list in the same format as supplied sorted into
      clockwise order of the point locations.
  """
  x = [p[0] for p in pts_ids]
  y = [p[1] for p in pts_ids]
  cx, cy = np.mean(x), np.mean(y)
  dx = [_ - cx for _ in x]
  dy = [_ - cy for _ in y]
  angles = [np.arctan2(dy, dx) for dx, dy in zip(dx, dy)]
  d = dict(zip(angles, pts_ids))
  return [pt_id for angle, pt_id in reversed(sorted(d.items()))]


def order_of_pts_cw_around_centre(pts:list[geom.Point], centre:geom.Point):
  """Returns the order of the supplied points clockwise relative to supplied 
  centre point, i.e. a list of the indices in clockwise order.

  Args:
      pts (list[geom.Point]): list of points to order.
      centre (geom.Point): centre relative to which CW order is determined.

  Returns:
      _type_: list of reordered points.
  """
  dx = [p.x - centre.x for p in pts]
  dy = [p.y - centre.y for p in pts]
  angles = [np.arctan2(dy, dx) for dx, dy in zip(dx, dy)]
  d = dict(zip(angles, range(len(pts))))
  return [i for angle, i in reversed(sorted(d.items()))]


def write_map_to_layers(gdf:gpd.GeoDataFrame, fname:str = "output.gpkg",
                        tile_var:str = "tile_id") -> None:
  """Writes supplied GeoDataFrame to a GPKG file with layers based on
  the tile_var attribute.

  Args:
    gdf (gpd.GeoDataFrame): the GeoDataFrame.
    fname (str, optional): filename to write.
    tile_var (str, optional): the attribute to use to separate
      output file into layers. Defaults to "tile_id".
  """
  grouped = gdf.groupby(tile_var, as_index = False)
  for e in pd.Series.unique(gdf[tile_var]):
    grouped.get_group(e).to_file(fname, layer = e, driver = "GPKG")


def get_collapse_distance(geometry:geom.Polygon) -> float:
  """Returns the distance under which the supplied polygon will shrink
  to nothing if negatively buffered by that distance.

  Performs a binary search between an upper bound based on the radius of
  the circle of equal area to the polygon, and 0.

  Args:
    geometry (geom.Polygon): the polygon.

  Returns:
    float: its collapse distance.
  """
  if is_convex(geometry):
    return get_apothem_length(geometry)
  radius = np.sqrt(geometry.area / np.pi)
  lower = 0
  upper = radius
  delta = 1e12
  stop_delta = radius / 1000
  while delta > stop_delta:
    new_r = (lower + upper) / 2
    if geometry.buffer(-new_r).area > 0:
      lower = new_r
      delta = upper - new_r
    else:
      upper = new_r
      delta = new_r - lower
  return new_r


def get_largest_polygon(polygons:gpd.GeoSeries) -> gpd.GeoSeries:
  """Returns the largest polygon in a GeoSeries as a GeoSeries of one polygon.

  Args:
    polygons (gpd.GeoSeries): the set of polygons to pick from.

  Returns:
    gpd.GeoSeries: the largest polygon.
  """
  actual_polygons = [p for p in polygons.geometry
                     if isinstance(p, (geom.Polygon, geom.MultiPolygon))]
  areas = [p.area for p in actual_polygons]
  max_area = max(areas)
  return gpd.GeoSeries([actual_polygons[areas.index(max_area)]])


def touch_along_an_edge(p1:geom.Polygon, p2:geom.Polygon) -> bool:
  """Tests if two polygons touch along an edge.

  Checks that the intersection area of the two polygons buffered by
  a small amount is large enough to indicate that they neighbour at more
  than a corner.

  Args:
    p1 (geom.Polygon): First polygon
    p2 (geom.Polygon): Second polygon

  Returns:
    bool: True if they neighbour along an edge
  """
  return p1.buffer(RESOLUTION, join_style = 2, cap_style = 3) \
    .intersection(p2.buffer(RESOLUTION, join_style = 2, cap_style = 3)) \
      .area > 16 * RESOLUTION ** 2


def get_width_height_left_bottom(gs:gpd.GeoSeries) -> tuple[float]:
  """Returns width, height, left and bottom limits of a GeoSeries

  Args:
    gs (geopandas.GeoSeries): GeoSeries for which limits are required.

  Returns:
    tuple: four float values of width, height, left and bottom of gs.
  """
  extent = gs.total_bounds
  return (extent[2] - extent[0], extent[3] - extent[1],
      extent[0], extent[1])


def get_bounding_ellipse(shapes:gpd.GeoSeries, 
                         mag:float = 1.0) -> gpd.GeoSeries:
  """Returns an ellipse containing the supplied shapes.

  The method used is to calculate the size of square that would contain
  the shapes, if they had an aspect ratio 1, then stretch the circle in
  the x, y directions according to the actual aspect ratio of the shapes.

  Args:
    shapes (gpd.GeoSeries): the shapes to be contained.
    mag (float, optional): optionally increase the size of the returned
      ellipse by this scale factor. Defaults to 1.0.

  Returns:
    gpd.GeoSeries: the set of shapes.
  """

  w, h, l, b = get_width_height_left_bottom(shapes)

  c = geom.Point(l + w / 2, b + h / 2)
  r = min(w, h) * np.sqrt(2)
  circle = [c.buffer(r)]
  return gridify(gpd.GeoSeries(circle, crs = shapes.crs).scale(
    w / r * mag / np.sqrt(2), h / r * mag / np.sqrt(2), origin = c))


def get_tiling_edges(tiles:gpd.GeoSeries) -> gpd.GeoSeries:
  """Returns linestring GeoSeries from supplied polygon GeoSeries.

  This is used to allow display of edges of tiles in legend when they are
  masked by an ellipse (if we instead clip polygons then the ellipse edge
  will also show in the result.)

  Args:
    shapes (gpd.GeoSeries): Polygons to convert.

  Returns:
    gpd.GeoSeries: LineStrings from the supplied Polygons.
  """
  tiles = tiles.explode(ignore_index = True)
  edges = [geom.LineString(p.exterior.coords) for p in tiles]
  return gpd.GeoSeries(edges, crs = tiles.crs)


def get_polygon_sector(polygon:geom.Polygon, start:float = 0.0,
                       end:float = 1.0) -> geom.Polygon:
  """Returns a sector of the provided Polygon.

  The returned sector is a section of the polygon boundary between the
  normalized start and end positions, and including the polygon centroid.
  Should (probably) only be applied to convex polygons.

  Args:
    shape (geom.Polygon): the Polygon.
    start (float): normalized start position along the boundary. Defaults to
      0.
    end (float): normalized start position along the boundary. Defaults to
      1.

  Returns:
    geom.Polygon: the requested polygon sector.
  """
  if start == end:
    # must return a null polygon since the calling context
    # expects to get something back... which most likely
    # is needed to align with other data
    return geom.Polygon()
  if start * end < 0:
    # either side of 0/1 so assume required sector includes 0
    e1, e2 = min(start, end), max(start, end)
    arc1 = shapely.ops.substring(geom.LineString(polygon.exterior.coords),
                                 np.mod(e1, 1), 1, normalized = True)
    arc2 = shapely.ops.substring(geom.LineString(polygon.exterior.coords),
                                 0, e2, normalized = True)
    sector = geom.Polygon([polygon.centroid] +
               list(arc1.coords) +
               list(arc2.coords)[1:])
  else:
    arc = shapely.ops.substring(geom.LineString(polygon.exterior.coords),
                                start, end, normalized = True)
    sector = geom.Polygon([polygon.centroid] + list(arc.coords))
  return gridify(sector)


def repair_polygon(
    polygon:Union[geom.Polygon, geom.MultiPolygon, gpd.GeoSeries],
    shrink_then_grow:bool = True) -> Union[geom.Polygon, gpd.GeoSeries]:
  """Convenience function to 'repair' a shapely polyon or GeoSeries by applying
  a negative buffer then the same positive buffer.

  Optionally the buffer may be applied in the opposite order (i.e. grow then
  shrink). This operation may also convert a MultiPolygon that has some 'stray'
  parts to a Polygon.

  This is method is often unofficially recommended (on stackexchange etc.)
  even in the shapely docs, to resolve topology issues and extraneous
  additional vertices appearing when spatial operations are repeatedly
  applied.

  Args:
    p (Union[geom.Polygon, gpd.GeoSeries]): Polygon or GeoSeries to clean.
    res (float, optional): buffer size to use. Defaults to 1e-3.
    shrink_then_grow (bool, optional): if True the negative buffer is
      applied first, otherwise the buffer operations are applied in
      reverse. Defaults to True.

  Returns:
    Union[geom.Polygon, gpd.GeoSeries]: the cleaned Polygon or GeoSeries.
  """
  if type(polygon) is gpd.GeoSeries:
    return gpd.GeoSeries([repair_polygon(p, shrink_then_grow) 
                          for p in polygon]).set_crs(polygon.crs)
  else:
    if polygon is None:
      return None
    if shrink_then_grow:
      return polygon.buffer(
        -RESOLUTION * 10, join_style = 2, cap_style = 3).buffer(
        RESOLUTION * 10, join_style = 2, cap_style = 3)
    else:
      return polygon.buffer(
        RESOLUTION * 10, join_style = 2, cap_style = 3).buffer(
        -RESOLUTION * 10, join_style = 2, cap_style = 3)


def safe_union(gs:gpd.GeoSeries,
               as_polygon:bool = False) -> Union[gpd.GeoSeries, geom.Polygon]:
  """Unions the supplied GeoSeries of Polygons while buffering them to avoid
  gaps and odd internal floating edges. Optionally returns a Polygon or a
  GeoSeries.

  Frequently when unioning polygons that are ostensibly adjacent 'rogue'
  internal boundaries remain in the result. We can avoid this by buffering the
  polygons before unioning them, then reversing the buffer on the unioned
  shape.

  Args:
    gs (gpd.GeoSeries): the Polygons to union.
    res (float, optional): size of the buffer to use. Defaults to 1e-3.
    as_polygon (bool, optional): if True returns a Polygon, otherwise
      returns a one Polygon GeoSeries. Defaults to False.

  Returns:
    Union[gpd.GeoSeries, geom.Polygon]: the resulting union of supplied
      polygons.
  """
  union = gpd.GeoSeries(
    [p.buffer(RESOLUTION * 10, join_style = 2, cap_style = 3) 
     for p in gs], crs = gs.crs) \
      .union_all().buffer(-RESOLUTION * 10, join_style = 2, cap_style = 3)
  if as_polygon:
    return gridify(union)
  else:
    return gridify(gpd.GeoSeries([union], crs = gs.crs))


def zigzag_between_points(
    p0:geom.Point, p1: geom.Point, n:int, h:float = 1.0,  smoothness: int = 0):

  template_steps = n * 2 + 1
  r = p0.distance(p1)
  
  x = np.linspace(0, n * np.pi, template_steps, endpoint = True)
  y = [np.sin(x) for x in x]
  s = interpolate.InterpolatedUnivariateSpline(x, y, k = 2)

  spline_steps = (n + smoothness) * 2 + 1
  xs = np.linspace(0, n * np.pi, spline_steps, endpoint = True)
  ys = s(xs)
  
  sfx = 1 / max(x) * r
  sfy = h * r / 2
  theta = np.arctan2(p1.y - p0.y, p1.x - p0.x)

  ls = geom.LineString([geom.Point(x, y) for x, y in zip(xs, ys)])
  ls = affine.translate(ls, 0, -(ls.bounds[1] + ls.bounds[3]) / 2)
  ls = affine.scale(ls, xfact = sfx, yfact = sfy, origin = (0, 0))
  ls = affine.rotate(ls, theta, (0, 0), use_radians = True)
  x0, y0 = list(ls.coords)[0]
  return affine.translate(ls, p0.x - x0, p0.y - y0)


def get_translation_transform(dx:float, dy:float) -> list[float]:
  """Returns the shapely affine transform tuple for a translation.

  Args:
      dx (float): translation distance in x direction.
      dy (float): translation distance in y direction.

  Returns:
    list[float]: a six item list of floats, per the shapely.affinity.
    affine_transform method, see 
      https://shapely.readthedocs.io/en/stable/manual.html#affine-transformations
  """
  return [1, 0, 0, 1, dx, dy]


def get_rotation_transform(
    angle:float, centre:tuple[float] = None) -> list[float]:
  """Returns the shapely affine transform tuple for a rotation, optionally
  about a supplied centre point.

  Args:
      angle (float): the angle of rotation (in degrees).
      centre (tuple[float], optional): An option centre location. Defaults to 
      None, which will in turn be converted to (0, 0).

  Returns:
    list[float]: a six item list of floats, per the shapely.affinity.
      affine_transform method, see 
        https://shapely.readthedocs.io/en/stable/manual.html#affine-transformations
  """
  if centre is None or np.allclose((0, 0), centre, atol = RESOLUTION):
    a = np.radians(angle)
    return [np.cos(a), -np.sin(a), np.sin(a), np.cos(a), 0, 0]
  
  dx, dy = centre
  t1 = get_translation_transform(-dx, -dy)
  r = get_rotation_transform(angle)
  t2 = get_translation_transform(dx, dy)
  return combine_transforms([t1, r, t2])


def get_reflection_transform(
    angle:float, centre:tuple[float] = None) -> list[float]:
  """Returns a shapely affine transform tuple that will reflect a shape
  in a line at the specified angle, optionally through a specified centre
  point.

  Args:
    angle (float): angle to the x-axis of the line of reflection.
    centre (tuple[float], optional): point through which the line of 
      reflection passes. Defaults to None, which
      will in turn be converted to (0, 0).

  Returns:
    list[float]: a six item list of floats, per the shapely.affinity.
      affine_transform method, see 
        https://shapely.readthedocs.io/en/stable/manual.html#affine-transformations
  """
  if centre is None or np.allclose((0, 0), centre, atol = RESOLUTION):
    A = 2 * np.radians(angle)
    return (np.cos(A), np.sin(A), np.sin(A), -np.cos(A), 0, 0)
  dx, dy = centre
  t1 = get_translation_transform(-dx, -dy)
  r = get_reflection_transform(angle)
  t2 = get_translation_transform(dx, dy)
  return combine_transforms([t1, r, t2])


def combine_transforms(transforms:list[list[float]]) -> list[float]:
  """Returns a shapely affine transform list that combines the listed
  sequence of transforms applied in order.

  Args:
    transforms (list[list[float]]): sequence of transforms to combine.

  Returns:
    list[float]: a transform tuple combining the supplied transforms applied
      in order, see 
        https://shapely.readthedocs.io/en/stable/manual.html#affine-transformations
  """
  result = np.identity(3)
  for t in transforms:
    result = as_numpy_matrix(t) @ result
  return as_shapely_transform(result)


def reverse_transform(transform:list[float]) -> list[float]:
  """Returns the inverse shapely affine transform of the supplied transform.

  Args:
    transform (list[float]): the transform for which the inverse is desired.

  Returns:
    list[float]: shapely affine transform tuple that will invert the supplied
      transform.
  """
  return as_shapely_transform(
    np.linalg.inv(as_numpy_matrix(transform)))


def as_shapely_transform(arr:np.array) -> list[float]:
  """Returns the shapely affine transform list equivalent to the supplied
  numpy matrix of a conventional augmented affine transform matrix.

  Args:
    arr (np.array): augmented affine transform matrix of the desired 
      transform.

  Returns:
    list[float]: desired shapely affine transform list of floats.
  """
  return [arr[0][0], arr[0][1], arr[1][0], arr[1][1], arr[0][2], arr[1][2]]


def as_numpy_matrix(transform:list[float]) -> np.array:
  """Converts the supplied shapely affine transform list to an augmented
  affine transform matrix in numpy array form. This makes combining transforms
  much easier.

  Args:
    transform (list[float]): the transform in shapely format.

  Returns:
    np.array: the transform in numpy matrix format.
  """
  return np.array([[transform[0], transform[1], transform[4]],
                   [transform[2], transform[3], transform[5]],
                   [           0,            0,            1]])


StraightLine = namedtuple("StraightLine", "A B C")
"""Named tuple to represent equation of a straight line in standard 
Ax + By + C = 0 form."""


def get_straight_line(p1:geom.Point, p2:geom.Point, 
                      perpendicular:bool = False) -> StraightLine:
  if perpendicular:
    ls = affine.rotate(geom.LineString([p1, p2]), 90)
    pts = [p for p in ls.coords]
    x1, y1 = pts[0]
    x2, y2 = pts[1]
  else:
    x1, y1 = p1.x, p1.y
    x2, y2 = p2.x, p2.y
  return StraightLine(y1 - y2, x2 - x1, x1 * y2 - x2 * y1)


def get_straight_line(p1:geom.Point, p2:geom.Point, 
                      perpendicular:bool = False) -> StraightLine:
  if perpendicular:
    ls = affine.rotate(geom.LineString([p1, p2]), 90)
    pts = [p for p in ls.coords]
    x1, y1 = pts[0]
    x2, y2 = pts[1]
  else:
    x1, y1 = p1.x, p1.y
    x2, y2 = p2.x, p2.y
  return StraightLine(y1 - y2, x2 - x1, x1 * y2 - x2 * y1)


def get_intersection(line1:StraightLine, 
                     line2:StraightLine) -> geom.Point:
  x_set, y_set = False, False
  denominator = line1.A * line2.B - line2.A * line1.B
  if np.isclose(line1.A, 0, atol = 1e-4, rtol = 1e-4):
    y = -line1.C / line1.B
    y_set = True
  elif np.isclose(line2.A, 0, atol = 1e-4, rtol = 1e-4):
    y = -line2.C / line2.B
    y_set = True
  if np.isclose(line1.B, 0, atol = 1e-4, rtol = 1e-4):
    x = -line1.C / line1.A
    x_set = True
  elif np.isclose(line2.B, 0, atol = 1e-4, rtol = 1e-4):
    x = -line2.C / line2.A
    x_set = True
  if np.isclose(denominator, 0, atol = 1e-4, rtol = 1e-4):
    return None
  x = x if x_set else (line1.B * line2.C - line2.B * line1.C) / denominator
  y = y if y_set else (line1.C * line2.A - line2.C * line1.A) / denominator
  return geom.Point(x, y)


def get_prototile_from_vectors(vecs:list[tuple[float]]):
  if len(vecs) == 2:
    v1, v2 = vecs
    v1 = [_ / 2 for _ in v1]
    v2 = [_ / 2 for _ in v2]
    return geom.Polygon([( v1[0] + v2[0],  v1[1] + v2[1]),
                         ( v2[0] - v1[0],  v2[1] - v1[1]),
                         (-v1[0] - v2[0], -v1[1] - v2[1]),
                         (-v2[0] + v1[0], -v2[1] + v1[1])])
  else:
    v1, v2, v3 = vecs
    v4 = [-_ for _ in v1]
    v5 = [-_ for _ in v2]
    v6 = [-_ for _ in v3]
    q1 = geom.Polygon([v1, v2, v4, v5]) 
    q2 = geom.Polygon([v2, v3, v5, v6]) 
    q3 = geom.Polygon([v3, v4, v6, v1])
    return gridify(q3.intersection(q2).intersection(q1))
