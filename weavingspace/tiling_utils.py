#!/usr/bin/env python
# coding: utf-8

"""Utility functions for use by various classes in the `weavingspace` 
package. Most are geometric convenience functions for commonly applied 
operations.
"""

from typing import Iterable, Union
import re
import string

from math import fsum
from math import isclose

import numpy as np

import geopandas as gpd
import pandas as pd
import shapely.geometry as geom
import shapely.affinity as affine
import shapely.wkt as wkt
import shapely.ops

# this causes a circular import error and is only needed for typing
# from weavingspace.tile_unit import TileUnit


def _parse_strand_label(s:str) -> list[str]:
    """Breaks a strand label specifiction in to a list of labels

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
    """Strands label string to split into strand labels.

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
    strand_ids = (strand_ids
                  if len(strand_ids) == 3
                  else strand_ids + [[""]])
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
    shape_c = shape.centroid.coords[0]
    return (target[0] - shape_c[0], target[1] - shape_c[1])


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
    return geom.Polygon(corners)


def is_convex(poly:geom.Polygon) -> bool:
    return isclose(poly.area, poly.convex_hull.area)


def incentre(poly:geom.Polygon) -> geom.Point:
    """A different polygon centre, which produces better results for some
    of the dual tilings where polygons are not regular... see
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
    
    It is possible that intersecting two angle bisectors of the polygon would 
    also work, but it is unclear if this is an equivalent construction (more
    reading required...) 
    
        Args:
        poly (geom.Polygon): the polygon.

    Returns:
        geom.Point: the incentre of the polygon.
    """
    poly = ensure_ccw(poly)
    corners = [geom.Point(p[0], p[1]) 
               for p in list(poly.exterior.coords)]
    c = poly.centroid
    # if corners are equidistant from the centroid, then return that
    d0 = corners[0].distance(c)
    if all([isclose(pt.distance(c), d0) for pt in corners[0:-1]]):
        return c
    else:  # find the incentre
        r = get_apothem(poly)
        # NOTE: for some reason a simple negative buffer here does not work
        e1 = geom.LineString(corners[:2]).parallel_offset(r, side = "left")
        e2 = geom.LineString(corners[1:3]).parallel_offset(r, side = "left")
        return e1.intersection(e2)  # TODO: add exception if no intersection

        # b1 = get_angle_bisector(poly, 0)
        # b2 = get_angle_bisector(poly, 1)
        # return b1.intersection(b2)


def get_apothem(poly:geom.Polygon) -> float:
    return 2 * poly.area / geom.LineString(poly.exterior.coords).length


def ensure_ccw(poly:geom.Polygon) -> geom.Polygon:
    if not geom.LinearRing(poly.exterior.coords).is_ccw:
        return geom.Polygon(list(reversed(poly.exterior.coords))[:-1])
    else:
        return poly


def get_angle_bisector(poly:geom.Polygon, v = 0) -> geom.LineString:
    pts = [geom.Point(p) for p in poly.exterior.coords][:-1]
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
    return affine.scale(geom.LineString([p0, m]), 100, 100, origin = p0)


def _get_interior_vertices(polys:gpd.GeoDataFrame) -> gpd.GeoSeries:
    """Returns points not on the outer boundary of the supplied set 
    of polygons.

    Args:
        polys (gpd.GeoDataFrame): a set of polygons.

    Returns:
        gpd.GeoSeries: interior vertices of the set of polygons.
    """
    polygons = gridify(polys.geometry)
    uu = safe_union(polygons, as_polygon = True)
    interior_pts = set()
    for poly in polygons:
        for pt in poly.exterior.coords:
            if uu.contains(
                geom.Point(pt).buffer(1e-3, resolution = 1, join_style = 2)):
                interior_pts.add(pt)
    return gpd.GeoSeries([geom.Point(p) for p in interior_pts])

    
def gridify(gs:gpd.GeoSeries, precision:int = 6) -> gpd.GeoSeries:
    """Returns the supplied GeoSeries rounded to the specified precision.
    
    Works by round-tripping through WKT, which seems like the easiest
    way to do this given the variety of ways in which coordinates are
    stored depending on the geometry.

    Args:
        gs (gpd.GeoSeries): geometries to gridify.
        precision (int, optional): digits of precision. Defaults to 6.

    Returns:
        gpd.GeoSeries: the rounded geometries.
    """
    return gpd.GeoSeries(
        list(gs.apply(
            wkt.dumps, rounding_precision = precision).apply(wkt.loads)))

    
def get_dual_tile_unit(t) -> gpd.GeoDataFrame:
    """Converts supplied TileUnit to a candidate GeoDataFrame of its dual 
    TileUnit.
    
    NOTE: this is complicated and not remotely guaranteed to work!
    a particular issue is that where to place the vertices of the faces
    of the dual with respect to the tiles in the original is ill-defined. 
    This is because the dual process is topologically not metrically defined, 
    so that exact vertex locations are ambiguous. Tiling duality is defined in 
    Section 4.2 of Grunbaum B, Shephard G C, 1987 _Tilings and Patterns_ (W. H. 
    Freeman and Company, New York)
    
    NOTE: In general, this method will work only if all supplied elements are 
    regular polygons. A known exception is if the only non-regular polygons are 
    triangles.
    
    NOTE: 'clean' polygons are required. If supplied polygons have messy 
    vertices with multiple points where there is only one proper point, bad 
    things are likely to happen! Consider using `clean_polygon()` on the
    element geometries.
    
    Because of the above limitations, we only return a GeoDataFrame 
    for inspection. However some `weavingspace.tile_unit.TileUnit` setup 
    methods in `weavingspace.tiling_geometries` use this method, where we are 
    confident the returned dual is valid.    

    Args:
        t (TileUnit): the tiling for which the dual is required.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame that could be the elements attribute for
            a TileUnit of the dual tiling.
    """
    # get a local patch of this Tiling
    local_patch = t.get_local_patch(r = 3, include_0 = True)     
    # Find the interior points of these tiles - these will be guaranteed
    # to have a sequence of surrounding tiles incident on them 
    interior_pts = _get_interior_vertices(local_patch)
    # Compile a list of the polygons incident on the interior points
    cycles = []
    for pt in interior_pts:
        cycles.append(
            set([poly_id for poly_id, p in enumerate(local_patch.geometry) 
                 if pt.distance(p) < t.fudge_factor]))
    # convert the polygon ID sequences to (centroid.x, centroid.y, ID) tuples
    dual_faces = []
    for cycle in cycles:
        ids, pts = [], []
        for poly_id in cycle:
            ids.append(local_patch.element_id[poly_id])
            poly = local_patch.geometry[poly_id]
            pts.append(incentre(poly))
        # sort them into CCW order so they are well formed
        sorted_coords = sort_ccw([(pt.x, pt.y, id) 
                                  for pt, id in zip(pts, ids)])
        dual_faces.append(
            (geom.Polygon([(pt_id[0], pt_id[1]) for pt_id in sorted_coords]), 
             "".join([pt_id[2] for pt_id in sorted_coords])))
    # ensure the resulting face centroids are inside the original tile
    # displaced a little to avoid uncertainties at corners/edges
    dual_faces = [(f, id) for f, id in dual_faces
                  if affine.translate(t.tile.geometry[0],
                                      t.fudge_factor, 
                                      t.fudge_factor).contains(f.centroid)]
    gdf = gpd.GeoDataFrame(
        data = {"element_id": [f[1] for f in dual_faces]}, crs = t.crs,
        geometry = gpd.GeoSeries([f[0] for f in dual_faces]))
    # ensure no duplicates
    gdf = gdf.dissolve(by = "element_id", as_index = False).explode(
        index_parts = False, ignore_index = True)
    
    gdf.element_id = relabel(gdf.element_id)
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
        

def sort_ccw(pts_ids:list[tuple[float, float, str]]):
    """Sorts supplied tuple of x, y, ID into counterclockwise order.

    Args:
        pts_ids (list[tuple[float, float, str]]): A tuple of a pair of 
            floats and a string.

    Returns:
        list: a list in the same format as supplied sorted into 
            counter-clockwise order of the point locations.
    """
    x = [p[0] for p in pts_ids]
    y = [p[1] for p in pts_ids]
    cx, cy = np.mean(x), np.mean(y)
    dx = [_ - cx for _ in x]
    dy = [_ - cy for _ in y]
    angles = [np.arctan2(dy, dx) for dx, dy in zip(dx, dy)]
    d = dict(zip(angles, pts_ids))
    return [pt_id for angle, pt_id in sorted(d.items())]


def write_map_to_layers(gdf:gpd.GeoDataFrame, fname:str = "output.gpkg",
                        element_var:str = "element_id") -> None:
    """Writes supplied GeoDataFrame to a GPKG file with layers based on 
    the element_var attribute.

    Args:
        gdf (gpd.GeoDataFrame): the GeoDataFrame.
        fname (str, optional): filename to write.
        element_var (str, optional): the attribute to use to separate
            output file into layers. Defaults to "element_id".
    """
    grouped = gdf.groupby(element_var, as_index = False)
    for e in pd.Series.unique(gdf[element_var]):
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
        return get_apothem(geometry)
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
    areas = list(polygons.area)
    max_area = max(areas)
    return gpd.GeoSeries([polygons[areas.index(max_area)]])


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
    return p1.buffer(1e-3, resolution = 1, join_style = 2).intersection(
        p2.buffer(1e-3, resolution = 1, join_style = 2)).area > 1e-5


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


def get_bounding_ellipse(
        shapes:gpd.GeoSeries, mag:float = 1.0) -> gpd.GeoSeries:
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
    return gpd.GeoSeries(circle, crs = shapes.crs).scale(
        w / r * mag / np.sqrt(2), h / r * mag / np.sqrt(2), origin = c)
    
    
def get_boundaries(shapes:gpd.GeoSeries) -> gpd.GeoSeries:
    """Returns linestring GeoSeries from supplied polygon GeoSeries.

    This is used to allow display of edges of tiles in legend when they are 
    masked by an ellipse (if we instead clip polygons then the ellipse edge 
    will also show in the result.)
    
    Args:
        shapes (gpd.GeoSeries): Polygons to convert.

    Returns:
        gpd.GeoSeries: LineStrings from the supplied Polygons.
    """
    bdys = [geom.LineString(p.exterior.coords) for p in shapes]
    return gpd.GeoSeries(bdys, crs = shapes.crs)


def get_polygon_sector(shape:geom.Polygon, start:float = 0.0, 
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
        arc1 = shapely.ops.substring(geom.LineString(shape.exterior.coords), 
                                     np.mod(e1, 1), 1, normalized = True) 
        arc2 = shapely.ops.substring(geom.LineString(shape.exterior.coords), 
                                     0, e2, normalized = True) 
        sector = geom.Polygon([shape.centroid] + 
                              list(arc1.coords) + 
                              list(arc2.coords)[1:])
    else:
        arc = shapely.ops.substring(geom.LineString(shape.exterior.coords), 
                                    start, end, normalized = True)
        sector = geom.Polygon([shape.centroid] + list(arc.coords))
    return clean_polygon(sector)


def clean_polygon(p:Union[geom.Polygon, gpd.GeoSeries], 
                      res:float = 1e-3, shrink_then_grow:bool = True
                  ) -> Union[geom.Polygon, gpd.GeoSeries]:
    """Convenience function to 'clean' a shapely polyon or GeoSeries by applying
    a negative buffer then the same positive buffer.

    Optionally the buffer may be applied in the opposite order (i.e. grow then 
    shrink)
        
    This is a procedure often unofficially recommended (on stackexchange etc.) 
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
    if shrink_then_grow:
        return p.buffer(-res, resolution = 1, join_style = 2).buffer(
                res, resolution = 1, join_style = 2)
    else:
        return p.buffer(res, resolution = 1, join_style = 2).buffer(
                -res, resolution = 1, join_style = 2)
    

def safe_union(gs:gpd.GeoSeries, res:float = 1e-3, 
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
    union = gs.buffer(res, resolution = 1, join_style = 2) \
                .unary_union \
                .buffer(-res, resolution = 1, join_style = 2) 
    if as_polygon:
        return union
    else:
        return gpd.GeoSeries([union], crs = gs.crs)