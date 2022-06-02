#!/usr/bin/env python
# coding: utf-8

from enum import Enum
import re
import string
import geopandas as gpd
import pandas as pd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine
import shapely.wkt as wkt


class TileShape(Enum):
    RECTANGLE = "rectangle"
    HEXAGON = "hexagon"
    TRIANGLE = "triangle"
    DIAMOND = "diamond"


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
            each direction and () designates combining labels into a single strand that will be sliced lengthwise. Example output:

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


def get_axis_from_label(label:str = "a", strands:str = "a|b|c"):
    index = strands.index(label)
    return strands[:index].count("|")


def get_colour_ramp(geometry, n = 10, a = 0):
    c = geometry.centroid
    g = affine.rotate(geometry, -a, origin = c)
    bb = g.bounds
    cuts = np.linspace(bb[0], bb[2], n + 1)
    slices = []
    for l, r in zip(cuts[:-1], cuts[1:]):
        slice = geom.Polygon([(l, bb[1] - 1), (r, bb[1] - 1), 
                              (r, bb[3] + 1), (l, bb[3] + 1)])
        slices.append(slice.intersection(g)) 
    return [affine.rotate(s, a, origin = c) for s in slices]


def get_regular_polygon(spacing, n) -> geom.Polygon: 
    R = spacing / np.cos(np.radians(180 / n)) / 2
    a0 = -90 + 180 / n
    a_diff = 360 / n
    angles = [a0 + a * a_diff for a in range(n)]
    corners = [(R * np.cos(np.radians(a)), 
                R * np.sin(np.radians(a))) for a in angles]
    return geom.Polygon(corners)


def incentre(tri:geom.Polygon) -> geom.Point:
    corners = [geom.Point(p[0], p[1]) 
               for p in list(tri.exterior.coords)]  # ABCA
    lengths = [p.distance(q) 
               for p, q in zip(corners[:-1], corners[1:])]  # AB, BC, CA
    lengths = lengths[1:] + lengths[:1]  # BC, CA, AB 
    perimeter = np.sum(lengths)
    incentre = [0, 0]
    for p, l in zip(corners[:-1], lengths):
        incentre[0] = incentre[0] + p.x * l / perimeter
        incentre[1] = incentre[1] + p.y * l / perimeter
    return geom.Point(incentre)
    

def get_interior_vertices(polys:gpd.GeoDataFrame) -> gpd.GeoSeries:
    polygons = gridify(polys.geometry)
    uu = polygons.unary_union.buffer(1e-3).buffer(-1e-3)
    interior_pts = set()
    for poly in polygons:
        for pt in poly.exterior.coords:
            if uu.contains(geom.Point(pt).buffer(1e-3)):
                interior_pts.add(pt)
    return gpd.GeoSeries([geom.Point(p) for p in interior_pts])

    
def gridify(gs, precision = 6) -> gpd.GeoSeries:
    return gpd.GeoSeries(
        list(gs.apply(
            wkt.dumps, rounding_precision = precision).apply(wkt.loads)))

    
# Converts the supplied TileUnit to a candidate GeoDataFrame of its dual
# BUT NOTE: this is complicated and not remotely guaranteed to work!
# a particular issue is that where to place the vertices of the faces
# of the dual with respect to the tiles in the original is ill-defined.
# This is essentially because the dual process is not metrically defined
# only topologically, so the vertex locations are arbitrary.
# Only return a GeoDataFrame because there is no easy or obvious way to 
# infer the output's tiling geometry (which we would need for a TileUnit)
def get_dual_tile_unit(t) -> gpd.GeoDataFrame:
    # get a local patch for a 3x scaled tile extent of this Tiling
    local_patch = t.get_local_patch(r = 3, include_0 = True)     
    # Find the interior points of these tiles - these will be guaranteed
    # to have a sequence of surrounding tiles incident on them 
    interior_pts = get_interior_vertices(local_patch)
    # Compile a list of the polygons incident on the interior points
    cycles = []
    for pt in interior_pts:
        cycles.append(
            set([i for i, p in enumerate(local_patch.geometry) 
                 if pt.distance(p) < t.fudge_factor]))
    # These can be used to construct the dual polygons
    dual_faces = []
    ids = []
    for cycle in cycles:
        id = []
        coords = []
        for i in cycle:
            id.append(local_patch.element_id[i])
            poly = local_patch.geometry[i]
            if len(poly.exterior.coords) == 4:
                centroid = incentre(poly)
            else:
                centroid = poly.centroid
            coords.append([centroid.x, centroid.y])
        # sort them into CCW order so they are well formed
        sorted_coords = sort_ccw([(p, i) for p, i in zip(coords, id)])
        dual_faces.append(geom.Polygon([pt_id[0] for pt_id in sorted_coords]))
        # a reasonable stab at element IDs is the sequence of element_id
        # values from the original tiling in the order encountered
        ids.append("".join([pt_id[1] for pt_id in sorted_coords]))
    # ensure the resulting faces actually contain the centroids of
    # the generating polygon...
    dual_faces = [(f, id) for f, id in zip(dual_faces, ids)
                  if affine.translate(t.tile.geometry[0],
                                      t.fudge_factor, 
                                      t.fudge_factor).contains(f.centroid)]
    gdf = gpd.GeoDataFrame(
        data = {"element_id": [f[1] for f in dual_faces]}, crs = t.crs,
        geometry = gpd.GeoSeries([f[0] for f in dual_faces]))
    # ensure no duplicates
    gdf = gdf.dissolve(by = "element_id", as_index = False).explode(
        index_parts = False, ignore_index = True)
    new_ids = {}
    id_count = 0
    for id in gdf.element_id:
        if not id in new_ids:
            new_ids[id] = string.ascii_letters[id_count]
            id_count = id_count + 1
    gdf.element_id = [new_ids[id] for id in gdf.element_id]
    return gdf

    
# sort supplied  points into CCW order - by measuring angles around
# their centroid
def sort_ccw(pts_ids):
    x = [p[0][0] for p in pts_ids]
    y = [p[0][1] for p in pts_ids]
    cx, cy = np.mean(x), np.mean(y)
    dx = [_ - cx for _ in x]
    dy = [_ - cy for _ in y]
    angles = [np.arctan2(dy, dx) for dx, dy in zip(dx, dy)]
    d = dict(zip(angles, pts_ids))
    return [pt_id for amgle, pt_id in sorted(d.items())]


def write_map_to_layers(gdf, fname = "output.gpkg", element_var = "element_id"):
    grouped = gdf.groupby(element_var, as_index = False)
    # for e in gdf[element_var].drop_duplicates():
    for e in pd.Series.unique(gdf[element_var]):
        grouped.get_group(e).to_fule(fname, layer = e, driver = "GPKG")
        

def get_insets(geometry:geom.Polygon, n = 25, equal_areas = True):
    radius = get_collapse_distance(geometry)
    bandwidths = range(1, n + 2) if equal_areas else [1] * (n + 1)
    distances = np.cumsum(bandwidths)
    distances = distances * radius / distances[-1]
    nested_geoms = [geometry.buffer(-d, join_style = 2) 
                    for d in distances]
    return [g1.difference(g2) 
            for g1, g2 in zip(nested_geoms[:-1], nested_geoms[1:])]


def get_collapse_distance(geometry):
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
