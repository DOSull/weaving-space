#!/usr/bin/env python
# coding: utf-8

from enum import Enum
from typing import Iterable, Union
import re
import string

import matplotlib
import matplotlib.colors
import numpy as np

import geopandas as gpd
import pandas as pd
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


def get_colour_ramp(geometry:geom.Polygon, 
                    n:int = 25, a:float = 0) -> list[geom.Polygon]:
    """Returns a series of 'slices' of the supplied polygon so that it
    can be used as a legend colour ramp. 
    
    NOTE: Perhaps belongs in WeaveUnit?

    Args:
        geometry (geom.Polygon): the polygon (usually a weave element)
        n (int, optional): the number of slices. Defaults to 25.
        a (float, optional): orientation of the polygon, it will sliced
            across this orientation. Defaults to 0.

    Returns:
        list[geom.Polygon]: list of geom.Polygons.
    """
    c = geometry.centroid
    g = affine.rotate(geometry, -a, origin = c)
    bb = g.bounds
    cuts = np.linspace(bb[0], bb[2], n + 1)
    # add a margin for error, since intersection operations can cause odd
    # effects when parallel lines or other near-misses are involved...
    # First, the y-axis
    bb[1] = bb[1] - 1
    bb[3] = bb[3] + 1
    # and then the x-axis
    cuts[0] = cuts[0] - 1
    cuts[-1] = cuts[-1] + 1
    slices = []
    for l, r in zip(cuts[:-1], cuts[1:]):
        slice = geom.Polygon([(l, bb[1]), (r, bb[1]), (r, bb[3]), (l, bb[3])])
        slices.append(slice.intersection(g)) 
    return [affine.rotate(s, a, origin = c) for s in slices]


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


def incentre(tri:geom.Polygon) -> geom.Point:
    """A different triangle centre, which produces better results for some
    of the dual tilings were triangles are not equilateral... see
    https://en.wikipedia.org/wiki/Incenter

    Args:
        tri (geom.Polygon): the triangle.

    Returns:
        geom.Point: the incentre of the triangle.
    """
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

    
def get_dual_tile_unit(t) -> gpd.GeoDataFrame:
    """Converts supplied TileUnit to a candidate GeoDataFrame of its dual 
    TileUnit.
    
    NOTE: this is complicated and not remotely guaranteed to work!
    a particular issue is that where to place the vertices of the faces
    of the dual with respect to the tiles in the original is ill-defined.
    This is because the dual process is topologically not metrically defineed, 
    so that exact vertex locations are ambiguous. 
    
    We therefore only return a GeoDataFrame for inspection. 
    
    However some of TileUnit setup methods in tiling_geometries.py use this
    method, where we are confident the dual TileUnit returned is valid.    

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
    interior_pts = get_interior_vertices(local_patch)
    # Compile a list of the polygons incident on the interior points
    cycles = []
    for pt in interior_pts:
        cycles.append(
            set([poly_id for poly_id, p in enumerate(local_patch.geometry) 
                 if pt.distance(p) < t.fudge_factor]))
    # convert the polygon ID sequences to (centroid, ID) tuples
    dual_faces = []
    for cycle in cycles:
        ids, coords = [], []
        for poly_id in cycle:
            ids.append(local_patch.element_id[poly_id])
            poly = local_patch.geometry[poly_id]
            if len(poly.exterior.coords) == 4:
                centroid = incentre(poly)
            else:
                centroid = poly.centroid
            coords.append([centroid.x, centroid.y])
        # sort them into CCW order so they are well formed
        sorted_coords = sort_ccw([(p, i) for p, i in zip(coords, ids)])
        dual_faces.append(
            (geom.Polygon([pt_id[0] for pt_id in sorted_coords]), 
             "".join([pt_id[1] for pt_id in sorted_coords])))
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
        

    
# sort supplied  points into CCW order - by measuring angles around
# their centroid
def sort_ccw(pts_ids:list[tuple[tuple[float, float], str]]):
    """Sorts supplied tuple of a pair of floats (points) and an ID into
    counterclockwise order.

    Args:
        pts_ids (list[tuple[tuple[float, float], str]]): A tuple of a pair of 
            floats and a string.

    Returns:
        list: a list in the same format as supplied sorted into 
            counter-clockwise order.
    """
    x = [p[0][0] for p in pts_ids]
    y = [p[0][1] for p in pts_ids]
    cx, cy = np.mean(x), np.mean(y)
    dx = [_ - cx for _ in x]
    dy = [_ - cy for _ in y]
    angles = [np.arctan2(dy, dx) for dx, dy in zip(dx, dy)]
    d = dict(zip(angles, pts_ids))
    return [pt_id for angle, pt_id in sorted(d.items())]


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


def plot_subsetted_gdf(ax, gdf, vars = {}, pals = {}, **kwargs):
    ids = pd.Series.unique(gdf.element_id)
    groups = gdf.groupby("element_id")
    for id in ids:
        subset = groups.get_group(id)
        # Handle custom color assignments via 'pals' parameter.
        # Result is setting 'cmap' variable used in plot command afterwards.
        if (isinstance(pals, 
                        (str, matplotlib.colors.Colormap,
                        matplotlib.colors.LinearSegmentedColormap,
                        matplotlib.colors.ListedColormap))):
            cmap=pals  # user wants one palette for all ids
        elif (len(pals) == 0):
            cmap = 'Reds'  # set a default... here, to Brewer's 'Reds'
        elif (id not in pals):
            cmap = 'Reds'  # id has no color specified in dict, use default
        elif (isinstance(pals[id], 
                            (str, matplotlib.colors.Colormap,
                            matplotlib.colors.LinearSegmentedColormap,
                            matplotlib.colors.ListedColormap))):
            cmap = pals[id]  # user specified colors for this id so use it
        elif (isinstance(pals[id], dict)):
            colormap_dict = pals[id]
            data_unique_sorted = subset[vars[id]].unique()
            data_unique_sorted.sort()
            cmap = matplotlib.colors.ListedColormap(
                [colormap_dict[x] for x in data_unique_sorted])
        else:
            raise Exception(f"Color map for '{id}' is not a known type, but is {str(type(pals[id]))}")

        subset.plot(ax = ax, column = vars[id], cmap=cmap, **kwargs)
    return ax


def get_largest_polygon(polygons:gpd.GeoSeries) -> gpd.GeoSeries:
    areas = list(polygons.area)
    max_area = max(areas)
    return gpd.GeoSeries([polygons[areas.index(max_area)]])
