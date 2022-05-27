#!/usr/bin/env python
# coding: utf-8

import string
import geopandas as gpd
import pandas as pd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine
import shapely.wkt as wkt


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
    n_ids = len(set(list(gdf.element_id)))
    gdf.element_id = list(string.ascii_letters[:n_ids])
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

