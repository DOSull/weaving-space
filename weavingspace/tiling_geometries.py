#!/usr/bin/env python
# coding: utf-8

from typing import TYPE_CHECKING
import copy
import itertools
import string

import geopandas as gpd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine

import tiling_utils
from tiling_utils import TileShape

if TYPE_CHECKING:
    from tile_units import TileUnit


def setup_none_tile(unit:"TileUnit") -> None: 
    setup_base_tile(unit, unit.tile_shape)
    unit.elements = gpd.GeoDataFrame(
        data = {"element_id": ["a"]}, crs = unit.crs,
        geometry = copy.deepcopy(unit.tile.geometry))


def setup_base_tile(unit, shape:TileShape) -> None: 
    unit.tile_shape = shape
    tile = tiling_utils.get_regular_polygon(
        unit.spacing, n = (4 
                        if unit.tile_shape in (TileShape.RECTANGLE, )
                        else (6 
                            if unit.tile_shape in (TileShape.HEXAGON, )
                            else 3)))
    unit.tile = gpd.GeoDataFrame(
        geometry = gpd.GeoSeries([tile]), crs = unit.crs)
    return


def setup_cairo(unit) -> None:
    setup_base_tile(unit, TileShape.RECTANGLE)
    d = unit.spacing        
    x = d / 2 / (np.cos(np.radians(15)) + np.cos(np.radians(75)))
    p1 = geom.Polygon([(x, 0),
                        (0, 0),
                        (0, x),
                        (x * np.sqrt(3) / 2, x + x / 2),
                        (x * (1 + np.sqrt(3)) / 2, 
                        x * (3 - np.sqrt(3)) / 2)])
    p1 = affine.rotate(p1, -15, (0, 0))
    p2 = affine.rotate(p1, 90, (0, 0))
    p3 = affine.rotate(p2, 90, (0, 0))
    p4 = affine.rotate(p3, 90, (0, 0))

    unit.elements = gpd.GeoDataFrame(
        data = {"element_id": list("abcd")}, crs = unit.crs,
        geometry = gpd.GeoSeries([p1, p2, p3, p4]))
    unit.regularised_tile = copy.deepcopy(unit.tile)
    unit.regularised_tile.geometry = gpd.GeoSeries([
        unit.elements.geometry.buffer(1e-3, join_style = 2) \
            .unary_union.buffer(-1e-3)])
    
        
    
def setup_hex_dissection(unit):
    setup_base_tile(unit, TileShape.HEXAGON)
    hex = tiling_utils.get_regular_polygon(unit.spacing, 6)
    v = list(hex.exterior.coords)
    m = [((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
                    for p1, p2 in zip(v[:-1], v[1:])]
    p = list(itertools.chain(*zip(v[:-1], m)))
    if unit.n == 2:
        s = (geom.Polygon([p[0], p[2], p[4], p[6], (0, 0)])
                if unit.dissection_offset == 0
                else geom.Polygon(p[1:8] + [(0, 0)]))
        slices = [affine.rotate(s, a, origin = (0, 0)) 
                    for a in range(0, 360, 180)]
    elif unit.n == 3:
        s = (geom.Polygon([p[0], p[2], p[4], (0, 0)])
                if unit.dissection_offset == 0
                else geom.Polygon([p[1], p[2], p[4], p[5], (0, 0)]))
        slices = [affine.rotate(s, a, origin = (0, 0)) 
                    for a in range(0, 360, 120)]
    elif unit.n == 4:
        s1 = (geom.Polygon([p[0], p[2], p[3], (0, 0)])
                if unit.dissection_offset == 0
                else geom.Polygon([p[1], p[2], p[4], (0, 0)]))
        s2 = (geom.Polygon([p[3], p[4], p[6], (0, 0)])
                if unit.dissection_offset == 0
                else geom.Polygon([p[4], p[6], p[7], (0, 0)]))
        slices = [s1, s2, 
                    affine.rotate(s1, 180, (0, 0)),
                    affine.rotate(s2, 180, (0, 0))]
    elif unit.n == 6:
        s = (geom.Polygon([p[0], p[2], (0, 0)])
                if unit.dissection_offset == 0
                else geom.Polygon([p[1], p[2], p[3], (0, 0)]))
        slices = [affine.rotate(s, a, origin = (0, 0)) 
                    for a in range(0, 360, 60)]
    elif unit.n == 12:
        s1 = geom.Polygon([p[0], p[1], (0, 0)])
        s2 = geom.Polygon([p[1], p[2], (0, 0)])
        slices1 = [affine.rotate(s1, a, origin = (0, 0)) 
                    for a in range(0, 360, 60)]
        slices2 = [affine.rotate(s2, a, origin = (0, 0)) 
                    for a in range(0, 360, 60)]
        slices = itertools.chain(slices1, slices2)
    
    unit.elements = gpd.GeoDataFrame(
        data = {"element_id": list(string.ascii_lowercase)[:unit.n]}, 
        crs = unit.crs,
        geometry = gpd.GeoSeries(slices))


def setup_laves(unit) -> None:
    if unit.code == "3.3.3.3.3.3":
        setup_base_tile(unit, TileShape.HEXAGON)
        setup_none_tile()
        return
    if unit.code == "3.3.3.3.6":
        setup_laves_33336(unit)
        return
    elif unit.code == "3.3.3.4.4":
        unit.n = 2
        unit.dissection_offset = 1
        setup_hex_dissection(unit)
        return
    elif unit.code == "3.3.4.3.4":
        setup_cairo(unit)
        return
    elif unit.code == "3.4.6.4":
        unit.n = 6
        unit.dissection_offset = 1
        setup_hex_dissection(unit)
        return
    elif unit.code == "3.6.3.6":
        unit.n = 3
        unit.dissection_offset = 0
        setup_hex_dissection(unit)
        return
    elif unit.code == "3.12.12":
        setup_laves_31212(unit)
        return
    elif unit.code == "4.4.4":
        setup_base_tile(unit, TileShape.RECTANGLE)
        setup_none_tile()
        return
    elif unit.code == "4.6.12":
        unit.n = 12
        setup_hex_dissection(unit)
        return
    elif unit.code == "4.8.8":
        setup_laves_488(unit)
        return
    elif unit.code == "6.6.6":
        setup_base_tile(unit, TileShape.TRIANGLE)
        setup_none_tile()
    else:
        print(f"[{unit.code}] is not a valid Laves code.")

    unit.tiling_type = None
    unit.setup_none_tile()
    return


def setup_laves_33336(unit) -> None:
    setup_base_tile(unit, TileShape.HEXAGON)
    offset_a = np.degrees(np.arctan(1 / 3 / np.sqrt(3)))
    sf = 1 / np.sqrt(7)
    tile = unit.tile.geometry[0]
    hex = affine.scale(tile, sf, sf)
    hex = affine.translate(hex, 0, hex.bounds[3] - hex.bounds[1])
    hex_p = [p for p in hex.exterior.coords]
    petal = geom.Polygon([(0, 0)] + hex_p[1:5])
    petals = [affine.rotate(petal, a + offset_a, origin = (0, 0))
                    for a in range(30, 360, 60)]
    unit.elements = gpd.GeoDataFrame(
        data = {"element_id": list("abcdef")}, 
        crs = unit.crs,
        geometry = gpd.GeoSeries(petals))
    unit.regularised_tile = copy.deepcopy(unit.tile)
    unit.regularised_tile.geometry = gpd.GeoSeries([
        unit.elements.geometry.buffer(1e-3, join_style = 2) \
            .unary_union.buffer(-1e-3)])
        

def setup_laves_488(unit) -> None:
    setup_base_tile(unit, TileShape.RECTANGLE)
    tile = unit.tile.geometry[0]
    pts = [p for p in tile.exterior.coords]
    tri1 = geom.Polygon([pts[0], pts[1], geom.Point(0, 0)])
    tri2 = affine.rotate(tri1, 90, (0, 0))
    tri3 = affine.rotate(tri2, 90, (0, 0))
    tri4 = affine.rotate(tri3, 90, (0, 0))
    unit.elements = gpd.GeoDataFrame(
        data = {"element_id": list("abcd")}, 
        crs = unit.crs,
        geometry = gpd.GeoSeries([tri1, tri2, tri3, tri4]))
    unit.regularised_tile = copy.deepcopy(unit.tile)
    unit.regularised_tile.geometry = gpd.GeoSeries([
        unit.elements.geometry.buffer(1e-3, join_style = 2) \
            .unary_union.buffer(-1e-3)])
        

def setup_laves_31212(unit):
    setup_base_tile(unit, TileShape.TRIANGLE)
    tri = unit.tile.geometry[0]
    pts = [p for p in tri.exterior.coords]
    t1 = geom.Polygon([pts[0], pts[1], tri.centroid])
    elements = [affine.rotate(t1, a, tri.centroid) 
                for a in range(0, 360, 120)]
    unit.elements = gpd.GeoDataFrame(
        data = {"element_id": list("abc")}, crs = unit.crs,
        geometry = gpd.GeoSeries(elements)
    )


def setup_archimedean(unit) -> None:
    if unit.code == "3.3.3.3.3.3":
        setup_none_tile(unit, TileShape.TRIANGLE)
        return
    if unit.code == "3.3.3.3.6":
        print(f"The code [{unit.code}] is unsupported.")
    elif unit.code == "3.3.3.4.4":
        print(f"The code [{unit.code}] is unsupported.")
    elif unit.code == "3.3.4.3.4":
        setup_laves(unit)
        unit.elements = tiling_utils.get_dual_tile_unit(unit)
        unit.regularised_tile = \
            gpd.GeoSeries([unit.elements.geometry.buffer(
                1e-3, 1).unary_union.buffer(-1e-3, 1)])
        return
    elif unit.code == "3.4.6.4":
        # our dual tiling code isn't correct for this one
        setup_archimedean_3464(unit)
        return
    elif unit.code == "3.6.3.6":
        setup_laves(unit)
        unit.elements = tiling_utils.get_dual_tile_unit(unit)
        unit.regularised_tile = copy.deepcopy(unit.tile)
        unit.regularised_tile.geometry = gpd.GeoSeries([
            unit.elements.geometry.buffer(1e-3, join_style = 2) \
                .unary_union.buffer(-1e-3)])
        return
    elif unit.code == "3.12.12":
        setup_laves(unit)
        unit._modify_tile()
        unit._modify_elements()
        unit.elements = tiling_utils.get_dual_tile_unit(unit)
        unit.regularised_tile = copy.deepcopy(unit.tile)
        unit.regularised_tile.geometry = gpd.GeoSeries([
            unit.elements.geometry.buffer(1e-3, join_style = 2) \
                .unary_union.buffer(-1e-3)])
        return
    elif unit.code == "4.4.4":
        setup_base_tile(unit, TileShape.RECTANGLE)
        setup_none_tile(unit)
        return
    elif unit.code == "4.6.12":
        setup_laves(unit)
        unit.elements = tiling_utils.get_dual_tile_unit(unit)
        unit.regularised_tile = copy.deepcopy(unit.tile)
        unit.regularised_tile.geometry = gpd.GeoSeries([
            unit.elements.geometry.buffer(1e-3, join_style = 2) \
                .unary_union.buffer(-1e-3)])
        return
    elif unit.code == "4.8.8":
        setup_laves(unit)
        unit.elements = tiling_utils.get_dual_tile_unit(unit)
        unit.regularised_tile = copy.deepcopy(unit.tile)
        unit.regularised_tile.geometry = gpd.GeoSeries([
            unit.elements.geometry.buffer(1e-3, join_style = 2) \
                .unary_union.buffer(-1e-3)])
        return
    elif unit.code == "6.6.6":
        setup_base_tile(unit, TileShape.HEXAGON)
        setup_none_tile(unit)
    else:
        print(f"[{unit.code}] is not a valid Laves code.")

    unit.tiling_type = None
    setup_none_tile(unit)
    return


def setup_archimedean_3464(unit) -> None:
    setup_base_tile(unit, TileShape.HEXAGON)
    sf = np.sqrt(3) / (1 + np.sqrt(3))
    tile = tiling_utils.get_regular_polygon(unit.spacing, 6)

    hex = affine.scale(tile, sf, sf)
    corners = [p for p in hex.exterior.coords]
    p1 = corners[1]
    p2 = corners[0]
    dx, dy = p2[0] - p1[0], p2[1] - p1[1]
    p3 = (p2[0] - dy, p2[1] + dx)
    p4 = (p3[0] - dx, p3[1] - dy)
    sq1 = geom.Polygon([p1, p2, p3, p4])
    sq2 = affine.rotate(sq1, 60, (0, 0))
    sq3 = affine.rotate(sq2, 60, (0, 0))
    p5 = [pt for pt in sq2.exterior.coords][2]
    tri1 = geom.Polygon([p1, p4, p5])
    tri2 = affine.rotate(tri1, 60, (0, 0))

    unit.elements = gpd.GeoDataFrame(
        data = {"element_id": list("abcdef")}, 
        crs = unit.crs,
        geometry = gpd.GeoSeries([hex, sq1, sq2, sq3, tri1, tri2]))
    unit.regularised_tile = copy.deepcopy(unit.tile)
    unit.regularised_tile.geometry = gpd.GeoSeries([
        unit.elements.geometry.buffer(1e-3, join_style = 2) \
            .unary_union.buffer(-1e-3)])
