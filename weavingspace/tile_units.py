#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
import itertools
from enum import Enum

import geopandas as gpd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine


class TileShape(Enum):
    RECTANGLE = "rectangle"
    HEXAGON = "hexagon"
    TRIANGLE = "triangle"
    TRIHEX = "tri-hex"
    TRIDIAMOND = "tri-diamond"


@dataclass
class TileUnit:
    elements:gpd.GeoDataFrame = None
    tile:gpd.GeoDataFrame = None
    tile_shape:TileShape = TileShape.RECTANGLE
    spacing:float = 1000.


    def __init__(self, shape:TileShape = TileShape.RECTANGLE, **kwargs) -> None:
        self.tile_shape = shape
        for k, v in kwargs.items():
            self.__dict__[k] = v
        self.tile = self.get_base_tile()
        self.elements = self.tile.overlay(
            self.make_elements())
        if self.tile_shape == TileShape.TRIANGLE:
            self._modify_elements()
            self._modify_tile()
    
    
    def get_base_tile(self) -> gpd.GeoDataFrame: 
        n = (4 
             if self.tile_shape == TileShape.RECTANGLE
             else (6 
                   if self.tile_shape == TileShape.HEXAGON
                   else 3))
        R = self.spacing / np.cos(np.radians(180 / n)) / 2
        a0 = -90 + 180 / n
        angles = [a0 + a for a in range(0, 360, 360 // n)]
        corners = [(R * np.cos(np.radians(a)), 
                    R * np.sin(np.radians(a))) for a in angles]
        tile = geom.Polygon(corners)
        return gpd.GeoDataFrame(geometry = [tile], crs = self.crs) 
    
    
    # TO CONSIDER: there is probably a refactoring to make this more 
    # elegant - at present the TileGrid class has a similar method...
    # they both operate in the same way, one on a single geometry GeoSeries
    # the other on a multi-element one
    def _modify_elements(self) -> None:
        elements = self.elements.geometry
        ids = list(self.elements.element_id)
        if self.to_hex:
            # rotate to point down
            elements = elements.rotate(180, origin = (0, 0))
        elements = elements.translate(0, -elements.total_bounds[1])
        twins = [elements.rotate(a, origin = (0, 0)) 
                 for a in range(0, 360, (60 if self.to_hex else 180))]
        twins = itertools.chain(*twins)
        self.elements = gpd.GeoDataFrame(
            data = {"element_id": ids * (6 if self.to_hex else 2)},
            geometry = gpd.GeoSeries(twins), crs = self.elements.crs)
        self.elements = self.elements.dissolve(
                                    by = "element_id", as_index = False)
        self.tile.geometry = self._modify_tile()
        self.tile_shape = (TileShape.TRIHEX 
                           if self.to_hex 
                           else TileShape.TRIDIAMOND)
        return None


    # make up a hexagon from the supplied tile because it is a triangle
    # and also change the grid_type to a hexagon
    # assume the triangle is point up with centroid at (0, 0)
    def _modify_tile(self) -> None:
        tile = self.tile.geometry[0]
        if self.to_hex:  # hexagon
            # rotate so point is down
            tile = affine.rotate(tile, 180)
        # translate to sit on x-axis
        tile = affine.translate(tile, 0, -tile.bounds[1])
        # make rotated copies
        twins = [affine.rotate(tile, a, origin = (0, 0)).buffer(0.1) 
                for a in range(0, 360, (60 if self.to_hex else 180))]
        merged_tile = gpd.GeoSeries(twins).unary_union.buffer(-0.1)
        return gpd.GeoSeries([merged_tile])
    

    def make_elements(self) -> gpd.GeoDataFrame:
        return gpd.GeoDataFrame(
            data = {"element_id": list("a")}, crs = self.crs,
            geometry = self.tile.geometry)
        # square = geom.Polygon([(0, 0), (d * 2, 0), (d * 2, d * 2), (0, d * 2)])
        # squares = [affine.rotate(square, a, origin = (0, 0))
        #         for a in range(45, 360, 90)]
        # return gpd.GeoDataFrame(
        #     data = {"element_id": list("abcd")},
        #     geometry = gpd.GeoSeries(squares), crs = crs)
    
    
    def normalize_elements(self) -> None:
        if self.tile_shape in (TileShape.RECTANGLE, ):
            pts = self.tile.geometry[0].exterior.coords
            vectors = [(pts[1][0] - pts[0][0], pts[1][1] - pts[0][1]), 
                       (pts[3][0] - pts[0][0], pts[3][1] - pts[0][1])]
            return vectors
        elif self.tile_shape in (TileShape.HEXAGON, TileShape.TRIHEX):
            pass
        else:
            pass  # TileShape.TRIDIAMOND
        return None
