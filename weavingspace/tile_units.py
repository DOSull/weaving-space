#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
import itertools
from enum import Enum
import copy

import geopandas as gpd
import pandas as pd
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
class Tileable:
    fudge_factor:float = 1e-3
    elements:gpd.GeoDataFrame = None
    tile:gpd.GeoDataFrame = None
    regularised_tile:gpd.GeoDataFrame = None
    vectors:list[tuple[float]] = None
    spacing:float = 1000.
    tile_shape:TileShape = TileShape.RECTANGLE
    crs:int = 3857
        
    def get_vectors(self) -> list[tuple[float]]:
        """
        Returns symmetry translation vectors as floating point pairs. Derived from the size and shape of the tile attribute. These are not the minimal translation vectors, but the 'face to face' vectors of the tile, such that a hexagonal tile will have 3 vectors, not the minimal parallelogram pair. Also supplies the inverse vectors.

        Returns:
            list[tuple[float]]: _description_
        """
        bb = self.tile.geometry[0].bounds
        w, h = bb[2] - bb[0], bb[3] - bb[1]
        if self.tile_shape in (TileShape.RECTANGLE, ):
            return [(dx, dy) 
                    for dx in (-w, 0, w) 
                    for dy in (-h, 0, h)
                    if (dx == 0 or dy == 0) and dx != dy]
        elif self.tile_shape in (TileShape.HEXAGON, TileShape.TRIHEX):
            angles = [np.pi * 2 * i / 12 for i in range(1, 12, 2)]
            return [(h * np.cos(a), h * np.sin(a)) for a in angles]
        else: # TRIDIAMOND
            return None
        
    
    def merge_fragments(
            self, fragments:list[geom.Polygon]) -> list[geom.Polygon]:
        """
        Merges a set of polygons based on testing if they touch when subjected to the translation vectors provided by the get_vectors() method. Called by regularise_elements() method to combine elements in a tile that may be fragmented as supplied but will combine when tiled into single elements. This step makes for more efficient implementation of the tiling of map regions.

        Args:
            fragments (list[geom.Polygon]): A set of polygons to merge.

        Returns:
            list[geom.Polygon]: A minimal list of merged polygons.
        """
        if len(fragments) == 1:
            return fragments
        changes_made = True
        tile = self.tile.geometry[0]
        reg_tile = self.regularised_tile.geometry[0]
        while changes_made:
            changes_made = False
            for v in self.vectors:
                next_frags = []
                t_fragments = [affine.translate(f, v[0], v[1]) 
                               for f in fragments]
                matches = set()
                for i, f1 in enumerate(fragments):
                    for j, f2, in enumerate(t_fragments):
                        if i != j and f1.distance(f2) < 1e-3:
                            matches.add((i, j))
                fragments_to_remove = set()
                for i, j in matches:
                    f1 = fragments[i]
                    f2 = t_fragments[j]
                    u1 = (f1.buffer(self.fudge_factor) | 
                          f2.buffer(self.fudge_factor))
                    u2 = affine.translate(u1, -v[0], -v[1])
                    if tile.intersection(u1).area > tile.intersection(u2).area:
                        next_frags.append(u1)
                        reg_tile = reg_tile | u1
                        reg_tile = reg_tile - u2
                    else:
                        next_frags.append(u2)
                        reg_tile = reg_tile | u2
                        reg_tile = reg_tile - u1
                    changes_made = True
                    fragments_to_remove.add(i)
                    fragments_to_remove.add(j)
                fragments = [f for i, f in enumerate(fragments) 
                            if not (i in fragments_to_remove)]
                fragments = next_frags + fragments
        self.regularised_tile.geometry[0] = reg_tile
        return fragments


    def regularise_elements(self) -> None:
        """
        Combines separate elements that share an element_id value into single elements, if they would end up touching when tiled. Also adjusts the regularised_tile attribute according.

        Returns:
            None: 
        """
        self.vectors = self.get_vectors()
        elements = []
        element_ids = []
        self.regularised_tile.geometry = \
            self.regularised_tile.geometry.buffer(self.fudge_factor)
        # This preserves order while finding uniques, unlike list(set()).
        # Reordering ids might cause confusion when colour palettes
        # are not assigned explicitly to each id, but in the order
        # encountered in the element_id Series of the GeoDataFrame.
        ids = list(pd.Series.unique(self.elements.element_id))
        for id in ids:
            fragment_set = list(
                self.elements[self.elements.element_id == id].geometry)
            merge_result = self.merge_fragments(fragment_set)
            elements.extend(merge_result)
            element_ids.extend([id] * len(merge_result))
        new_elements = gpd.GeoDataFrame(
            data = {"element_id": element_ids}, crs = self.crs,
            geometry = gpd.GeoSeries(elements)
        )
        self.elements = new_elements
        self.regularised_tile.geometry = \
            self.regularised_tile.geometry.buffer(-self.fudge_factor)
        return None



@dataclass
class TileUnit(Tileable):

    def __init__(self, **kwargs) -> None:
        unit = self.get_tile_unit(**kwargs)
        self.elements = unit["elements"]
        self.tile = unit["tile"]
        self.regularised_tile = unit["regularised_tile"]
        for k, v in kwargs.items():
            self.__dict__[k] = v
        if self.tile_shape == TileShape.TRIANGLE:
            self._modify_elements()
            self._modify_tile()
        if self.regularised_tile is None:
            self.regularised_tile = copy.copy(self.tile)
            self.regularise_elements()
    
    
    def get_base_tile(self, shape, spacing, crs) -> gpd.GeoDataFrame: 
        n = (4 
             if shape == TileShape.RECTANGLE
             else (6 
                   if shape == TileShape.HEXAGON
                   else 3))
        R = spacing / np.cos(np.radians(180 / n)) / 2
        a0 = -90 + 180 / n
        angles = [a0 + a for a in range(0, 360, 360 // n)]
        corners = [(R * np.cos(np.radians(a)), 
                    R * np.sin(np.radians(a))) for a in angles]
        tile = geom.Polygon(corners)
        return gpd.GeoDataFrame(geometry = [tile], crs = crs) 
    
    # change the element in a triangle tile to either a hex or
    # a diamond. Assumes the triangle is point up with centroid
    # at (0, 0) 
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
        self.elements = self.elements.dissolve(by = "element_id", 
                                               as_index = False)
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
    
    
    def get_tile_unit(self, tiling_type:str = "cairo", spacing:float = 10000,
                      tile_shape:TileShape = TileShape.RECTANGLE, to_hex = False, crs:int = 3857) -> dict:
        return self.get_cairo(spacing, tile_shape, crs)


    def get_cairo(self, spacing, tile_shape, crs):
        d = spacing
        x = d / 2 * (np.sqrt(3) - 1) / 2 / np.sqrt(3)

        p1 = geom.Polygon([(-d / 2 + x, 0),
                           (-d / 4, -d / 4), 
                           (0, -x), 
                           (0, x),
                           (-d / 4, d / 4)])
        p1 = affine.translate(p1, x, x) 
        p2 = affine.rotate(p1, 90, origin = p1.exterior.coords[1])
        p3 = affine.rotate(p1, 180, origin = p1.exterior.coords[1])
        p4 = affine.rotate(p1, 270, origin = p1.exterior.coords[1])
        p5 = affine.rotate(p4, 90, origin = p4.exterior.coords[4]) 
        p6 = affine.rotate(p4, 270, origin = p4.exterior.coords[4])
        p7 = affine.rotate(p6, 270, origin = p6.exterior.coords[1])
        p8 = affine.rotate(p7, 180, origin = p7.exterior.coords[4])

        elements = gpd.GeoDataFrame(
            data = {"element_id": list("abcdacbd")}, crs = crs,
            geometry = gpd.GeoSeries([p1, p2, p3, p4, p5, p6, p7, p8]))
        tile = self.get_base_tile(tile_shape, spacing, crs)
        reg_tile = copy.copy(tile)
        reg_tile.geometry = gpd.GeoSeries([elements.geometry.unary_union])
        
        return {"elements": elements, 
                "tile": tile, 
                "regularised_tile": reg_tile}