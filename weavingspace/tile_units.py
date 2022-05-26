#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
import itertools
from enum import Enum
import copy
import string

import geopandas as gpd
import pandas as pd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine

import tiling_utils


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
        
    def get_vectors(self, return_values:bool = True) -> list[tuple[float]]:
        """
        Returns symmetry translation vectors as floating point pairs. Derived from the size and shape of the tile attribute. These are not the minimal translation vectors, but the 'face to face' vectors of the tile, such that a hexagonal tile will have 3 vectors, not the minimal parallelogram pair. Also supplies the inverse vectors.

        Returns:
            list[tuple[float]]: _description_
        """
        bb = self.tile.geometry[0].bounds
        w, h = bb[2] - bb[0], bb[3] - bb[1]
        vec_dict = ({(0, 0): (0, 0)} 
                    if self.tile_shape == TileShape.RECTANGLE
                    else {(0, 0, 0): (0, 0)})
        if self.tile_shape in (TileShape.RECTANGLE, ):
            vec_dict[(1, 0)] = (w, 0)
            vec_dict[(0, 1)] = (0, h)
            vec_dict[(-1, 0)] = (-w, 0)
            vec_dict[(0, -1)] = (0, -h)
        elif self.tile_shape in (TileShape.HEXAGON, TileShape.TRIHEX):
            # hex grid coordinates associated with each of the vectors
            i = [0, 1, 1, 0, -1, -1]
            j = [1, 0, -1, -1, 0, 1]
            k = [-1, -1, 0, 1, 1, 0]
            angles = [np.pi * 2 * i / 12 for i in range(1, 12, 2)]
            vecs = [(h * np.cos(a), h * np.sin(a)) for a in angles]
            vec_dict = {(i, j, k): v for i, j, k, v in zip(i, j, k, vecs)}
        else: # TRIDIAMOND
            vec_dict[(1, 0)] = (w / 2, h / 2)
            vec_dict[(0, 1)] = (-w / 2, h / 2)
            vec_dict[(-1, 0)] = (-w / 2. -h / 2)
            vec_dict[(0, -1)] = (w / 2, -h / 2)
        return (list(vec_dict.values())
                if return_values
                else vec_dict)
        
    
    def merge_fragments(
            self, fragments:list[geom.Polygon]) -> list[geom.Polygon]:
        """
        Merges a set of polygons based on testing if they touch when subjected to the translation vectors provided by the get_vectors() method. Called by regularise_elements() method to combine elements in a tile that may be fragmented as supplied but will combine when tiled into single elements. This step makes for more efficient implementation of the tiling of map regions.

        Args:
            fragments (list[geom.Polygon]): A set of polygons to merge.

        Returns:
            list[geom.Polygon]: A minimal list of merged polygons.
        """
        tile = self.tile.geometry[0]
        reg_tile = self.regularised_tile.geometry[0]
        if len(fragments) == 1:
            return fragments
        changes_made = True
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
    
    
    def get_local_patch(self, r:int = 1, 
                        include_0:bool = False) -> gpd.GeoDataFrame:
        ids = []
        tiles = []
        vecs = ({(0, 0): (0, 0)} 
                if self.tile_shape == TileShape.RECTANGLE
                else {(0, 0, 0): (0, 0)})
        last_vecs = copy.deepcopy(vecs)
        vectors = self.get_vectors(return_values = False)
        for i in range(r):
            new_vecs = {}
            for k1, v1 in last_vecs.items():
                for k2, v2 in vectors.items():
                    new_key = tuple([k1[i] + k2[i] for i in range(len(k1))])
                    new_val = (v1[0] + v2[0], v1[1] + v2[1])
                    if not new_val in vecs: 
                        new_vecs[new_key] = new_val
            vecs = vecs | new_vecs
            last_vecs = new_vecs
        if not include_0:
            vecs.pop((0, 0) if self.tile_shape == TileShape.RECTANGLE
                     else (0, 0, 0))
        for v in vecs.values():
            ids.extend(self.elements.element_id)
            tiles.extend(self.elements.geometry.apply(
                affine.translate, xoff = v[0], yoff = v[1]))
        return gpd.GeoDataFrame(
            data = {"element_id": ids}, crs = self.crs,
            geometry = gpd.GeoSeries(tiles)
        )
        

@dataclass
class TileUnit(Tileable):
    tiling_type:str = None
    dissection_offset:int = 1
    laves_code:int = "3.3.4.3.4"
    to_hex:bool = False
        
    def __init__(self, **kwargs) -> None:
        for k, v in kwargs.items():
            self.__dict__[k] = v
        self.setup_tile_unit()
        if self.tile_shape == TileShape.TRIANGLE:
            self._modify_elements()
            self._modify_tile()
        if self.regularised_tile is None: 
            self.regularised_tile = copy.deepcopy(self.tile)
            self.regularise_elements()
    
    
    def get_base_tile(self) -> gpd.GeoDataFrame: 
        tile = tiling_utils.get_regular_polygon(
            self.spacing, n = (4 
                          if self.tile_shape in (TileShape.RECTANGLE, )
                          else (6 
                                if self.tile_shape in (TileShape.HEXAGON, 
                                                       TileShape.TRIHEX)
                                else 3)))
        return gpd.GeoDataFrame(geometry = [tile], crs = self.crs) 
 
    
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
    
    
    def setup_tile_unit(self) -> None:
        if self.tiling_type == "cairo":
            self.setup_cairo()
        elif self.tiling_type == "hex-dissection":
            self.setup_hex_dissection()
        elif self.tiling_type == "laves":
            self.setup_laves_tiles()
        else:
            self.setup_none_tile()
        return
            
    
    def setup_none_tile(self) -> None: 
        self.tile = self.get_base_tile()
        self.elements = gpd.GeoDataFrame(
            data = {"element_id": ["a"]}, crs = self.crs,
            geometry = copy.deepcopy(self.tile.geometry))


    def setup_cairo(self) -> None:
        d = self.spacing
        # x = d / 2 * (np.sqrt(3) - 1) / 2 / np.sqrt(3)
        
        # p1 = geom.Polygon([(-d / 2 + x, 0),
        #                    (-d / 4, -d / 4), 
        #                    (0, -x), 
        #                    (0, x),
        #                    (-d / 4, d / 4)])
        # p1 = affine.translate(p1, x, x) 
        # p2 = affine.rotate(p1, 90, origin = p1.exterior.coords[1])
        # p3 = affine.rotate(p1, 180, origin = p1.exterior.coords[1])
        # p4 = affine.rotate(p1, 270, origin = p1.exterior.coords[1])
        # p5 = affine.rotate(p4, 90, origin = p4.exterior.coords[4]) 
        # p6 = affine.rotate(p4, 270, origin = p4.exterior.coords[4])
        # p7 = affine.rotate(p6, 270, origin = p6.exterior.coords[1])
        # p8 = affine.rotate(p7, 180, origin = p7.exterior.coords[4])

        # self.elements = gpd.GeoDataFrame(
        #     data = {"element_id": list("abcdacbd")}, crs = self.crs,
        #     geometry = gpd.GeoSeries([p1, p2, p3, p4, p5, p6, p7, p8]))
        
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

        self.elements = gpd.GeoDataFrame(
            data = {"element_id": list("abcd")}, crs = self.crs,
            geometry = gpd.GeoSeries([p1, p2, p3, p4]))
        self.tile_shape = TileShape.RECTANGLE
        self.tile = self.get_base_tile()
        self.regularised_tile = copy.deepcopy(self.tile)
        self.regularised_tile.geometry = gpd.GeoSeries(
            [self.elements.geometry.buffer(1e-3).unary_union.buffer(-1e-3)])
        
            
        
    def setup_hex_dissection(self):
        self.tile_shape = TileShape.HEXAGON
        hex = tiling_utils.get_regular_polygon(self.spacing, 6)
        v = list(hex.exterior.coords)
        m = [((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
                     for p1, p2 in zip(v[:-1], v[1:])]
        p = list(itertools.chain(*zip(v[:-1], m)))
        if self.n == 2:
            s = (geom.Polygon([p[0], p[2], p[4], p[6], (0, 0)])
                 if self.dissection_offset == 0
                 else geom.Polygon(p[1:8] + [(0, 0)]))
            slices = [affine.rotate(s, a, origin = (0, 0)) 
                      for a in range(0, 360, 180)]
        elif self.n == 3:
            s = (geom.Polygon([p[0], p[2], p[4], (0, 0)])
                 if self.dissection_offset == 0
                 else geom.Polygon([p[1], p[2], p[4], p[5], (0, 0)]))
            slices = [affine.rotate(s, a, origin = (0, 0)) 
                      for a in range(0, 360, 120)]
        elif self.n == 4:
            s1 = (geom.Polygon([p[0], p[2], p[3], (0, 0)])
                 if self.dissection_offset == 0
                 else geom.Polygon([p[1], p[2], p[4], (0, 0)]))
            s2 = (geom.Polygon([p[3], p[4], p[6], (0, 0)])
                 if self.dissection_offset == 0
                 else geom.Polygon([p[4], p[6], p[7], (0, 0)]))
            slices = [s1, s2, 
                      affine.rotate(s1, 180, (0, 0)),
                      affine.rotate(s2, 180, (0, 0))]
        elif self.n == 6:
            s = (geom.Polygon([p[0], p[2], (0, 0)])
                 if self.dissection_offset == 0
                 else geom.Polygon([p[1], p[2], p[3], (0, 0)]))
            slices = [affine.rotate(s, a, origin = (0, 0)) 
                      for a in range(0, 360, 60)]
        elif self.n == 12:
            s1 = geom.Polygon([p[0], p[1], (0, 0)])
            s2 = geom.Polygon([p[1], p[2], (0, 0)])
            slices1 = [affine.rotate(s1, a, origin = (0, 0)) 
                       for a in range(0, 360, 60)]
            slices2 = [affine.rotate(s2, a, origin = (0, 0)) 
                       for a in range(0, 360, 60)]
            slices = itertools.chain(slices1, slices2)
        
        self.elements = gpd.GeoDataFrame(
            data = {"element_id": list(string.ascii_lowercase)[:self.n]}, 
            crs = self.crs,
            geometry = gpd.GeoSeries(slices))
        self.tile_shape = TileShape.HEXAGON
        self.tile = self.get_base_tile()


    def setup_laves_tiles(self) -> None:
        if self.laves_code == "3.3.3.3.3.3":
            self.tile_shape = TileShape.HEXAGON
            self.setup_none_tile()
            return
        if self.laves_code == "3.3.3.3.6":
            self.tile_shape = TileShape.HEXAGON
            self.setup_laves_33336()
            return
        elif self.laves_code == "3.3.3.4.4":
            self.tile_shape = TileShape.HEXAGON
            self.n = 2
            self.dissection_offset = 1
            self.setup_hex_dissection()
            return
        elif self.laves_code == "3.3.4.3.4":
            self.tile_shape = TileShape.RECTANGLE
            self.setup_cairo()
            return
        elif self.laves_code == "3.4.6.4":
            self.tile_shape = TileShape.HEXAGON
            self.n = 6
            self.dissection_offset = 1
            self.setup_hex_dissection()
            return
        elif self.laves_code == "3.6.3.6":
            self.tile_shape = TileShape.HEXAGON
            self.n = 3
            self.dissection_offset = 0
            self.setup_hex_dissection()
            return
        elif self.laves_code == "3.12.12":
            print(f"The code [{self.laves_code}] is unsupported.")
        elif self.laves_code == "4.4.4":
            self.tile_shape = TileShape.RECTANGLE
            self.setup_none_tile()
            return
        elif self.laves_code == "4.6.12":
            self.tile_shape = TileShape.HEXAGON
            self.n = 12
            self.setup_hex_dissection()
            return
        elif self.laves_code == "4.8.8":
            print(f"The code [{self.laves_code}] is unsupported.")
        elif self.laves_code == "6.6.6":
            self.tile_shape = TileShape.TRIANGLE
            self.setup_none_tile()
        else:
            print(f"[{self.laves_code}] is not a valid Laves code.")

        self.tiling_type = None
        self.setup_none_tile()
        return
    
    
    def setup_laves_33336(self) -> None:
        offset_a = np.degrees(np.arctan(1 / 3 / np.sqrt(3)))
        sf = 1 / np.sqrt(7)

        tile = tiling_utils.get_regular_polygon(self.spacing, 6)

        hex = affine.scale(tile, sf, sf)
        hex = affine.translate(hex, 0, hex.bounds[3] - hex.bounds[1])
        hex_p = [p for p in hex.exterior.coords]
        petal = geom.Polygon([(0, 0)] + hex_p[1:5])
        petals = [affine.rotate(petal, a + offset_a, origin = (0, 0))
                     for a in range(30, 360, 60)]

        self.elements = gpd.GeoDataFrame(
            data = {"element_id": list("abcdef")}, 
            crs = self.crs,
            geometry = gpd.GeoSeries(petals))
        self.tile_shape = TileShape.HEXAGON
        self.tile = self.get_base_tile()
        self.regularised_tile = gpd.GeoSeries(
            [self.elements.geometry.buffer(1e-3).unary_union.buffer(-1e-3)])
            
