#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
import copy
import string

import geopandas as gpd
import pandas as pd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine
import shapely.wkt as wkt

import tiling_utils
from tiling_utils import TileShape

import tiling_geometries


@dataclass
class Tileable:
    fudge_factor:float = 1e-3
    spacing:float = 1000.
    crs:int = 3857
    tile_shape:TileShape = TileShape.RECTANGLE
    tile:gpd.GeoDataFrame = None
    vectors:list[tuple[float]] = None
    elements:gpd.GeoDataFrame = None
    regularised_tile:gpd.GeoDataFrame = None
    margin:float = 0.
        
        
    def get_vectors(self, return_values:bool = True) -> list[tuple[float]]:
        """
        Returns symmetry translation vectors as floating point pairs. Derived from the size and shape of the tile attribute. These are not the minimal translation vectors, but the 'face to face' vectors of the tile, such that a hexagonal tile will have 3 vectors, not the minimal parallelogram pair. Also supplies the inverse vectors.

        Returns:
            list[tuple[float]]: _description_
        """
        bb = self.tile.geometry[0].bounds
        w, h = bb[2] - bb[0], bb[3] - bb[1]
        vec_dict = {}
        if self.tile_shape in (TileShape.RECTANGLE, ):
            vec_dict[(1, 0)] = (w, 0)
            vec_dict[(0, 1)] = (0, h)
            vec_dict[(-1, 0)] = (-w, 0)
            vec_dict[(0, -1)] = (0, -h)
        elif self.tile_shape in (TileShape.HEXAGON, ):
            # hex grid coordinates associated with each of the vectors
            i = [0, 1, 1, 0, -1, -1]
            j = [1, 0, -1, -1, 0, 1]
            k = [-1, -1, 0, 1, 1, 0]
            angles = [np.pi * 2 * i / 12 for i in range(1, 12, 2)]
            vecs = [(h * np.cos(a), h * np.sin(a)) for a in angles]
            vec_dict = {(i, j, k): v for i, j, k, v in zip(i, j, k, vecs)}
        else: # DIAMOND
            vec_dict[(1, 0)] = (w / 2, h / 2)
            vec_dict[(0, 1)] = (-w / 2, h / 2)
            vec_dict[(-1, 0)] = (-w / 2, -h / 2)
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
                    u1 = (f1.buffer(self.fudge_factor, join_style = 2) | 
                          f2.buffer(self.fudge_factor, join_style = 2))
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
            self.regularised_tile.geometry.buffer(
                self.fudge_factor, join_style = 2)
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
            self.regularised_tile.geometry.buffer(
                -self.fudge_factor, join_style = 2)
        self.regularised_tile.geometry = \
            self.regularised_tile.geometry.explode(index_parts = False,
                                                   ignore_index = True)[:1]
        return None
    
    
    def get_local_patch(self, r:int = 1, 
                        include_0:bool = False) -> gpd.GeoDataFrame:
        ids = []
        tiles = []
        vecs = ({(0, 0, 0): (0, 0)} 
                if self.tile_shape in (TileShape.HEXAGON, )
                else {(0, 0): (0, 0)})
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
            vecs.pop((0, 0, 0) if self.tile_shape in (TileShape.HEXAGON, )
                     else (0, 0))
        for v in vecs.values():
            ids.extend(self.elements.element_id)
            tiles.extend(self.elements.geometry.apply(
                affine.translate, xoff = v[0], yoff = v[1]))
        return gpd.GeoDataFrame(
            data = {"element_id": ids}, crs = self.crs,
            geometry = gpd.GeoSeries(tiles)
        )
    
    
    def fit_elements_to_tile(self, centre_tile:int = 0) -> None:
        dxy = self.elements.geometry[centre_tile].centroid
        self.elements.geometry = self.elements.translate(-dxy.x, -dxy.y)
        patch = (self.get_local_patch(r = 2, include_0 = True)
                 if self.tile_shape in (TileShape.RECTANGLE, )
                 else self.get_local_patch(r = 1, include_0 = True))
        self.elements = patch.clip(self.tile)
        self.elements.geometry = self.elements.geometry \
            .buffer(-self.fudge_factor) \
            .buffer(self.fudge_factor, join_style = 2)
        self.elements = self.elements[self.elements.geometry.area > 0]
        self.regularised_tile = copy.deepcopy(self.tile)
        return
    
        
    def plot(self, ax = None, show_tile:bool = True, show_reg_tile:bool = True, 
             show_ids = True, show_vectors:bool = False, r:int = 0,
             tile_edgecolor:str = "k", reg_tile_edgcolor:str = "r", 
             facecolor:str = "#00000000", cmap:list[str] = None, 
             figsize:tuple[float] = (8, 8), **kwargs) -> None:
        w = self.tile.geometry[0].bounds[2] - self.tile.geometry[0].bounds[0] 
        n_cols = len(set(self.elements.element_id))
        if cmap is None:
            cm = "Dark2" if n_cols <= 8 else "Paired"
        else:
            cm = cmap
        if ax is None:
            ax = self.elements.plot(column = "element_id", cmap = cm, 
                                    figsize = figsize, **kwargs)
        else:
            self.elements.plot(ax = ax, column = "element_id", cmap = cm, 
                               figsize = figsize, **kwargs)
        if show_ids:
            for id, tile in zip(self.elements.element_id, 
                                self.elements.geometry):
                ax.annotate(id, (tile.centroid.x, tile.centroid.y),
                            ha = "center", va = "center",
                            bbox = {"lw": 0, "fc": "#ffffff40"})
        if r > 0:
            self.get_local_patch(r = r).plot(ax = ax, column = "element_id",
                                             alpha = 0.25, cmap = cm, **kwargs)
        if show_tile:
            self.tile.plot(ax = ax, edgecolor = tile_edgecolor, lw = 0.5,
                           facecolor = facecolor, **kwargs) 
        if show_vectors:
            for v in self.vectors[:len(self.vectors) // 2]:
                ax.arrow(0, 0, v[0], v[1], color = "k", width = w * 0.002,
                         head_width = w * 0.05, length_includes_head = True)
        if show_reg_tile:
            self.regularised_tile.plot(ax = ax, edgecolor = reg_tile_edgcolor,
                                       facecolor = facecolor, linewidth = 2, 
                                       **kwargs)
        return ax
    
    

@dataclass
class TileUnit(Tileable):
    tiling_type:str = None
    dissection_offset:int = 1
    code:str = "3.3.4.3.4"
        
    def __init__(self, **kwargs) -> None:
        for k, v in kwargs.items():
            self.__dict__[k] = v
        self.setup_tile_unit()
        self.vectors = self.get_vectors()
        if self.tile_shape == TileShape.TRIANGLE:
            self._modify_tile()
            self._modify_elements()
        if self.regularised_tile is None: 
            self.regularised_tile = copy.deepcopy(self.tile)
            self.regularise_elements()
        if self.margin > 0:
            self.regularised_tile = self.regularised_tile.scale(
                xfact = 1 - self.margin, yfact = 1 - self.margin)
            self.elements.geometry = self.elements.geometry.scale(
                xfact = 1 - self.margin, yfact = 1 - self.margin,
                origin = self.regularised_tile.geometry[0].centroid)


    def setup_tile_unit(self) -> None:
        if self.tiling_type == "cairo":
            tiling_geometries.setup_cairo(self)
        elif self.tiling_type == "hex-dissection":
            tiling_geometries.setup_hex_dissection(self)
        elif self.tiling_type == "laves":
            tiling_geometries.setup_laves(self)
        elif self.tiling_type == "archimedean":
            tiling_geometries.setup_archimedean(self)
        else:
            tiling_geometries.setup_none_tile(self)
        return

    
    # change the element in a triangle tile to either a hex or
    # a diamond. Assumes the triangle is point up with centroid
    # at (0, 0) 
    def _modify_elements(self) -> None:
        elements = self.elements.geometry
        ids = list(self.elements.element_id)
        new_ids = list(string.ascii_letters[:(len(ids) * 2)])
        elements = elements.translate(0, -elements.total_bounds[1])
        twins = [affine.rotate(element, a, origin = (0, 0)) 
                 for element in elements
                 for a in range(0, 360, 180)]
        self.elements = gpd.GeoDataFrame(
            data = {"element_id": new_ids},
            geometry = gpd.GeoSeries(twins), crs = self.elements.crs)
        self.tile.geometry = self._modify_tile()
        return None


    # make up a diamond from the supplied tile because it is a triangle
    # and also change the grid_type to a diamond
    # assume the triangle is point up with centroid at (0, 0)
    def _modify_tile(self) -> None:
        tile = self.tile.geometry[0]
        # translate to sit on x-axis
        tile = affine.translate(tile, 0, -tile.bounds[1])
        # make rotated copies
        twins = [affine.rotate(tile, a, origin = (0, 0)).buffer(
                                self.fudge_factor, join_style = 2)
                 for a in range(0, 360, 180)]
        merged_tile = \
            gpd.GeoSeries(twins).unary_union.buffer(
                -self.fudge_factor, join_style = 2)
        self.tile_shape = TileShape.DIAMOND
        return gpd.GeoSeries([merged_tile])
                

    def plot_legend(self, ax, vars:dict[str:str], pals:dict[str:str], 
                    data:dict[str:list], zoom:float = 1., 
                    map_rotation:float = 0, rotate_text:bool = False, **kwargs):
        ax.set_axis_off()
        tiles, ids, vals = [], [], []
        for id, t in zip(self.elements.element_id, 
                         self.elements.geometry):
            data_vals = data[id]
            data_vals.sort()
            n = len(data_vals)
            tiles.extend(tiling_utils.get_insets(t, n))
            ids.extend([id] * n)
            vals.extend(data_vals)

        gdf = gpd.GeoDataFrame(
            data = {"val": vals, "id": ids}, crs = self.crs, 
            geometry = gpd.GeoSeries(tiles))
        gdf.geometry = gdf.geometry.rotate(map_rotation, origin = (0, 0)) 
        
        bb = [x / zoom for x in gdf.geometry.total_bounds]
        ax.set_xlim(bb[0], bb[2])
        ax.set_ylim(bb[1], bb[3])
        ax.axhspan(bb[1], bb[3], fc = "lightgrey", lw = 0)

        self.get_local_patch(r = 2) \
            .geometry.rotate(map_rotation, origin = (0, 0)).plot(
                ax = ax, fc = "w", ec = "grey", lw = 0.5)

        groups = gdf.groupby("id")
        for id in pd.Series.unique(gdf.id):
            item = groups.get_group(id)
            item.plot(ax = ax, column = "val", 
                      cmap = pals[id], lw = 0, **kwargs)
        
        rotated_unit = copy.deepcopy(self)
        rotated_unit.elements.geometry = \
            rotated_unit.elements.geometry.rotate(map_rotation, origin = (0, 0))
        for id, ele in zip(rotated_unit.elements.element_id, 
                           rotated_unit.elements.geometry):
            c = ele.centroid
            rot = (0 if not rotate_text or c.x == 0
                else (np.degrees(np.arctan2(c.y, c.x)) + 90) % 180 - 90)
            ax.annotate(vars[id], xy = (c.x, c.y),
                        ha = ("left" if c.x >= 0 else "right"), va = "center",
                        rotation = rot, rotation_mode = "anchor",
                        bbox = {"lw": 0, "fc": "#ffffff40"})
        return None
    