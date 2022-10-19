#!/usr/bin/env python
# coding: utf-8

"""Implements `weavingspace.tileable.TileShape` and 
`weavingspace.tileable.Tileable` the base classes for 
`weavingspace.tile_unit.TileUnit` and `weavingspace.weave_unit.WeaveUnit`.

`Tileable` should not be called directly, but is instead accessed from the 
`weavingspace.tile_unit.TileUnit` or `weavingspace.weave_unit.WeaveUnit` 
constructor. 

Several methods of `weavingspace.tileable.Tileable` are generally useful and 
can be accessed through its subclasses.
"""

from enum import Enum
from typing import Union
from dataclasses import dataclass
import copy

import matplotlib.pyplot as pyplot

import geopandas as gpd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine

import weavingspace.tiling_utils as tiling_utils


class TileShape(Enum):
    """The available base tile shapes.
    
    NOTE: the TRIANGLE type does not persist, but should be converted to a
    DIAMOND or HEXAGON type during `Tileable` construction.
    """
    RECTANGLE = "rectangle"
    HEXAGON = "hexagon"
    TRIANGLE = "triangle"
    DIAMOND = "diamond"


@dataclass
class Tileable:
    """Class to represent a tileable set of geometries.

    Args:
        elements (gpd.GeoDataFrame): the geometries with associated
            element_id attribute encoding their different colouring.
        tile (gpd.GeoDataFrame): the tileable polygon (a rectangle, a
            hexagon, or a diamond).
        spacing (float): the tile spacing -- effectively the 'resolution' of
            the tiling. Defaults to `1000`.
        tile_shape (TileShape): the tile shape. Defaults to 
            `TileShape.RECTANGLE`.
        vectors (list[tuple[float]]): translation vector symmetries of the 
            tiling.
        regularised_tile (gpd.GeoDataFrame): a polygon containing the 
            elements of this `weavingspace.tileable.Tileable` --- most 
            often a union of those polygons.
        crs (int): the coordinate reference system of the tile. Most often 
            an EPSG code, but any valid geopandas CRS specification is
            valid. Defaults to 3857 (i.e. Web Mercator).
        fudge_factor (float): a distance in units of self.crs to be used in
            geometry clean ups (for example this buffer distance is applied
            before unioning polygons.) Defaults to `1e-3`.
        debug (bool, optional): if True prints debug messages. Defaults to
            False.
    """
    elements:gpd.GeoDataFrame = None
    tile:gpd.GeoDataFrame = None
    spacing:float = 1000.
    tile_shape:TileShape = TileShape.RECTANGLE
    vectors:list[tuple[float]] = None
    regularised_tile:gpd.GeoDataFrame = None
    crs:int = 3857
    fudge_factor:float = 1e-3
    rotation:float = 0.0
    debug:bool = False
    
    # Tileable constructor called by subclasses - should not be used directly
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            self.__dict__[k] = v
        if self.debug:
            print(f"""Debugging messages enabled for Tileable (but there aren't
                  any at the moment...)""")
        self._setup_tile_and_elements()
        self.setup_vectors()
        if self.regularised_tile is None:
            self.regularise_elements()
        return
        
        
    def setup_vectors(self) -> None:
        """Sets up the translation vectors of the tiling.
        """
        self.vectors = self.get_vectors()
        return None
    
        
    def get_vectors(self, as_dict:bool = False
        ) -> Union[ dict[tuple[int], tuple[float]],list[tuple[float]] ]:
        """
        Returns symmetry translation vectors as floating point pairs. 
        
        Derived from the size and shape of the tile attribute. These are not the minimal translation vectors, but the 'face to face' vectors of the tile, such that a hexagonal tile will have 3 vectors, not the minimal parallelogram pair. Also supplies the inverse vectors.
        
        Optionally returns the vectors in a dictionary indexed by their 
        coordinates, e.g.
        
            {
                ( 1,  0): ( 100, 0), ( 0,  1): (0,  100),
                (-1,  0): (-100, 0), ( 0, -1): (0, -100)
            }
            
        For a tileable of type `TileShape.HEXAGON`, the indexing tuples 
        have three components. See https://www.redblobgames.com/grids/hexagons/ 
        
        Args:
            as_dict (bool): If `False` returns the vectors only. If
                `True` returns a dictionary of the vectors indexed by tuples
                in the grid coordinate system. Defaults to `False`.
        
        Returns:
            Union[ dict[tuple[int],tuple[float]], list[tuple[float]] ]: 
                either the vectors as a list of float tuples, or a dictionary 
                of those vectors indexed by integer coordinate tuples. 
        """
        t = self.tile.geometry[0]
        pts = [p for p in t.exterior.coords][:-1]
        n_pts = len(pts)
        vec_dict = {}
        if n_pts == 4:
            vecs = [(q[0] - p[0], q[1] - p[1]) 
                    for p, q in zip(pts, pts[1:] + pts[:1])]
            vec_dict[(1, 0)] = vecs[0]  # (w, 0)
            vec_dict[(0, 1)] = vecs[1]  # (0, h)
            vec_dict[(-1, 0)] = vecs[2]  # (-w, 0)
            vec_dict[(0, -1)] = vecs[3]  # (0, -h)
        elif n_pts == 6:
            vecs = [(q[0] - p[0], q[1] - p[1]) 
                    for p, q in zip(pts, pts[2:] + pts[:2])]
            # hex grid coordinates associated with each of the vectors
            i = [0, 1, 1, 0, -1, -1]
            j = [1, 0, -1, -1, 0, 1]
            k = [-1, -1, 0, 1, 1, 0]
            vec_dict = {(i, j, k): v for i, j, k, v in zip(i, j, k, vecs)}
        else: # garbage
            pass
        return (vec_dict
                if as_dict
                else list(vec_dict.values()))
        
        
    # Make up a regularised tile by carefully unioning the elements
    def setup_regularised_tile_from_elements(self) -> None:
        """Sets the regularised tile to a union of the elements.
        """
        self.regularised_tile = copy.deepcopy(self.tile)
        self.regularised_tile.geometry = tiling_utils.safe_union(
            self.elements.geometry, self.spacing * 1e-6)
        # This simplification seems very crude but fixes all kinds of issues...
        self.regularised_tile.geometry[0] = \
            self.regularised_tile.geometry[0].simplify(self.spacing * 1e-6)
        return
        
    
    def merge_fragments(
                self, fragments:list[geom.Polygon]) -> list[geom.Polygon]:
        """
        Merges a set of polygons based on testing if they touch when subjected 
        to the translation vectors provided by `get_vectors()`. 
        
        Called by `regularise_elements()` to combine elements in a tile that 
        may be fragmented as supplied but will combine when tiled into single 
        elements. This step makes for more efficient implementation of the 
        tiling of map regions.

        Args:
            fragments (list[geom.Polygon]): A set of polygons to merge.

        Returns:
            list[geom.Polygon]: A minimal list of merged polygons.
        """
        tile = self.tile.geometry[0]
        reg_tile = self.regularised_tile.geometry[0].buffer(
            self.fudge_factor, resolution = 1, join_style = 2)
        if len(fragments) == 1:
            return fragments
        changes_made = True
        while changes_made:
            changes_made = False
            for v in self.vectors:
                # empty list to collect the new fragments
                # assembled in this iteration
                next_frags = []
                t_fragments = [affine.translate(f, v[0], v[1]) 
                               for f in fragments]
                # build a set of any near matching pairs of 
                # fragments and their translated copies
                matches = set()
                for i, f1 in enumerate(fragments):
                    for j, f2, in enumerate(t_fragments):
                        if i != j and tiling_utils.touch_along_an_edge(f1, f2):
                            matches.add((i, j))
                # determine which of these when unioned has the 
                # larger area in common with the base tile
                fragments_to_remove = set()
                for i, j in matches:
                    f1 = fragments[i]
                    f2 = t_fragments[j]
                    u1 = (f1.buffer(self.fudge_factor, 
                                    resolution = 1, join_style = 2) | 
                          f2.buffer(self.fudge_factor, 
                                    resolution = 1, join_style = 2))
                    u2 = affine.translate(u1, -v[0], -v[1])
                    if tile.intersection(u1).area > tile.intersection(u2).area:
                        u1 = u1.buffer(-self.fudge_factor,
                                       resolution = 1, join_style = 2)
                        next_frags.append(u1)
                        reg_tile = (reg_tile | u1) - u2
                    else:
                        u2 = u2.buffer(-self.fudge_factor,
                                       resolution = 1, join_style = 2)
                        next_frags.append(u2)
                        reg_tile = (reg_tile | u2) - u1
                    changes_made = True
                    fragments_to_remove.add(i)
                    fragments_to_remove.add(j)
                fragments = [f for i, f in enumerate(fragments) 
                            if not (i in fragments_to_remove)]
                fragments = next_frags + fragments
        self.regularised_tile.geometry[0] = \
            reg_tile.buffer(-self.fudge_factor, resolution = 1, join_style = 2)
        return [f for f in fragments if not f.is_empty] # don't return any duds


    def regularise_elements(self) -> None:
        """Combines separate elements that share an element_id value into 
        single elements, if they would end up touching when tiled. 
        
        Also adjusts the `weavingspace.tileable.Tileable.regularised_tile` 
        attribute accordingly.
        """
        self.regularised_tile = copy.deepcopy(self.tile)
        # This preserves order while finding uniques, unlike list(set()).
        # Reordering ids might cause confusion when colour palettes
        # are not assigned explicitly to each id, but in the order
        # encountered in the element_id Series of the GeoDataFrame.
        elements, element_ids = [], []
        ids = list(self.elements.element_id.unique())
        for id in ids:
            fragment_set = list(
                self.elements[self.elements.element_id == id].geometry)
            merge_result = self.merge_fragments(fragment_set)
            elements.extend(merge_result)
            element_ids.extend([id] * len(merge_result))

        self.elements = gpd.GeoDataFrame(
            data = {"element_id": element_ids}, crs = self.crs,
            geometry = gpd.GeoSeries(elements))

        self.regularised_tile.geometry = \
            self.regularised_tile.geometry.explode(index_parts = False,
                                                   ignore_index = True)
        if self.regularised_tile.geometry.shape[0] > 1:
            self.regularised_tile.geometry = \
                tiling_utils.get_largest_polygon(self.regularised_tile.geometry)
        self.regularised_tile.geometry[0] = \
            self.regularised_tile.geometry[0].simplify(self.spacing / 100)
        return None
    
    
    def get_local_patch(
            self, r:int = 1, include_0:bool = False) -> gpd.GeoDataFrame:
        """Returns a GeoDataFrame with translated copies of the Tileable.
        
        The geodataframe takes the same form as the `Tileable.tile` attribute.

        Args:
            r (int, optional): the number of translation vector steps required. 
                Defaults to `1`.
            include_0 (bool, optional): If True includes the Tileable itself at 
                (0, 0). Defaults to `False`.

        Returns:
            gpd.GeoDataFrame: A GeoDataframe of the elements extended
                by a number of translation vectors.
        """
        # a dictionary of all the vectors we need, starting with (0, 0)
        vecs = ({(0, 0, 0): (0, 0)} 
                if self.tile_shape in (TileShape.HEXAGON, )
                else {(0, 0): (0, 0)})
        # a dictionary of the last 'layer' of added vectors
        last_vecs = copy.deepcopy(vecs)
        # get the translation vectors in a dictionary indexed by coordinates
        # we keep track of the sum of vectors using the (integer) coordinates
        # to avoid duplication of moves due to floating point inaccuracies 
        vectors = self.get_vectors(as_dict = True)
        for i in range(r):
            new_vecs = {}
            for k1, v1 in last_vecs.items():
                for k2, v2 in vectors.items():
                    # add the coordinates to make a new key...
                    new_key = tuple([k1[i] + k2[i] for i in range(len(k1))])
                    # ... and the vector components to make a new value
                    new_val = (v1[0] + v2[0], v1[1] + v2[1])
                    # if we haven't reached here before store it
                    if not new_val in vecs: 
                        new_vecs[new_key] = new_val
            # extend the vectors and set the last layer to the set just added
            vecs = vecs | new_vecs
            last_vecs = new_vecs
        if not include_0:  # throw away the identity vector
            vecs.pop((0, 0, 0) if self.tile_shape in (TileShape.HEXAGON, )
                     else (0, 0))
        ids, tiles = [], []
        for v in vecs.values():
            ids.extend(self.elements.element_id)
            tiles.extend(self.elements.geometry.apply(
                affine.translate, xoff = v[0], yoff = v[1]))
        return gpd.GeoDataFrame(
            data = {"element_id": ids}, crs = self.crs,
            geometry = gpd.GeoSeries(tiles)
        )
    
    
    def fit_elements_to_tile(self, centre_element:int = 0) -> None:
        """Fits the tile elements so they sit inside the tile boundary.
        
        If elements in a tile project outside the boundaries of the tile, this
        method will clip them so that they don't. This may result in 
        'fragmented' elements, i.e. elements that would form a single element
        after tiling which are separated into distinct fragments.

        Args:
            centre_element (int, optional): the index position of the central 
                element. Defaults to `0`.
        """
        dxy = self.elements.geometry[centre_element].centroid
        self.elements.geometry = self.elements.translate(-dxy.x, -dxy.y)
        # use r = 2 because rectangular tiles may need diagonal neighbours
        patch = (self.get_local_patch(r = 2, include_0 = True)
                 if self.tile_shape in (TileShape.RECTANGLE, )
                 else self.get_local_patch(r = 1, include_0 = True))
        self.elements = patch.clip(self.tile)
        # repair any weirdness...
        self.elements.geometry = tiling_utils.clean_polygon(
            self.elements.geometry)
        self.elements = self.elements[self.elements.geometry.area > 0]
        self.regularised_tile = copy.deepcopy(self.tile)
        return None
    

    # applicable to both TileUnits and WeaveUnits
    def inset_elements(self, inset:float = 0) -> "Tileable":
        """Returns a new Tileable with an inset applied around the edges of the 
        tiling elements.

        Works by applying a negative buffer of specfied size to all elements.
        Elements that collapse to zero area are removed and the element_id
        attribute updated accordingly.

        NOTE: this method is likely to not preserve the relative area of 
        elements.

        Args:
            inset (float, optional): The distance to inset. Defaults to `0`.
        
        Returns:
            "Tileable": the new inset Tileable.
        """
        inset_elements, inset_ids = [], []
        for p, id in zip(self.elements.geometry, self.elements.element_id):
            b = p.buffer(-inset, resolution = 1, join_style = 2)
            if not b.area <= 0:
                inset_elements.append(b)
                inset_ids.append(id)
        result = copy.deepcopy(self)
        result.elements = gpd.GeoDataFrame(
            data = {"element_id": inset_ids}, crs = self.crs,
            geometry = gpd.GeoSeries(inset_elements))
        return result
    

    def plot(self, ax = None, show_tile:bool = True, show_reg_tile:bool = True, 
             show_ids:str = "element_id", show_vectors:bool = False, r:int = 0,
             tile_edgecolor:str = "k", reg_tile_edgcolor:str = "r", 
             r_alpha:float = 0.3, cmap:list[str] = None, 
             figsize:tuple[float] = (8, 8), **kwargs) -> pyplot.axes:
        """Plots a representation of the Tileable on the supplied axis. **kwargs
        are passed on to matplotlib.plot()

        Args:
            ax (_type_, optional): matplotlib axis to draw to. Defaults to None.
            show_tile (bool, optional): if `True` show the tile outline. 
                Defaults to `True`.
            show_reg_tile (bool, optional): if `True` show the regularised tile 
                outline. Defaults to `True`.
            show_ids (str, optional): if `element_id` show the element_ids. If
                `id` show index number. If None or `""` don't label elements.
                Defaults to `element_id`.
            show_vectors (bool, optional): if `True` show the translation 
                vectors (not the minimal pair, but those used by 
                `get_local_patch()`). Defaults to `False`.
            r (int, optional): passed to `get_local_patch()` to show context if 
                greater than 0. Defaults to `0`.
            r_alpha (float, optional): alpha setting for units other than the 
                central one. Defaults to 0.3.
            tile_edgecolor (str, optional): outline colour for the tile. 
                Defaults to `"k"`.
            reg_tile_edgcolor (str, optional): outline colour for the 
                regularised. Defaults to `"r"`.
            cmap (list[str], optional): colour map to apply to the central 
                tiles. Defaults to `None`.
            figsize (tuple[float], optional): size of the figure. 
                Defaults to `(8, 8)`.
        """
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
        if show_ids != None and show_ids != "":
            do_label = True
            if show_ids == "element_id" or show_ids == True:
                labels = self.elements.element_id
            elif show_ids == "id":
                labels = [str(i) for i in range(self.elements.shape[0])]
            else:
                do_label = False
            if do_label:
                for id, tile in zip(labels, self.elements.geometry):
                    ax.annotate(id, (tile.centroid.x, tile.centroid.y),
                                ha = "center", va = "center",
                                bbox = {"lw": 0, "fc": "#ffffff40"})
        if r > 0:
            self.get_local_patch(r = r).plot(ax = ax, column = "element_id",
                                             alpha = r_alpha, cmap = cm,
                                             **kwargs)
        if show_tile:
            self.tile.plot(ax = ax, ec = tile_edgecolor, lw = 0.5,
                           fc = "#00000000", **kwargs) 
        if show_vectors:  # note that arrows in mpl are dimensioned in plotspace
            for v in self.vectors[:len(self.vectors) // 2]:
                ax.arrow(0, 0, v[0], v[1], color = "k", width = w * 0.002,
                         head_width = w * 0.05, length_includes_head = True)
        if show_reg_tile:
            self.regularised_tile.plot(ax = ax, ec = reg_tile_edgcolor, 
                                       fc = "#00000000", lw = 2, **kwargs)
        return ax


    def _get_legend_elements(self):
        """Returns the elements augmented by a rotation column.
        
        This base implementation may be overridden by specific tile unit types.
        In particular see 
        `weavingspace.weave_unit.WeaveUnit._get_legend_elements()`.
        """
        elements = copy.deepcopy(self.elements)
        elements["rotation"] = 0
        return elements
    
    
    def transform_scale(self, xscale:float = 1.0, 
                        yscale:float = 1.0) -> "Tileable":
        """Transforms tileable by scaling.

        Args:
            xscale (float, optional): x scale factor. Defaults to 1.0.
            yscale (float, optional): y scale factor. Defaults to 1.0.

        Returns:
            Tileable: the transformed Tileable.
        """
        result = copy.deepcopy(self)
        result.elements.geometry = self.elements.geometry.scale(
            xscale, yscale, origin = (0, 0))
        result.tile.geometry = self.tile.geometry.scale(
            xscale, yscale, origin = (0, 0))
        result.regularised_tile.geometry = self.regularised_tile.geometry.scale(
            xscale, yscale, origin = (0, 0))
        result.vectors = result.get_vectors()
        return result
    
    
    def transform_rotate(self, angle:float = 0.0) -> "Tileable":
        """Transforms tiling by rotation.

        Args:
            angle (float, optional): angle to rotate by. Defaults to 0.0.

        Returns:
            Tileable: the transformed Tileable.
        """
        result = copy.deepcopy(self)
        result.elements.geometry = self.elements.geometry.rotate(
            angle, origin = (0, 0))
        result.tile.geometry = self.tile.geometry.rotate(
            angle, origin = (0, 0))
        result.regularised_tile.geometry = \
            self.regularised_tile.geometry.rotate(angle, origin = (0, 0))
        result.vectors = result.get_vectors()
        result.rotation = result.rotation + angle
        return result
    
    
    def transform_skew(self, xa:float = 0.0, ya:float = 0.0) -> "Tileable":
        """Transforms tiling by skewing

        Args:
            xa (float, optional): x direction skew. Defaults to 0.0.
            ya (float, optional): y direction skew. Defaults to 0.0.

        Returns:
            Tileable: the transformed Tileable.
        """
        result = copy.deepcopy(self)
        result.elements.geometry = self.elements.geometry.skew(
            xa, ya, origin = (0, 0))
        result.tile.geometry = self.tile.geometry.skew(
            xa, ya, origin = (0, 0))
        result.regularised_tile.geometry = \
            self.regularised_tile.geometry.skew(xa, ya, origin = (0, 0))
        result.vectors = result.get_vectors()
        return result
