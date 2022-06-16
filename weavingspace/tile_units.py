#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
import copy
import string

import geopandas as gpd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine

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
    # margin:float = 0.
    
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            self.__dict__[k] = v
        self.setup_tile_and_elements()
        self.setup_vectors()
        if self.regularised_tile is None:
            self.regularise_elements()
        return
        
        
    def setup_vectors(self) -> None:
        self.vectors = self.get_vectors()
        return None
    
        
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
        
        
    # Make up a regularised tile by carefully unioning the elements
    def make_regularised_tile_from_elements(self) -> None:
        self.regularised_tile = copy.deepcopy(self.tile)
        self.regularised_tile.geometry = tiling_utils.safe_union(
            self.elements.geometry)
        # This simplification seems very crude but fixes all kinds of issues...
        self.regularised_tile.geometry[0] = \
            self.regularised_tile.geometry[0].simplify(self.spacing / 100)
        return
        
    
    def merge_fragments(self, 
                        fragments:list[geom.Polygon]) -> list[geom.Polygon]:
        """
        Merges a set of polygons based on testing if they touch when subjected to the translation vectors provided by the get_vectors() method. Called by regularise_elements() method to combine elements in a tile that may be fragmented as supplied but will combine when tiled into single elements. This step makes for more efficient implementation of the tiling of map regions.

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
                        next_frags.append(u1)
                        reg_tile = (reg_tile | u1) - u2
                    else:
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
        return fragments


    def regularise_elements(self) -> None:
        """
        Combines separate elements that share an element_id value into single elements, if they would end up touching when tiled. Also adjusts the regularised_tile attribute according.
        """
        self.regularised_tile = copy.deepcopy(self.tile)
        self.regularised_tile.geometry = \
            self.regularised_tile.geometry.buffer(
                self.fudge_factor, resolution = 1, join_style = 2)
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
        # This simplification seems very crude but fixes all kinds of issues...
        self.regularised_tile.geometry[0] = \
            self.regularised_tile.geometry[0].simplify(self.spacing / 100)
        return None
    
    
    def get_local_patch(self, r:int = 1, 
                        include_0:bool = False) -> gpd.GeoDataFrame:
        """Returns a GeoDataFrame with translated copies of the Tileable.

        Args:
            r (int, optional): the number of translation vector steps required. 
                Defaults to 1.
            include_0 (bool, optional): If True includes the Tileable itself at 
                (0, 0). Defaults to False.

        Returns:
            gpd.GeoDataFrame: A GeoDataframe of the Tileable's elements extended
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
        vectors = self.get_vectors(return_values = False)
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

        Args:
            centre_element (int, optional): the index position of the central 
                element. Defaults to 0.
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
    def inset_elements(self, inset:float = 1) -> None:
        inset_elements, inset_ids = [], []
        for p, id in zip(self.elements.geometry, self.elements.element_id):
            b = p.buffer(-inset, resolution = 1, join_style = 2)
            if not b.area <= 0:
                inset_elements.append(b)
                inset_ids.append(id)
        self.elements = gpd.GeoDataFrame(
            data = {"element_id": inset_ids}, crs = self.crs,
            geometry = gpd.GeoSeries(inset_elements))
        return
    

    def plot(self, ax = None, show_tile:bool = True, show_reg_tile:bool = True, 
             show_ids = True, show_vectors:bool = False, r:int = 0,
             tile_edgecolor:str = "k", reg_tile_edgcolor:str = "r", 
             cmap:list[str] = None, figsize:tuple[float] = (8, 8), 
             **kwargs) -> None:
        """Plots a representation of the Tileable on the supplied axis. **wargs
        are passed on to matplotlib.plot()

        Args:
            ax (_type_, optional): matplotlib axis to draw to. Defaults to None.
            show_tile (bool, optional): if True show the tile outline. 
                Defaults to True.
            show_reg_tile (bool, optional): if True show the regularised tile 
                outline. Defaults to True.
            show_ids (bool, optional): if True show the element_ids. 
                Defaults to True.
            show_vectors (bool, optional): if true show the translation vectors 
                (not the minimal pair, but those used by get_local_patch). 
                Defaults to False.
            r (int, optional): passed to get_local_patch to show context if 
                greater than 0. Defaults to 0.
            tile_edgecolor (str, optional): outline colour for the tile. 
                Defaults to "k".
            reg_tile_edgcolor (str, optional): outline colour for the 
                regularised. Defaults to "r".
            cmap (list[str], optional): colour map to apply to the central 
                tiles. Defaults to None.
            figsize (tuple[float], optional): size of the figure. 
                Defaults to (8, 8).

        Returns:
            _type_: the axis.
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
        if show_ids:
            for id, tile in zip(self.elements.element_id, 
                                self.elements.geometry):
                ax.annotate(id, (tile.centroid.x, tile.centroid.y),
                            ha = "center", va = "center",
                            bbox = {"lw": 0, "fc": "#ffffff40"})
        if r > 0:
            self.get_local_patch(r = r).plot(ax = ax, column = "element_id",
                                             alpha = 0.3, cmap = cm, **kwargs)
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
        """
        elements = copy.deepcopy(self.elements)
        elements["rotation"] = 0
        return elements
        

@dataclass
class TileUnit(Tileable):
    tiling_type:str = None
    dissection_offset:int = 1
    n:int = 3
    code:str = "3.3.4.3.4"
        
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        if self.tile_shape == TileShape.TRIANGLE:
            self._modify_tile()
            self._modify_elements()
            self.setup_vectors()
            self.make_regularised_tile_from_elements()
        if self.regularised_tile is None:
            self.make_regularised_tile_from_elements()


    def setup_tile_and_elements(self) -> None:
        """Delegates setup of the unit to various functions depending
        on self.tiling_type.
        """
        if self.tiling_type == "cairo":
            tiling_geometries.setup_cairo(self)
        elif self.tiling_type == "hex-dissection":
            tiling_geometries.setup_hex_dissection(self)
        elif self.tiling_type == "laves":
            tiling_geometries.setup_laves(self)
        elif self.tiling_type == "archimedean":
            tiling_geometries.setup_archimedean(self)
        elif self.tiling_type in ("hex-colouring", "hex-coloring"):
            tiling_geometries.setup_hex_colouring(self)
        elif self.tiling_type in ("square-colouring", "square-coloring"):
            tiling_geometries.setup_square_colouring(self)
        else:
            tiling_geometries.setup_none_tile(self)
        return

    
    def _modify_elements(self) -> None:
        """It is not trivial to tile a triangle, so this function augments
        augments the elements of a triangular tile to a diamond by 180 degree 
        rotation. Operation is 'in place'.
        """
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
        return None


    def _modify_tile(self) -> None:
        """It is not trivial to tile a triangular tile so this function 
        changes the tile to a diamond by copying and joining a 180 degree
        rotated copy. Operation is 'in-place'.
        """
        tile = self.tile.geometry[0]
        # translate to sit on x-axis
        tile = affine.translate(tile, 0, -tile.bounds[1])
        # make rotated copies
        # buffering applied to ensure union 'sticks'
        twins = [affine.rotate(tile, a, origin = (0, 0)) \
                    .buffer(self.fudge_factor, resolution = 1, join_style = 2) 
                 for a in range(0, 360, 180)]
        # and here we undo the buffer
        merged_tile = gpd.GeoSeries(twins).unary_union.buffer(
                -self.fudge_factor, resolution = 1, join_style = 2)
        self.tile_shape = TileShape.DIAMOND
        self.tile.geometry = gpd.GeoSeries([merged_tile])
        return None

    
    def _get_legend_key_shapes(self, polygon:geom.Polygon, 
                               counts:int = 25, angle:float = 0, 
                               categorical:bool = False) -> list[geom.Polygon]:
        """Returns a set of shapes that can be used to make a legend key 
        symbol for the supplied polygon. In TileUnit this is a set of 'nested
        polygons.

        Args:
            polygon (geom.Polygon): the polygon to symbolise.
            n (int, optional): number of steps. Defaults to 25.
            rot (float, optional): rotation that may have to be applied.  
                Not used in the TileUnit case. Defaults to 0.

        Returns:
            list[geom.Polygon]: a list of nested polygons.
        """
        if not categorical:
            n = sum(counts)
            bandwidths = list(np.cumsum(counts))
            bandwidths = [bw / n for bw in bandwidths]
            bandwidths = [bw if bw > 0.05 else 0.05 for bw in bandwidths]
            n = sum(bandwidths)
            bandwidths = [0] + [bw / n for bw in bandwidths]
            # # make buffer widths that will yield approx equal area 'annuli'
            # bandwidths = range(n_steps + 1)
            # sqrt exaggerates outermost annuli, which can otherwise disappear
            bandwidths = [np.sqrt(bw) for bw in bandwidths]
            distances = np.cumsum(bandwidths)
            # get the negative buffer distance that will 'collapse' the polygon
            radius = tiling_utils.get_collapse_distance(polygon)
            distances = distances * radius / distances[-1]
            nested_polys = [polygon.buffer(-d, resolution = 1, 
                                           join_style = 2) for d in distances]
            # return converted to annuli (who knows someone might set alpha < 1)
            return [g1.difference(g2) for g1, g2 in 
                    zip(nested_polys[:-1], nested_polys[1:])]      
        else:
            n = sum(counts)
            slice_posns = list(np.cumsum(counts))
            slice_posns = [0] + [p / n for p in slice_posns]
            return [tiling_utils.get_polygon_sector(polygon, i, j) 
                    for i, j in zip(slice_posns[:-1], slice_posns[1:])]


    def inset_tile(self, d:float = 0) -> None:
        self.elements = self.elements.clip(
            self.regularised_tile.geometry.buffer(
                -d, resolution = 1, join_style = 2))
        return
    
    
    def scale_elements(self, sf:float = 1) -> None:
        self.elements.geometry = self.elements.geometry.scale(
            sf, sf, origin = (0, 0))
        return


