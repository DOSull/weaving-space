#!/usr/bin/env python
# coding: utf-8

"""The `WeaveUnit` subclass of `weavingspace.tileable.Tileable` implements
tileable geometric patterns constructed by specifying 2- or 3-axial weaves. 
            
Examples:
    Explain usage here...
"""

import logging
import itertools
from dataclasses import dataclass
from typing import Iterable, Union

import pandas as pd
import geopandas as gpd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine
import shapely.ops

import weavingspace.weave_matrices as weave_matrices
import weavingspace.tiling_utils as tiling_utils

from weavingspace._loom import Loom
from weavingspace._weave_grid import WeaveGrid

from weavingspace.tileable import TileShape
from weavingspace.tileable import Tileable


@dataclass
class WeaveUnit(Tileable):
    """Extends Tileable to allow for tiles that appear like woven patterns.

    Args:
        weave_type (str, optional): the type of weave pattern, one of 
            "plain",  "twill", "basket", "this", "cube" or "hex". Defaults
            to "plain".
        aspect (float, optional): width of strands relative to the spacing. 
            Defaults to 1.
        n (Union[int,tuple[int]]): number of over-under strands in biaxial 
            weaves. Only one item is required in a plain weave. Twill and 
            basket patterns expect an even number of elements in the tuple. 
            Defaults to (2, 2).
        strands (str, optional): specification of the strand labels 
            along each axis. Defaults to "a|b|c".
        debug (bool, optional): if True prints debug messages. Defaults to
            False.
    """  
    weave_type:str = "plain"
    aspect:float = 1.
    n:Union[int, tuple[int]] = (2, 2)
    strands:str = "a|b|c"
    _tie_up:np.ndarray = None
    _tr:np.ndarray = None
    _th:np.ndarray = None
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.weave_type = self.weave_type.lower()
    

    def _setup_tile_and_elements(self) -> None:
        """Returns dictionary with weave unit and tile GeoDataFrames based on
        parameters already supplied to the constructor.
        """
        self._parameter_info()

        if self.weave_type in ("hex", "cube"):
            self._setup_triaxial_weave_unit()
            self.tile_shape = TileShape.HEXAGON
        else:
            self._setup_biaxial_weave_unit()
            self.tile_shape = TileShape.RECTANGLE
        return


    def _parameter_info(self) -> None:
        """Outputs logging message concerning the supplied aspect settings.
        """    
        
        if self.aspect == 0:
            logging.info("Setting aspect to 0 is probably not a great plan.")

        if self.aspect < 0 or self.aspect > 1:
            logging.warning("""Values of aspect outside the range 0 to 1 won't produce tiles that will look like weaves, but they might be pretty anyway! Values less than -1 seem particularly promising, especially with opacity set less than 1.""")

        return None


    def _setup_biaxial_weave_unit(self) -> None:
        """Returns weave elements GeoDataFrame and tile GeoDataFrame in a dictionary based on paramters already supplied to constructor.
        """    
        warp_threads, weft_threads, _ = \
            tiling_utils.get_strand_ids(self.strands)
        
        if self.weave_type == "basket" and isinstance(self.n, (list, tuple)):
            self.n = self.n[0]
        
        p = weave_matrices.get_weave_pattern_matrix(
            weave_type = self.weave_type, n = self.n, warp = warp_threads,
            weft = weft_threads, tie_up = self._tie_up, tr = self._tr, 
            th = self._th)

        self._make_shapes_from_coded_weave_matrix(
            Loom(p), strand_labels = [weft_threads, warp_threads, []])


    def _get_triaxial_weave_matrices(self, 
            strands_1:Union[list[str],tuple[str]] = ["a"], 
            strands_2:Union[list[str],tuple[str]] = ["b"], 
            strands_3:Union[list[str],tuple[str]] = ["c"]
        ) -> Loom:
        """Returns encoded weave pattern matrix as Loom of three biaxial matrices.

        Allowed weave_types: "cube" or "hex".
        
        "hex" is not flexible and will fail with any strand label lists that are not length 3 and include more than one non-blank "-" item. You can generate the "hex" weave with the default settings in any case!
        
        Strand lists should be length 3 or length 1. "cube" tolerates more options than "hex" for the items in the strand lists.
        
        Defaults will produce 'mad weave'.
    
        Args:
            strands_1 (Union[list[str],tuple[str]], optional): list of labels 
                for warp strands. Defaults to ["a"].
            strands_2 (Union[list[str],tuple[str]], optional): list of labels 
                for weft strands. Defaults to ["b"].
            strands_3 (Union[list[str],tuple[str]], optional): list of labels 
                for weft strands. Defaults to ["c"].

        Returns:
            Loom: which combines the three biaxial weaves 12, 23 and 31 implied 
                by the strand label lists.
        """
        if self.weave_type == "hex":
            loom = Loom(
                weave_matrices.get_weave_pattern_matrix(
                    weave_type = "this", tie_up = np.ones((6, 6)), 
                    warp = strands_1, weft = strands_2),
                weave_matrices.get_weave_pattern_matrix(
                    weave_type = "this", tie_up = np.ones((6, 6)), 
                    warp = strands_2, weft = strands_3),
                weave_matrices.get_weave_pattern_matrix(
                    weave_type = "this", tie_up = np.ones((6, 6)), 
                    warp = strands_3, weft = strands_1),
            )
        else: # "cube"
            loom = Loom(  
                # Note n = (1,2,1,2) is required here to force 6x6 twill
                weave_matrices.get_weave_pattern_matrix(
                    weave_type = "twill", n = (1, 2, 1, 2), 
                    warp = strands_1, weft = strands_2),
                weave_matrices.get_weave_pattern_matrix(
                    weave_type = "twill", n = (1, 2, 1, 2), 
                    warp = strands_2, weft = strands_3),
                weave_matrices.get_weave_pattern_matrix(
                    weave_type = "twill", n = (1, 2, 1, 2), 
                    warp = strands_3, weft = strands_1),
            )
        return loom


    def _setup_triaxial_weave_unit(self) -> None:
        """Returns weave elements GeoDataFrame and tile GeoDataFrame in a dictionary based on parameters already supplied to constructor.
        """    
        strands_1, strands_2, strands_3 = \
            tiling_utils.get_strand_ids(self.strands)
        
        loom = self._get_triaxial_weave_matrices(
            strands_1 = strands_1, strands_2 = strands_2, strands_3 = strands_3)
        
        self._make_shapes_from_coded_weave_matrix(
            loom, strand_labels = [strands_1, strands_2, strands_3])


    # builds the geometric elements associated with a given weave supplied as
    # 'loom' containing the coordinates in an appropriate grid (Cartesian or
    # triangular) and the orderings of the strands at each coordinate location
    def _make_shapes_from_coded_weave_matrix(self, loom:Loom, 
            strand_labels:list[list[str]] = [["a"], ["b"], ["c"]]) -> None:
        """Returns weave elements GeoDataFrame and tile GeoDataFrame in a dictionary

        Args:
            loom (Loom): matrix or stack of matrices representing the weave 
                pattern.
            strand_labels (list[list[str]], optional): list of lists of labels
                for strands in each direction. Defaults to [["a"], ["b"], ["c"]]
        """  
        grid = WeaveGrid(loom.n_axes, loom.orientations, self.spacing)
        # expand the list of strand labels if needed in each direction
        labels = []
        for dim, thread in zip(loom.dimensions, strand_labels):
            labels.append(thread * int(np.ceil(dim // len(thread))))
        weave_polys, cells, strand_ids = [], [], []
        for coords, strand_order in zip(loom.indices, loom.orderings):
            ids = [thread[coord] for coord, thread in zip(coords, labels)]
            cells.append(grid.get_grid_cell_at(coords))
            if strand_order is None: continue  # No strands present
            if strand_order == "NA": continue  # Inconsistency in layer order
            n_slices = [len(id) for id in ids]
            next_polys = grid.get_visible_cell_strands(
                width = self.aspect, coords = coords,
                strand_order = strand_order, n_slices = n_slices)
            weave_polys.extend(next_polys)
            next_labels = [list(ids[i]) for i in strand_order]  # list of lists
            next_labels = list(itertools.chain(*next_labels))  # flatten 
            strand_ids.extend(next_labels)
        # sometimes empty polygons make it to here, so 
        # filter those out along with the associated IDs
        real_polys = [not p.is_empty for p in weave_polys]
        weave_polys = [p for p, b in zip(weave_polys, real_polys) if b]
        strand_ids = [id for id, b in zip(strand_ids, real_polys) if b]
        # note that the approx_tile is important for the 
        # biaxial case, which is not necessarily centred on (0, 0)
        approx_tile = shapely.ops.unary_union(cells)
        tile = grid.get_tile_from_cells(approx_tile)
        atc = approx_tile.centroid
        shift = (-atc.x, -atc.y)
        self.elements = self._get_weave_elements_gdf(
            weave_polys, strand_ids, shift)
        self.tile = gpd.GeoDataFrame(
            geometry = gpd.GeoSeries([tile]), crs = self.crs)
        return None


    def _get_weave_elements_gdf(
            self, polys:list[geom.Polygon], strand_ids:list[str], 
            offset:tuple[float]) -> gpd.GeoDataFrame:
        """Makes a GeoDataFrame from weave element polygons, labels, tile, etc.

        Args:
            polys (list[Polygon | MultiPolygon]): list of weave element 
                polygons.
            strand_ids (list[str]): list of strand labels.
            offset (tuple[float]): offset to centre the weave elements on the 
                tile.

        Returns:
            geopandas.GeoDataFrame: GeoDataFrame clipped to the tile, with 
                margin applied.
        """
        weave = gpd.GeoDataFrame(
            data = {"element_id": strand_ids},
            geometry = gpd.GeoSeries(
                [affine.translate(p, offset[0], offset[1]) for p in polys]))
        weave = weave[weave.element_id != "-"]
        # weave.geometry = tiling_utils.clean_polygon(weave.geometry)
        weave = weave.dissolve(by = "element_id", as_index = False)
        weave = weave.explode(index_parts = False, ignore_index = True)
        weave.geometry = tiling_utils.clean_polygon(weave.geometry)
        return weave.set_crs(self.crs)


    def _get_axis_from_label(self, label:str = "a", strands:str = None):
        """Determines the axis of an element_id from the strands spec string.

        Args:
            label (str, optional): the element_id. Defaults to "a".
            strands (str, optional): the strand spec. Defaults to the WeaveUnit
                strands attribute.

        Returns:
            _type_: the axis in which the supplied element_is found.
        """
        if strands == None:
            strands = self.strands
        index = strands.index(label)
        return strands[:index].count("|")


    def _get_legend_elements(self) -> gpd.GeoDataFrame:
        """Returns elements suitable for use in a legend representation.
        
        One element for each element_id value will be chosen, close to the
        centre of the tile extent, and not among the smallest elements present
        (for example not a short length of strand mostly hidden by other 
        elements)

        Returns:
            gpd.GeoDataFrame: the chosen elements.
        """
        angles = ((0, 240, 120) 
                  if self.weave_type in ("hex", "cube") 
                  else (90, 0))
        element_ids = pd.Series.unique(self.elements.element_id)
        groups = self.elements.groupby("element_id")
        elements, rotations = [], []
        for id in element_ids:
            candidates = groups.get_group(id)
            axis = self._get_axis_from_label(id, self.strands)
            elements.append(self._get_most_central_large_element(
                candidates, elements))
            rotations.append(-angles[axis] + self.rotation)
        return gpd.GeoDataFrame(
            data = {"element_id": element_ids, "rotation": rotations}, 
            crs = self.crs,
            geometry = gpd.GeoSeries(elements)
        )
        

    def _get_most_central_large_element(self, elements:gpd.GeoDataFrame,
                                        other_elements:list[geom.Polygon],
                                        ) -> geom.Polygon:
        """Gets a large element close to the centre of the WeaveUnit.

        Args:
            elements (gpd.GeoDataFrame): the set of elements to choose from.

        Returns:
            geom.Polygon: the chosen, large central element.
        """
        areas = [g.area for g in elements.geometry]
        min_area, max_area = min(areas), max(areas)
        if min_area / max_area > 0.5:
            geoms = list(elements.geometry)
        else:
            mean_log_a = np.mean(np.log(areas))
            geoms = [g for g, a in zip(elements.geometry, areas)
                            if np.log(a) > mean_log_a]
        if len(other_elements) == 0 or self.weave_type in ("cube", "hex"):
            d = [g.centroid.distance(geom.Point(0, 0)) for g in geoms]
        else:
            c = geom.MultiPolygon(other_elements).centroid
            d = [geom.MultiPolygon([g] + other_elements).centroid.distance(c)
                 for g in geoms]
        return geoms[d.index(min(d))]


    def _get_legend_key_shapes(self, polygon:geom.Polygon, 
                               counts:Iterable, angle:float = 0,
                               radial:bool = False) -> list[geom.Polygon]:
        """Returns a list of polygons obtained by slicing the supplied polygon
        across its length inton n slices. Orientation of the polygon is 
        indicated by the angle.
        
        The returned list of polygons can be used to form a colour ramp in a 
        legend.

        Args:
            polygon (geom.Polygon): the weave strand polygon to slice.
            n (int, optional): the number of slices required. Defaults to 25.
            angle (float, optional): orientation of the polygon. Defaults to 0.
            categorical (bool, optional): ignored by WeaveUnit.

        Returns:
            list[geom.Polygon]: a list of polygons.
        """
        c = polygon.centroid
        g = affine.rotate(polygon, -angle, origin = c)
        width, height, left, bottom = \
            tiling_utils.get_width_height_left_bottom(gpd.GeoSeries([g]))
        total = sum(counts)
        cuts = list(np.cumsum(counts))
        cuts = [0] + [c / total for c in cuts]
        cuts = [left + c * width for c in cuts]
        # add margin to avoid weird effects intersecting almost parallel lines.
        cuts[0] = cuts[0] - 1
        cuts[-1] = cuts[-1] + 1
        bottom = bottom - 1
        top = bottom + height + 1
        slices = []
        for l, r in zip(cuts[:-1], cuts[1:]):
            slice = geom.Polygon([(l, bottom), (r, bottom),
                                  (r, top), (l, top)])
            slices.append(slice.intersection(g)) 
        return [affine.rotate(s, angle, origin = c) for s in slices]


