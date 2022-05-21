#!/usr/bin/env python
# coding: utf-8

import logging
import itertools
from dataclasses import dataclass
from typing import Union
import copy

import geopandas as gpd
import numpy as np
import shapely.affinity as affine
import shapely.geometry as geom
import shapely.ops

import weave_matrices
import weave_utils
from loom import Loom
from weave_grids import WeaveGrid

from tile_units import TileShape
from tile_units import Tileable


@dataclass
class WeaveUnit(Tileable):
    """ Small data class containing elements of a weave unit.
    
    Attributes:
        elements: a GeoDataFrame of strand geometries.
        tile: a GeoDataFrame of the weave_unit tileable polygon (either a
            rectangle or a hexagon).
    """  
    weave_type:str = "plain"
    aspect:float = 1.
    margin:float = 0.
    n:Union[int, tuple[int]] = (2, 2)
    strands:str = "a|b|c"
    tie_up:np.ndarray = None
    tr:np.ndarray = None
    th:np.ndarray = None
    
    def __init__(self, **kwargs):
        """Constructor for WeaveUnit. Parameters are passed through to get_weave_unit function and also stored as instance attributes.
        
        Args:
            weave_type (str, optional): the type of weave pattern, one of 
                "plain",  "twill", "basket", "this", "cube" or "hex". Defaults
                to "plain".
            spacing (float, optional): spacing of threads in the weave in units 
                of the CRS. Defaults to 10000.
            aspect (float, optional): width of strands relative to the spacing. 
                Defaults to 1.
            margin (float, optional): margin applied to 'shrink' strand 
                elements, relative to the spacing. Defaults to 0.
            n (tuple of ints): number of over-under strands in biaxial weaves. 
                Only one item is required in a plain weave. Twill and basket patterns expect an even number of elements in the tuple. Defaults to (2, 2).
            strands (str, optional): specification of the strand labels 
                along each axis. Defaults to "a|b|c".
            tie_up (numpy.ndarray, optional): used when type is "this" to
                specify a desired weave pattern. See: Glassner A, 2002, “Digital weaving. 1” IEEE Computer Graphics and Applications 22(6) 108–118 DOI: 10.1109/MCG.2002.1046635. Defaults to None.
            tr (numpy.ndarray, optional): used when type is "this" to specify 
                the treadling matrix. See: Glassner 2002. Defaults to None.
            th (numpy.ndarray, optional): used when type is "this" to specify
                the threading matrix. See: Glassner 2002. Defaults to None.
            crs (int, optional): coordinate reference system. Usually an integer
                EPSG code, but any CRS object interpretable by geopandas will
                work. Defaults to 3857 (for Web Mercator).
        """
        unit = self._get_weave_unit(**kwargs)
        self.elements = unit["weave_unit"]
        self.tile = unit["tile"]
        self.regularised_tile = copy.copy(self.tile)
        for k, v in kwargs.items():
            self.__dict__[k] = v
        self.tile_shape = (TileShape.HEXAGON 
                           if self.weave_type in ("hex", "cube")
                           else TileShape.RECTANGLE) 
        self.vectors = self.get_vectors()
        self.regularise_elements()


    def _get_weave_unit(self, weave_type:str = "plain", spacing:float = 10000, 
            aspect:float = 1, margin:float = 0, n:tuple[int] = (2, 2), 
            strands:str = "a|b|c", tie_up:np.ndarray = None, 
            tr:np.ndarray = None, th:np.ndarray = None, crs:int = 3857
        ) -> dict:
        """Returns dictionary with weave unit and tile GeoDataFrames
        
        Args:
            weave_type (str, optional): the type of weave pattern, one of 
                "plain",  "twill", "basket", "this", "cube" or "hex". Defaults
                to "plain".
            spacing (float, optional): spacing of threads in the weave in units 
                of the CRS. Defaults to 10000.
            aspect (float, optional): width of strands relative to the spacing. 
                Defaults to 1.
            margin (float, optional): margin applied to 'shrink' strand 
                elements, relative to the spacing. Defaults to 0.
            n (tuple of ints): number of over-under strands in biaxial weaves. 
                Only one item is required in a plain weave. Twill and basket patterns expect an even number of elements in the tuple. Defaults to (2, 2).
            strands (str, optional): specification of the strand labels 
                along each axis. Defaults to "a|b|c".
            tie_up (numpy.ndarray, optional): used when type is "this" to
                specify a desired weave pattern. See: Glassner A, 2002, “Digital weaving. 1” IEEE Computer Graphics and Applications 22(6) 108–118 DOI: 10.1109/MCG.2002.1046635. Defaults to None.
            tr (numpy.ndarray, optional): used when type is "this" to specify 
                the treadling matrix. See: Glassner 2002. Defaults to None.
            th (numpy.ndarray, optional): used when type is "this" to specify
                the threading matrix. See: Glassner 2002. Defaults to None.
            crs (int, optional): coordinate reference system. Usually an integer
                EPSG code, but any CRS object interpretable by geopandas will
                work. Defaults to 3857 (for Web Mercator).

        Returns:
            dict: dictionary with contents {"weave_unit": GeoDataFrame of weave 
                elements, "tile": GeoDataFrame of the tile}.
        """
        self._parameter_info(margin, aspect)

        if weave_type in ("hex", "cube"):
            unit = get_triaxial_weave_unit(weave_type = weave_type,
                                           spacing = spacing, aspect = aspect,
                                           margin = margin, strands = strands, 
                                           crs = crs)
        else:
            unit = get_biaxial_weave_unit(weave_type = weave_type, n = n,
                                          spacing = spacing, aspect = aspect,
                                          margin = margin, strands = strands, 
                                          crs = crs, tie_up = tie_up, )
        return unit


    def _parameter_info(self, margin: float, aspect: float) -> None:
        """Outputs logging message concerning the supplied margin and aspect settings.

        Args:
            margin (float): weave unit margin.
            aspect (float): weave unit aspect.
        """    
        
        if aspect == 0:
            logging.info("Setting aspect to 0 is probably not a great plan.")

        if aspect < 0 or aspect > 1:
            logging.warning("""Values of aspect outside the range 0 to 1 won't 
                            produce tiles that will look like weaves, but they might be pretty anyway! Values less than -1 seem particularly promising, especially with opacity set less than 1.""")

        # maximum margin that will produce a weave-able tile
        max_margin = (1 - aspect) / 2
        if margin > max_margin:
            logging.warning(f"""With aspect set to {aspect:.3f} the largest margin 
                            that will work is {max_margin:.3f}. Lower values are required to produce proper tileable weaves. Specifically, with too wide a margin, strands in adjacent tiles will not 'join up' when tiled. Higher values will make nice tilings with broken strands, which aren't 'proper' weaves. The best alternative is to make the weave unit with margin = 0, then apply a negative buffer after you have tiled your map.""")   
        return None


def get_biaxial_weave_unit(
        spacing:float = 10_000., aspect:float = 1.,
        margin:float = 0., weave_type:str = "twill", crs:int = 3857, 
        n:Union[int, tuple[int]] = (2, 2), strands:str = "ab|cd", 
        tie_up:np.ndarray = weave_matrices.make_twill_matrix((2, 2)),
        tr:np.ndarray = None, th:np.ndarray = None
    ) -> dict:
    """Returns weave elements GeoDataFrame and tile GeoDataFrame in a dictionary

    Args:
        spacing (float, optional): spacing of strands. Defaults to 10_000.
        aspect (float, optional): width of strands relative to spacing.     
            Defaults to 1.0.
        margin (float, optional): inset margin relative to spacing. 
            Defaults to 0.0.
        weave_type (str, optional): one of "plain", "twill", "basket" or
            "this". Defaults to "twill".
        n (int | tuple[int], optional): over under pattern. See 
            make_over_under_row() for details. Defaults to 2.
        strands (str, optional): specification of strand labels. See 
            weaving_space_utils.get_strand_ids() for details. 
            Defaults to "ab|cd".
        crs (int, optional): CRS usually an EPSG int, but any geopandas CRS 
            object will work. Defaults to 3857.
        tie_up (np.ndarray, optional): a weave pattern matrix to pass thru in 
            the "this" case. Defaults to make_twill_matrix((2, 2)).
        tr (np.ndarray, optional): treadling matrix for the "this" case.        
            Defaults to None.
        th (np.ndarray, optional): threading matrix for the "this" case.        
            Defaults to None.

    Returns:
        dict: dictionary with contents {"weave_unit": GeoDataFrame of weave 
            elements, "tile": GeoDataFrame of the tile}.
    """    
    warp_threads, weft_threads, drop = weave_utils.get_strand_ids(strands)
    
    if weave_type == "basket" and isinstance(n, (list, tuple)):
        n = n[0]
    
    p = weave_matrices.get_weave_pattern_matrix(weave_type = weave_type, n = n, 
            warp = warp_threads, weft = weft_threads, 
            tie_up = tie_up, tr = tr, th = th)

    return make_shapes_from_coded_weave_matrix(
        Loom(p), spacing = spacing, width = aspect, margin = margin, 
        strand_labels = [weft_threads, warp_threads, []], crs = crs)


def get_triaxial_weave_matrices(
        weave_type:str = "cube", 
        strands_1:Union[list[str], tuple[str]] = ["a"], 
        strands_2:Union[list[str], tuple[str]] = ["b"], 
        strands_3:Union[list[str], tuple[str]] = ["c"]
    ) -> Loom:
    """Returns encoded weave pattern matrix as Loom of three biaxial matrices.
    
    See encode_biaxial_weave() for the encoding.

    Allowed weave_types: "cube" or "hex".
    
    "hex" is not flexible and will fail with any strand label lists that are not length 3 and include more than one non-blank "-" item. You can generate the "hex" weave with the default settings in any case!
    
    Strand lists should be length 3 or length 1. "cube" tolerates more options than "hex" for the items in the strand lists.
    
    Defaults will produce 'mad weave'.
 
    Args:
        weave_type (str, optional): one of "cube" or "hex". Defaults to "cube".
        strands_1 (list[str] | tuple[str], optional): list of labels for warp 
            strands. Defaults to ["a"].
        strands_2 (list[str] | tuple[str], optional): list of labels for weft 
            strands. Defaults to ["b"].
        strands_3 (list[str] | tuple[str], optional): list of labels for weft 
            strands. Defaults to ["c"].

    Returns:
        Loom: which combines the three biaxial weaves 12, 23 and 31 implied by 
            the strand label lists.
    """
    if weave_type == "hex":
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


def get_triaxial_weave_unit(
        spacing:float = 10_000., aspect:float = 1., 
        margin:float = 0., strands:str = "a|b|c", weave_type:str = "cube", crs:int = 3857
    ) -> dict:
    """Returns weave elements GeoDataFrame and tile GeoDataFrame in a dictionary

    Args:
        spacing (float, optional): spacing of strands. Defaults to 10_000.
        aspect (float, optional): width of strands relative to spacing.     
            Defaults to 1.0.
        margin (float, optional): inset margin relative to spacing. 
            Defaults to 0.0.
        weave_type (str, optional): one of "plain", "twill", "basket" or
            "this". Defaults to "twill".
        n (int | tuple[int], optional): over under pattern. See 
            make_over_under_row() for details. Defaults to 2.
        strands (str, optional): specification of strand labels. See 
            weaving_space_utils.get_strand_ids() for details. 
            Defaults to "a|b|c".
        crs (int, optional): CRS usually an EPSG int, but any geopandas CRS 
            object will work. Defaults to 3857.

    Returns:
        dict: dictionary with contents {"weave_unit": GeoDataFrame of weave 
            elements, "tile": GeoDataFrame of the tile}.
    """    
    strands_1, strands_2, strands_3 = weave_utils.get_strand_ids(strands)
    
    loom = get_triaxial_weave_matrices(
        weave_type = weave_type,
        strands_1 = strands_1, 
        strands_2 = strands_2,
        strands_3 = strands_3)
    
    return make_shapes_from_coded_weave_matrix(
        loom, spacing = spacing,
        width = aspect, margin = margin,
        strand_labels = [strands_1, strands_2, strands_3],
        crs = crs)


def centre_offset(shape: geom.Polygon, 
                    target:tuple[float] = (0, 0)) -> tuple[float]:
    """Returns vector required to move centroid of polygon to target. 

    Args:
        shape (Polygon): polygon to move.
        target (tuple[float], optional): target to move to. Defaults to (0, 0).

    Returns:
        tuple[float]: tuple of x, y movement required.
    """  
    shape_c = shape.centroid.coords[0]
    return (target[0] - shape_c[0], target[1] - shape_c[1])


# builds the sf associate with a given weave supplied as 'loom' which is a list
# containing the coordinates in an appropriate grid (Cartesian or triangular)
# and the orderings of the strands at each coordinate location
def make_shapes_from_coded_weave_matrix(
        loom:Loom, spacing:float = 1.0, width:float = 1.0, margin:float = 0.1,
        strand_labels:list[list[str]] = [["a"], ["b"], ["c"]], crs:int = 3857
    ) -> dict:
    """Returns weave elements GeoDataFrame and tile GeoDataFrame in a dictionary

    Args:
        loom (Loom): matrix or stack of matrices representing the weave pattern.
        spacing (float, optional): spacing of weave strands. Defaults to 1.0.
        width (float, optional): width of weave strands relative to spacing.
            Defaults to 1.0.
        margin (float, optional): inset margin applied to strands. Defaults 
            to 0.1.
        strand_labels (list[list[str]], optional): list of lists of labels
            for strands in each direction. Defaults to [["a"], ["b"], ["c"]].
        crs (int, optional): usually an EPSG code, but any object interpretable 
            by geopandas as a CRS. Defaults to 3857.

    Returns:
        dict: dictionary with contents {"weave_unit": GeoDataFrame of weave 
            elements, "tile": GeoDataFrame of the tile}
    """  
    grid = WeaveGrid(loom.n_axes, loom.orientations, spacing)
    # expand the list of strand labels if needed in each direction
    labels = []
    for dim, thread in zip(loom.dimensions, strand_labels):
        labels.append(thread * int(np.ceil(dim // len(thread))))
    weave_polys = []
    bb_polys = []
    strand_ids = []
    for coords, strand_order in zip(loom.indices, loom.orderings):
        ids = [thread[coord] for coord, thread in zip(coords, labels)]
        cell = grid.get_grid_cell_at(coords)
        bb_polys.append(grid._gridify(cell))
        # print(f"strand_order: {strand_order} ids: {ids}")
        if strand_order is None: continue  # No strands present
        if strand_order == "NA": continue  # Inconsistency in layer order
        n_slices = [len(id) for id in ids]
        next_polys = grid.get_visible_cell_strands(width, coords, 
                                                   strand_order, n_slices)
        weave_polys.extend(next_polys)
        next_labels = [list(ids[i]) for i in strand_order]  # a list of lists
        next_labels = list(itertools.chain(*next_labels))  # flatten 
        # print(f"n: {len(next_polys)} labels: {labels}")
        strand_ids.extend(next_labels)
    tile = shapely.ops.unary_union(bb_polys)
    shift = centre_offset(tile)
    tile = affine.translate(tile, shift[0], shift[1])
    return {
        "weave_unit": make_weave_gdf(weave_polys, strand_ids, tile, 
                                     shift, spacing, margin, crs),
        "tile": gpd.GeoDataFrame(geometry = gpd.GeoSeries([tile]),
                                       crs = crs)
    }


def make_weave_gdf(polys:list[Union[geom.Polygon, geom.MultiPolygon]], 
                   strand_ids:list[str], bb:geom.Polygon, 
                   offset:tuple[float], spacing:float, margin:float, 
                   crs:int) -> gpd.GeoDataFrame:
    """Makes a GeoDataFrame from weave element polygons, labels, tile, etc.

    Args:
        polys (list[Polygon | MultiPolygon]): list of weave element polygons.
        strand_ids (list[str]): list of strand labels.
        bb (Polygon): bounding shape tile.
        offset (tuple[float]): offset to centre the weave elements on the tile.
        spacing (float): spacing of weave strands.
        margin (float): inset margin to apply to strands.
        crs (int): CRS usually an EPSG code integer, but any geopandas
            acceptable CRS object will work.

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame clipped to the tile, with margin applied.
    """    
    weave = gpd.GeoDataFrame(
        data = {"element_id": strand_ids},
        geometry = gpd.GeoSeries(
            [affine.translate(p, offset[0], offset[1]) for p in polys]))
    weave = weave[weave.element_id != "-"]
    weave = weave.dissolve(by = "element_id", as_index = False)
    weave = weave.explode(index_parts = False, ignore_index = True)
    # this buffer operation cleans up some geometry issues
    weave.geometry = weave.buffer(-1e-4 * spacing)
    weave.geometry = weave.buffer((1e-4 - margin) * spacing)
    return weave.clip(bb).set_crs(crs)
