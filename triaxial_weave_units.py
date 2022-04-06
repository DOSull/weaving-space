#!/usr/bin/env python
# coding: utf-8

import numpy as np

from loom import Loom
from biaxial_weave_units import get_weave_pattern_matrix
from render_weave_grids import make_shapes_from_coded_weave_matrix
from weaving_space_utils import get_strand_ids

def get_triaxial_weave_matrices(*,
                                weave_type:str = "cube", 
                                strands_1:list[str]|tuple[str] = ["a"], 
                                strands_2:list[str]|tuple[str] = ["b"], 
                                strands_3:list[str]|tuple[str] = ["c"]
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
            get_weave_pattern_matrix(weave_type = "this", 
                                     tie_up = np.ones((6, 6)), 
                                     warp = strands_1, weft = strands_2),
            get_weave_pattern_matrix(weave_type = "this", 
                                     tie_up = np.ones((6, 6)), 
                                     warp = strands_2, weft = strands_3),
            get_weave_pattern_matrix(weave_type = "this", 
                                     tie_up = np.ones((6, 6)), 
                                     warp = strands_3, weft = strands_1),
        )
    else: # "cube"
        loom = Loom(  
            # Note n = (1,2,1,2) is required here to force 6x6 twill
            get_weave_pattern_matrix(weave_type = "twill", n = (1, 2, 1, 2), 
                                     warp = strands_1, weft = strands_2),
            get_weave_pattern_matrix(weave_type = "twill", n = (1, 2, 1, 2), 
                                     warp = strands_2, weft = strands_3),
            get_weave_pattern_matrix(weave_type = "twill", n = (1, 2, 1, 2), 
                                     warp = strands_3, weft = strands_1),
        )
    return loom


def get_triaxial_weave_unit(*, spacing:float = 10_000., aspect:float = 1., 
                            margin:float = 0., strands:str = "a|b|c", weave_type:str = "cube", crs:int = 3857) -> dict:
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
    strands_1, strands_2, strands_3 = get_strand_ids(strands)
    
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
