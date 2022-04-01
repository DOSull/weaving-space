#!/usr/bin/env python
# coding: utf-8

import numpy as np

from loom import Loom
from biaxial_weave_units import get_weave_pattern_matrix
from render_weave_grids import make_shapes_from_coded_weave_matrix
from weaving_space_utils import get_strand_ids

def get_triaxial_weave_matrices(*,
    weave_type = "cube", strands_1 = ["a"], strands_2 = ["b"], strands_3 = ["c"]):
    
    if weave_type == "hex":
        loom = Loom([
            get_weave_pattern_matrix(weave_type = "this", tie_up = np.ones((6, 6)), 
                                     warp = strands_1, weft = strands_2),
            get_weave_pattern_matrix(weave_type = "this", tie_up = np.ones((6, 6)), 
                                     warp = strands_2, weft = strands_3),
            get_weave_pattern_matrix(weave_type = "this", tie_up = np.ones((6, 6)), 
                                     warp = strands_3, weft = strands_1),
        ])
    else: # "cube"
        loom = Loom([
            get_weave_pattern_matrix(weave_type = "twill", n = (1, 2, 1, 2), 
                                     warp = strands_1, weft = strands_2),
            get_weave_pattern_matrix(weave_type = "twill", n = (1, 2, 1, 2), 
                                     warp = strands_2, weft = strands_3),
            get_weave_pattern_matrix(weave_type = "twill", n = (1, 2, 1, 2), 
                                     warp = strands_3, weft = strands_1),
        ])
    return loom


def get_triaxial_weave_unit(spacing = 10000, aspect = 1, margin = 0,
                            strands = "a|b|c", weave_type = "cube", crs = 3857):
    
    strand_ids = get_strand_ids(strands)
    strands_1 = strand_ids[0]
    strands_2 = strand_ids[1]
    strands_3 = strand_ids[2]
    
    loom = get_triaxial_weave_matrices(weave_type = weave_type,
                                       strands_1 = strands_1, 
                                       strands_2 = strands_2,
                                       strands_3 = strands_3)
    
    return make_shapes_from_coded_weave_matrix(loom, spacing = spacing,
                                               width = aspect, margin = margin,
                                               axis1_threads = strands_1,
                                               axis2_threads = strands_2,
                                               axis3_threads = strands_3, crs = crs)
