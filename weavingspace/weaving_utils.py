#!/usr/bin/env python
# coding: utf-8

# Some strand label parser stuff

import re

import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine


def _parse_strand_label(s:str) -> list[str]:
    """Breaks a strand label specifiction in to a list of labels

    Args:
        s (str): see get_strand_ids() for details
    Returns:
        list[str]: list of strand labels.
    """    
    clean_s = re.sub("[(]+", "(", re.sub("[)]+", ")", s))
    is_combo = False
    output = []
    current = ""
    for c in clean_s:
        if is_combo:
            if c == ")":
                output.append(current)
                current = ""
                is_combo = False
            else:
                current = current + c
        else:
            if c == "(":
                is_combo = True
            else:
                output.append(c)
    return output


def get_strand_ids(strands_spec: str) -> tuple[list[str]]:
    """Strands label string to split into strand labels.

    Args:
        strands_spec (str): string format "a|bc|(de)f" | separates strands in
            each direction and () designates combining labels into a single strand that will be sliced lengthwise. Example output:

                "a|bc|(de)f" -> (["a"], ["b", "c"], ["de", "f"])
            
            Superflous parentheses are removed, but no other error-checks are
            applied.

    Returns:
        tuple[str]: tuple of lists of labels for each set of strands.
    """    
    strand_ids = [_parse_strand_label(s) for s in strands_spec.split("|")]
    strand_ids = (strand_ids
                  if len(strand_ids) == 3
                  else strand_ids + [[""]])
    return tuple(strand_ids)


def centre_offset(shape: geom.Polygon, 
                  target:tuple[float] = (0, 0)) -> tuple[float]:
    """Returns vector required to move centroid of polygon to target. 

    Args:
        shape (Polygon): polygon to move.
        target (tuple[float], optional): target to move to. 
            Defaults to (0, 0).

    Returns:
        tuple[float]: tuple of x, y movement required.
    """  
    shape_c = shape.centroid.coords[0]
    return (target[0] - shape_c[0], target[1] - shape_c[1])


def get_axis_from_label(label:str = "a", strands:str = "a|b|c"):
    index = strands.index(label)
    return strands[:index].count("|")


def get_colour_ramp(geometry, n = 10, a = 0):
    c = geometry.centroid
    g = affine.rotate(geometry, -a, origin = c)
    bb = g.bounds
    cuts = np.linspace(bb[0], bb[2], n + 1)
    slices = []
    for l, r in zip(cuts[:-1], cuts[1:]):
        slice = geom.Polygon([(l, bb[1]), (r, bb[1]), (r, bb[3]), (l, bb[3])])
        slices.append(slice.intersection(g)) 
    return [affine.rotate(s, a, origin = c) for s in slices]
