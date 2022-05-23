#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass

import geopandas as gpd
import pandas as pd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine


def get_regular_polygon(spacing, n) -> geom.Polygon: 
    R = spacing / np.cos(np.radians(180 / n)) / 2
    a0 = -90 + 180 / n
    a_diff = 360 / n
    angles = [a0 + a * a_diff for a in range(n)]
    corners = [(R * np.cos(np.radians(a)), 
                R * np.sin(np.radians(a))) for a in angles]
    return geom.Polygon(corners)


def incentre(tri:geom.Polygon) -> geom.Point:
    corners = [geom.Point(p[0], p[1]) 
               for p in list(tri.exterior.coords)]  # ABCA
    lengths = [p.distance(q) 
               for p, q in zip(corners[:-1], corners[1:])]  # AB, BC, CA
    lengths = lengths[1:] + lengths[:1]  # BC, CA, AB 
    perimeter = np.sum(lengths)
    incentre = [0, 0]
    for p, l in zip(corners[:-1], lengths):
        incentre[0] = incentre[0] + p.x * l / perimeter
        incentre[1] = incentre[1] + p.y * l / perimeter
    return geom.Point(incentre)
    
