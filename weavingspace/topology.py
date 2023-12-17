#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
from collections import defaultdict

import geopandas as gpd
import shapely.geometry as geom
import numpy as np

from weavingspace.tileable import Tileable
import weavingspace.tiling_utils as tiling_utils

@dataclass
class Topology:
  """Class to represent topology of a Tileable.
  """

  unit: Tileable = None
  _patch: gpd.GeoDataFrame = None
  vertices: tuple[geom.Point] = None
  vertex_neighbours: tuple[int] = None
  vertex_polygons: tuple[int] = None
  edges: tuple[tuple[int]] = None
  lefts: dict[int, int] = None
  rights: dict[int, int] = None
  polygons: tuple[geom.Polygon] = None
  polygon_edges: tuple[tuple[int]] = None
  polygon_vertices: tuple[tuple[int]] = None
  centres: tuple[geom.Point] = None

  def __init__(self, unit:Tileable):
    self.unit = unit
    self._patch = self.unit.get_local_patch(r = 2, include_0 = True)
    self.polygons = tuple(self._patch.geometry)
    self.centres = tuple([tiling_utils.incentre(p) for p in self.polygons])
    self._set_vertices()
    self._set_edges_and_polygons()
    self._set_vertex_neighbours_and_polygons()


  def _set_vertices(self) -> None:
    vertices = []
    for p in self.polygons:
      corners = tiling_utils.get_corners(p, repeat_first = False)
      for c in corners:
        if not c in vertices:
          vertices.append(c)
    self.vertices = tuple(vertices)


  def _set_edges_and_polygons(self) -> None:
    edges = []
    lefts, rights = {}, {}
    polygon_edges = []
    p_vertices = []
    for p in self._patch.geometry:
      near_vertex_ids = [i for i, v in enumerate(self.vertices)
                         if p.distance(v) <= tiling_utils.RESOLUTION]
      near_vertices = [self.vertices[i] for i in near_vertex_ids]
      p_edges = tiling_utils.get_edges(p)
      p_edge_ids = []
      p_vertex_ids = []
      for e in p_edges:
        e_vertex_ids = [(self.vertices.index(geom.Point(v))) for v in e.coords]
        inserts = [idx for idx, v in zip(near_vertex_ids, near_vertices)
                   if idx not in e_vertex_ids and
                   e.distance(v) <= tiling_utils.RESOLUTION]
        x_along = [e.line_locate_point(self.vertices[idx], normalized = True)
                  for idx in inserts]
        sorted_inserts = [inserts[x_along.index(p)] for p in sorted(x_along)]
        vertex_sequence = \
          e_vertex_ids[:1] + sorted_inserts + e_vertex_ids[-1:]
        edge_segments = [(x, y) for x, y
                         in  zip(vertex_sequence[:-1], vertex_sequence[1:])]
        for segment in edge_segments:
          if segment in edges:
            p_edge_ids.append(edges.index(segment))
            rights[edges.index(segment)] = len(polygon_edges)
          elif tuple(reversed(segment)) in edges:
            # use the topojson convention to indicate reversal of direction
            rev_segment = tuple(reversed(segment))
            p_edge_ids.append(-edges.index(rev_segment) - 1)
            lefts[edges.index(rev_segment)] = len(polygon_edges)
          else:
            edges.append(segment)
            p_edge_ids.append(edges.index(segment))
            rights[edges.index(segment)] = len(polygon_edges)
          p_vertex_ids.append(segment[0])
      polygon_edges.append(tuple(p_edge_ids))
      p_vertices.append(tuple(p_vertex_ids))
    self.edges = tuple(edges)
    self.polygon_edges = tuple(polygon_edges)
    self.polygon_vertices = tuple(p_vertices)
    self.lefts = lefts
    self.rights = rights


  def _set_vertex_neighbours_and_polygons(self) -> None:
    vertex_neighbours = defaultdict(list)
    for i, e in enumerate(self.edges):
      vertex_neighbours[self.vertices[e[0]]].append(self.vertices[e[1]])
      vertex_neighbours[self.vertices[e[1]]].append(self.vertices[e[0]])
    ordered_neighbours = defaultdict(list)
    for p0, neighbours in vertex_neighbours.items():
      order = tiling_utils.order_pts_cw_relative_to_centre(neighbours, p0)
      ordered_neighbours[self.vertices.index(p0)] = \
        [self.vertices.index(neighbours[i]) for i in order]
    self.vertex_neighbours = tuple(tuple([ordered_neighbours[i] 
                                   for i in range(len(self.vertices))]))
    incident_polygons = []
    for i, neighbours in enumerate(self.vertex_neighbours):
      polys = []
      for j in neighbours:
        if (i, j) in self.edges:
          polys.append(self.rights[self.edges.index((i, j))])
        elif (j, i) in self.edges:
          if self.edges.index((j, i)) in self.lefts:
            polys.append(self.lefts[self.edges.index((j, i))])
      incident_polygons.append(polys)
    self.vertex_polygons = tuple(incident_polygons)
