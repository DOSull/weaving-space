#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
from collections import defaultdict
from typing import Union

import geopandas as gpd
import shapely.geometry as geom
import shapely.affinity as affine

import matplotlib.pyplot as pyplot

from weavingspace import TileUnit
from weavingspace import TileShape
import weavingspace.tiling_utils as tiling_utils

@dataclass
class Topology:
  """Class to represent topology of a Tileable.
  """

  unit: TileUnit = None
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

  def __init__(self, unit: TileUnit):
    self.unit = unit
    rad = 1 if unit.tile_shape == TileShape.HEXAGON else 2
    self._patch = self.unit.get_local_patch(r = rad, include_0 = True)
    self.polygons = tuple(self._patch.geometry)
    self.centres = tuple([tiling_utils.incentre(p) for p in self.polygons])
    self._set_vertices()
    self._set_edges_and_polygons()
    self._set_vertex_neighbours_and_polygons()


  def _set_vertices(self) -> None:
    """Sets up the list of vertices to include only unique vertices, i.e. 
    vertices shared between polygons are stored only once.
    """
    vertices = []
    for p in self.polygons:
      corners = tiling_utils.get_corners(p, repeat_first = False)
      for c in corners:
        # only add to vertices if not close to an existing vertex
        if not any ([c.distance(v) <= 2 * tiling_utils.RESOLUTION 
                     for v in vertices]):
          vertices.append(c)
    self.vertices = tuple(vertices)


  def _get_closest_vertex(self, v:Union[int, geom.Point], 
                          return_id:bool = True) -> Union[int, geom.Point]:
    """Finds the closest vertex to that supplied from the Topology's vertex 
    list.

    Args:
      v (Union[int, geom.Point]): vertex to find closest other vertex from, 
        supplied either as an integer index or a Point geometry.
      return_id (bool, optional): if True returns the closest vertex id, if  
        False the point. Defaults to True.

    Returns:
      Union[int, geom.Point]: the closest vertex in self.vertices as an index
        or Point geometry.
    """
    dists_to_vertices = sorted([(v.distance(vertex), i) 
                                for i, vertex in enumerate(self.vertices)])
    if return_id:
      return dists_to_vertices[0][1]
    else:
      return self.vertices[dists_to_vertices[0][1]]


  def _set_edges_and_polygons(self) -> None:
    """Sets up the edge and polygon collection.
    """
    # empty collections as lists and dictionaries so we can build them
    edges = []              # tiling edges as vertex indices 
    lefts, rights = {}, {}  # left and right neighboring tiles of each edge
    poly_edges = []         # polygons as sequences of edge indices
                            # edges traversed backwards encoded as -idx - 1
    poly_vertices = []      # polygons as sequences of vertex indices
    for p in self._patch.geometry:
      # get all the vertices close to this polygon - most will be its corners
      # but some will be corners of neighbouring tiles that lie on its edges
      # convenient to have the ids and the points in separate lists
      near_ids = [i for i, v in enumerate(self.vertices)
                  if p.distance(v) <= tiling_utils.RESOLUTION]
      near_vs = [self.vertices[i] for i in near_ids]
      p_edges = tiling_utils.get_edges(p)
      p_edge_ids, p_vertex_ids = [], []
      for e in p_edges:
        ends = [geom.Point(p) for p in e.coords]
        edge_v_ids = [self._get_closest_vertex(v) for v in ends] 
        # now find any potential induced vertices to insert along the edge
        inserts = [idx for idx, v in zip(near_ids, near_vs)
                   if idx not in edge_v_ids and
                   e.distance(v) <= tiling_utils.RESOLUTION]
        # might (rarely) be more than one so sort them by distance along edge
        dist_along = [e.line_locate_point(self.vertices[idx], normalized = True)
                      for idx in inserts]
        inserts = [inserts[dist_along.index(p)] for p in sorted(dist_along)]
        # and insert them
        v_seq = edge_v_ids[:1] + inserts + edge_v_ids[-1:]
        edge_segments = [(x, y) for x, y in zip(v_seq[:-1], v_seq[1:])]
        for segment in edge_segments:
          # add first vertex of segment to current polygon vertex list
          p_vertex_ids.append(segment[0])
          if segment in edges:
            # already in edge list so only add to current polygon edge list
            p_edge_ids.append(edges.index(segment))
            # rights[edges.index(segment)] = len(poly_edges)
          elif tuple(reversed(segment)) in edges:
            # reverse direction is in edge list so add it appropriately coded 
            # to current polygon edge list and record left neighbour
            rev_segment = tuple(reversed(segment))
            # use topojson convention to indicate reversal of direction
            p_edge_ids.append(-edges.index(rev_segment) - 1)
            lefts[edges.index(rev_segment)] = len(poly_edges)
          else: 
            # new edge so add everywhere!
            edges.append(segment)
            p_edge_ids.append(edges.index(segment))
            rights[edges.index(segment)] = len(poly_edges)
      # convert to tuples and add to the polygon lists
      poly_edges.append(tuple(p_edge_ids))
      poly_vertices.append(tuple(p_vertex_ids))
    self.edges = tuple(edges)
    self.polygon_edges = tuple(poly_edges)
    self.polygon_vertices = tuple(poly_vertices)
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


  def get_dual_tiles(self) -> gpd.GeoDataFrame:
    """Alternative approach to creating the dual tiie unit for a tile unit.

    Args:
        t (TileUnit): _description_

    Returns:
        gpd.GeoDataFrame: _description_
    """
    dual_vertices = [i for i, vps in enumerate(self.vertex_polygons) 
                    if len(vps) > 2]
    dual_faces = gpd.GeoSeries(
      [geom.Polygon([self.centres[i] for i in p])
      for p in [self.vertex_polygons[i] for i in dual_vertices]])
    
    # TODO: label and select the faces to include so that they are tileable...
    return gpd.GeoSeries(dual_faces)


  def plot(self, show_original_tiles: bool = True, 
                 show_tile_centres: bool = True,
                 show_all_vertices: bool = False,
                 show_all_edges: bool = False,
                 show_dual_polygons: bool = False):
    """Plots the topology.

    Args:
        show_original_tiles (bool, optional): _description_. Defaults to True.
        show_tile_centres (bool, optional): _description_. Defaults to True.
        show_all_vertices (bool, optional): _description_. Defaults to False.
    """
    fig = pyplot.figure(figsize = (10, 10))
    ax = fig.add_axes(111)
    
    if show_original_tiles:
      gpd.GeoSeries(self.polygons).plot(ax = ax, fc = "#80808000", 
                                        ec = "k", lw = 0.35)

    if show_tile_centres:
      for i, c in enumerate(self.centres):
        ax.annotate(text = i, xy = (c.x, c.y), fontsize = 10, 
                    color = "b", ha = "center", va = "center")
    
    tiling_vertex_ids = [i for i, v in enumerate(self.vertices)
                         if len(self.vertex_polygons) > 2]
    tiling_vertices = [self.vertices[i] for i in tiling_vertex_ids]
    if show_all_vertices:
      gpd.GeoSeries(self.vertices).plot(ax = ax, color = "r", markersize = 5)
      for i, v in enumerate(self.vertices):
        ax.annotate(text = i, xy = (v.x, v.y), fontsize = 8, color = "r",
                    ha = "right", va = "bottom")
    else:
      gpd.GeoSeries(tiling_vertices).plot(ax = ax, color = "r", markersize = 5)
      for i, v in zip(tiling_vertex_ids, tiling_vertices):
        ax.annotate(text = i, xy = (v.x, v.y), fontsize = 8, color = "r",
                    ha = "right", va = "bottom")
        
    edges = gpd.GeoSeries([affine.scale(
      geom.LineString(
        [self.vertices[self.edges[i][0]], 
         self.vertices[self.edges[i][1]]]).parallel_offset(35, side = "right"), 0.9, 0.9) for i in range(len(self.edges))]) 
    edges = gpd.GeoDataFrame({"id": range(len(self.edges))}, geometry = edges)
    if not show_all_edges:
      lr = set(self.lefts.keys()).intersection(set(self.rights.keys()))
      edges = edges.iloc[sorted(list(lr))]
    edges.plot(ax = ax, color = "g", lw = 0.75, ls = "dashed")
    for i, e in zip(edges.id, edges.geometry):
      ax.annotate(text = i, xy = (e.centroid.x, e.centroid.y), fontsize = 8,
                  ha = "center", va = "center")
      
    if show_dual_polygons:
      self.get_dual_tiles().buffer(-15).plot(ax = ax, fc = "dodgerblue", 
                                             alpha = 0.2)

    pyplot.axis("off")
    return ax
