#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
from collections import defaultdict
from typing import Union

import numpy as np
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

  tile_unit: TileUnit = None
  _patch: gpd.GeoDataFrame = None

  points: list[geom.Point] = None
  point_tiles: list[int] = None
  
  # corners: set = None
  # vertices: set = None
  # point_types: list[str] = None

  tiles: list[geom.Polygon] = None

  # outer_vertices: list[int] = None
  # vertices: list[int] = None
  # corners: list[int] = None
  vertex_neighbours: list[int] = None
  vertex_tiles: list[int] = None

  edges: list[list[int]] = None
  edge_lefts: dict[int, int] = None
  edge_rights: dict[int, int] = None

  tile_centres: list[geom.Point] = None
  tile_edges: list[tuple[int]] = None
  tile_vertices: list[tuple[int]] = None

  # TODO: think about how to add these in a coherent fashion!
  # sides: tuple[tuple[int]] = None
  # tile_corners: tuple[tuple[int]] = None
  # tile_sides: tuple[tuple[int]] = None

  def __init__(self, unit: TileUnit):
    self.tile_unit = unit
    rad = 1 if unit.base_shape == TileShape.HEXAGON else 2
    self._patch = self.tile_unit.get_local_patch(r = rad, include_0 = True)
    self.tiles = tuple(self._patch.geometry)
    self.tile_centres = tuple([tiling_utils.incentre(p) for p in self.tiles])
    self._setup_points()
    self._setup_tiles()


  def _setup_points(self) -> None:
    """Sets up the list of vertices to include only unique vertices, i.e. 
    vertices shared between polygons are stored only once.
    """
    self.points = []
    for p in self.tiles:
      corners = tiling_utils.get_corners(p, repeat_first = False)
      for c in corners:
        # only add if not close to an existing point
        if not any ([c.distance(v) <= 2 * tiling_utils.RESOLUTION 
                     for v in self.points]):
          self.points.append(c)


  def _get_closest_point(self, v:Union[int, geom.Point], 
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
    distances_to_points = sorted([(v.distance(vertex), i) 
                                  for i, vertex in enumerate(self.points)])
    if return_id:
      return distances_to_points[0][1]
    else:
      return self.points[distances_to_points[0][1]]


  def _setup_tiles(self) -> None:
    """Sets up the edge and tile collections.
    """
    self.point_tiles = defaultdict(list) # tiles incident on each point
    self.corners = set()                 # ids of corners
    self.edges = []                      # edges as tuples of point ids
    self.edge_lefts = {}                 # left and right tiles of each edge  
    self.edge_rights = {}
    self.tile_vertices = []              # tiles as tuples of vertces (these 
                                         # will not admit correct drawing)
    self.tile_corners = []               # tile corners as tuples of point ids
    for n_tiles, tile_geom in enumerate(self._patch.geometry):
      # get all the points close to this polygon - most will be its corners
      # but some will be corners of neighbouring tiles that lie on its edges
      near_ids = [i for i, v in enumerate(self.points)
                  if tile_geom.distance(v) <= tiling_utils.RESOLUTION]
      for id in near_ids:
        self.point_tiles[id].append(n_tiles)
        
      near_points = [self.points[i] for i in near_ids]
      tile_sides = tiling_utils.get_sides(tile_geom)
      tile_corner_ids, tile_vertex_ids = [], []
      for line_segment in tile_sides:
        ends = [geom.Point(p) for p in line_segment.coords]
        line_segment_ids = [self._get_closest_point(p) for p in ends] 
        tile_corner_ids.append(line_segment_ids[0])

        # now find any potential induced vertices to insert along the edge
        inserts = [idx for idx, p in zip(near_ids, near_points)
                   if idx not in line_segment_ids and
                   line_segment.distance(p) <= tiling_utils.RESOLUTION]
        # might (rarely) be more than one so sort them by distance along edge
        dist_along = [line_segment.line_locate_point(
          self.points[idx], normalized = True) for idx in inserts]
        inserts = [inserts[dist_along.index(p)] for p in sorted(dist_along)]
        # and insert them
        tile_vertex_ids.extend(line_segment_ids[:1] + inserts)
      # convert to tuples and add to the tile edge and vertex lists
      self.tile_corners.append(tuple(tile_corner_ids))
      self.tile_vertices.append(tile_vertex_ids)


  def get_tile_geometry(self, i) -> geom.Polygon:
    """Returns shapely geometry of tile with any tiling vertices that are
    'in line' along a tile side removed."""
    return geom.Polygon([self.points[j] for j in self.tile_corners[i]])


  def get_tile_geometries(self) -> gpd.GeoSeries:
    return gpd.GeoSeries([self.get_tile_geometry(i) 
                          for i, t in enumerate(self.tiles)
                          if len(self.tile_vertices[i]) > 2])
  
  
  def get_tile_topology(self, i) -> geom.Polygon:
    return geom.Polygon([self.points[j] for j in self.tile_vertices[i]])
    

  def get_tile_topologies(self) -> gpd.GeoSeries:
    return gpd.GeoSeries([self.get_tile_topology(i) 
                          for i, t in enumerate(self.tiles)])
  
  
  def plot(self):
    """Plots the topology.

    Args:
        show_original_tiles (bool, optional): _description_. Defaults to True.
        show_tile_centres (bool, optional): _description_. Defaults to True.
    """
    fig = pyplot.figure(figsize = (10, 10))
    ax = fig.add_axes(111)
    
    gpd.GeoSeries(self.get_tile_geometries()).plot(
      ax = ax, fc = "dodgerblue", alpha = 0.35, lw = 1, ec = "k")

    gpd.GeoSeries(self.get_tile_topologies()).buffer(-10, cap_style = 3).plot(
      ax = ax, fc = "grey", alpha = 0.5, lw = 0)
    
    for i, c in enumerate(self.tile_centres):
      ax.annotate(text = i, xy = (c.x, c.y), fontsize = 9, 
                  color = "b", ha = "center", va = "center")
        
    gpd.GeoSeries(self.points).plot(ax = ax, color = "r")
    for i, p in enumerate(self.points):
      ax.annotate(text = i, xy = (p.x, p.y), color = "r")

    pyplot.axis("off")
    return ax
  

  def get_dual_tiles(self) -> gpd.GeoDataFrame:
    """Alternative approach to creating the dual tiie unit for a tile unit.

    Args:
        t (TileUnit): _description_

    Returns:
        gpd.GeoDataFrame: _description_
    """
    dual_vertices = [(i, self.points[i]) 
                     for i, vps in enumerate(self.vertex_tiles) 
                     if len(vps) > 2]
    dual_faces = [geom.Polygon([self.tile_centres[i] for i in p])
      for p in [self.vertex_tiles[i] for i, v in dual_vertices]]
    dual_faces = [
      f for f in dual_faces 
      if affine.translate(self.tile_unit.prototile.geometry[0],
                          tiling_utils.RESOLUTION, tiling_utils.RESOLUTION) \
                            .intersects(f.centroid)]
        
    # TODO: label and select the faces to include so that they are tileable...
    return gpd.GeoSeries(dual_faces)
