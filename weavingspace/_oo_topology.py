#!/usr/bin/env python
# coding: utf-8
from __future__ import annotations

from dataclasses import dataclass
from typing import Union
from typing import Self

import geopandas as gpd
import shapely.geometry as geom

from weavingspace import TileUnit
from weavingspace import TileShape
import weavingspace.tiling_utils as tiling_utils


@dataclass
class Vertex:
  """Represents tiling vertices relative to other elements."""
  topology: Topology = None
  point: geom.Point = None
  incident_tiles: list[Tile] = None
  incident_edges: list[Edge] = None
  neighbours: list[Vertex] = None
  
  def __init__(self, topology: Topology, point:geom.Point):
    self.topology = topology
    self.point = point
    self.incident_tiles = []
    self.incident_edges = []
    self.neighbours = []
    pass
  
  def __eq__(self, other:Vertex) -> bool:
    return self.point.distance(other.point) <= 2 * tiling_utils.RESOLUTION
  
  def is_vertex(self) -> bool:
    return len(self.neighbours) > 2
    pass
  
  def add_incident_tile(self, tile:Tile) -> None:
    self.incident_tiles.append(tile)
    if len(self.incident_tiles) > 1:
       cw_order = tiling_utils.order_of_pts_cw_around_centre(
         [t.centre for t in self.incident_tiles], tile.centre)
       self.incident_tiles = [self.incident_tiles[i] for i in cw_order]
  
  def add_edge_to(self, other:Self) -> None:
    pass
  
  def add_edge_from(self, other:Self) -> None:
    pass

class Corner:
  point: geom.Point = None
  previous_corner: Union[Corner, Vertex] = None
  next_corner: Union[Corner, Vertex] = None
  
@dataclass
class Tile:
  """Represents tiling tiles relative to other elements.""" 
  topology: Topology = None
  corners: list[Corner] = None
  sides: list[Side] = None
  centre: geom.Point = None

  def __init__(self, topology, polygon:geom.Polygon):
    self.topology = topology
    self.vertices = []
    self.edges = []
    self.sides = []
    self.centre = tiling_utils.incentre(polygon)
    pass
  
  def get_adjacents(self) -> list[Self]:
    return [e.other_side(self) for e in self.edges]
  
  def get_neighbours(self) -> list[Self]:
    pass


@dataclass
class Edge:
  """Represents tiling edges relative to other elements."""
  topology: Topology = None
  vertices: list[Vertex] = None
  left: Tile = None
  right: Tile = None
  
  def __init__(self, topology:Topology):
    self.topology = topology
    self.vertices = []
    self.left = None
    self.right = None
    pass
  
  def get_start(self) -> Vertex:
    return self.vertices[0]
  
  def get_end(self) -> Vertex:
    return self.vertices[-1]
  
  def other_side(self, tile:Tile) -> Tile:
    if self.left == tile:
      return self.right
    elif self.right == tile:
      return self.left
    else:
      return None


@dataclass
class Topology:
  vertices:tuple[Vertex] = None
  tiles:tuple[Tile] = None
  edges:tuple[Edge] = None
  
  def __init__(self, unit:TileUnit) -> None:
    r = 1 if unit.base_shape == TileShape.HEXAGON else 2
    patch = unit.get_local_patch(r = 1, include_0 = True)
    self._initialise_vertices(patch)
    pass
  
  def _initialise_vertices(self, patch:gpd.GeoDataFrame) -> None:
    self.vertices = []
    self.tiles = []
    for shape in patch.geometry:
      tile = Tile(self, shape)
      self.tiles.append(tile)
      for corner in tiling_utils.get_corners(shape, repeat_first = False):
        v = self.get_vertex_at(corner)
        v.add_incident_tile(tile)
        tile.add_vertex(v)
  
  def get_vertex_at(self, pt:geom.Point) -> Union[None, Vertex]:
    new_v = Vertex(self, pt)
    old_vs = [old_v for old_v in self.vertices if old_v == new_v]
    if len(old_vs) == 0:
      self.vertices.append(new_v)
      return new_v
    else:
      return old_vs[0]
