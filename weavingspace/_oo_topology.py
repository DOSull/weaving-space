#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
from collections import defaultdict
from typing import Union
from typing import Iterable
import copy
import string

import geopandas as gpd
import shapely.geometry as geom
import matplotlib.pyplot as pyplot

from weavingspace import TileUnit
from weavingspace import Symmetries
import weavingspace.tiling_utils as tiling_utils

LETTERS = string.ascii_letters.upper()
letters = string.ascii_letters.lower()


@dataclass
class Tile(object):
  shape: geom.Polygon
  ID: int
  group: int
  centre: geom.Point
  corners: list[int]
  points: list[int]
  vertex_labels: list[str]
  edge_labels: list[str]
  
  def __init__(self, shape:geom.Polygon, ID:int):
    self.set_shape(tiling_utils.get_clean_polygon(shape))
    self.ID = ID

  def set_shape(self, shape):
    self.shape = shape
    self.centre = tiling_utils.incentre(self.shape)
    
  def offset_corners(self, offset):
    self.shape = tiling_utils.offset_polygon_corners(self.shape, offset)
    self.corners = self.corners[offset:] + self.corners[:offset]
    self.points = self.points[offset:] + self.points[:offset]
    self.edges = self.edges[offset:] + self.edges[:offset]
    
  def get_vertex_label_positions(self) -> list[geom.Point]:
    d = (self.shape.area ** 0.5) / 10
    pts = tiling_utils.get_corners(self.shape, repeat_first = False)
    return [
      geom.LineString([p, self.centre]).line_interpolate_point(d)
      for p in pts]
    
  def get_edge_label_positions(self) -> list[geom.Point]:
    d = (self.shape.area ** 0.5) / 10
    pts = tiling_utils.get_sides(self.shape)
    return [
      geom.LineString([p.centroid, self.centre]).line_interpolate_point(d)
      for p in pts]


@dataclass
class Vertex:
  point: geom.Point
  ID: int
  tiles: list[int]
  neighbours: list[int]
  edges: list[int]
  label: str = ""
  on_boundary: bool = False
  is_tiling_vertex: bool = False
  
  def __init__(self, point:geom.Point, ID:int):
    self.point = point
    self.ID = ID
    self.tiles = []
    self.neighbours = []
    self.edges = []
    
  def add_tile(self, tile_id):
    if not tile_id in self.tiles:
      self.tiles.append(tile_id)
      
  def add_neighbour(self, vertex_id):
    if not vertex_id in self.neighbours:
      self.neighbours.append(vertex_id)


@dataclass
class Edge:
  ID: int
  vertices: list[int]
  line: geom.LineString
  left_polygon: int
  right_polygon: int
  label: str = ""
  
  def __init__(self, vs:list[Vertex], ID:int):
    self.vertices = [v.ID for v in vs]
    self.line = geom.LineString([v.point for v in vs])
    self.ID = ID
  

class Topology:
  """Class to represent topology of a Tileable.
  """
  tile_unit: TileUnit
  _patch: gpd.GeoDataFrame
  tiles: list[Tile]
  example_tile_shapes: list[geom.Polygon]
  points: dict[int, Vertex]
  edges: dict[int, Edge]
  dual_tiles: list[geom.Polygon]
  n_tile_groups: int = 0

  def __init__(self, unit: TileUnit):
    """Builds topology for a supplied TileUnit. WeaveUnits will also work
    but results may be a bit strange, especially if aspect < 1.

    Args:
        unit (TileUnit): the weavingspace.TileUnit whose topology we want.
    """
    self.tile_unit = unit # keep this for reference
    self._patch = self.tile_unit.get_local_patch(r = 2, include_0 = True)
    self.tiles = [Tile(g, i) for i, g in enumerate(self._patch.geometry)]

    # build basic topology of tiles, vertices, and edges
    self._setup_unique_points_and_remake_tiles()
    self._setup_point_and_tile_relations()
    self._setup_boundary_points()
    self._assign_vertices()
    self._setup_edges()
    self._reorder_vertex_incident_tiles()
    self._order_vertex_incident_edges()
    # # # identify groups of distinct tiles
    self._identify_distinct_tiles()
    self._label_tiles()
    self._label_vertices()
    self._label_edges()
    self._generate_dual()
    self._trim_to_r1()

  def get_points(self) -> list[Vertex]:
    return list(self.points.values())

  def _setup_unique_points_and_remake_tiles(self) -> None:
    """Sets up the list of unique point locations in the tiling.
    """
    self.points = {}
    n = 0
    for tile in self.tiles:
      tile.corners = []
      corners = tiling_utils.get_corners(tile.shape, repeat_first = False)
      for corner in corners:
        # only add if not close to an existing point
        if not any ([corner.distance(p.point) <= 2 * tiling_utils.RESOLUTION 
                     for p in self.points.values()]):
          v = Vertex(corner, n)
          self.points[n] = v
          tile.corners.append(v.ID)
          n = n + 1
        else:
          tile.corners.append(self._get_closest_vertex(corner).ID)
      # replace the tile with this version
      tile.set_shape(geom.Polygon([self.points[i].point for i in tile.corners]))

  def _get_closest_vertex(self, pt:geom.Point) -> Vertex:
    """Finds the closest vertex to that supplied from the Topology's vertex 
    list.

    Args:
      v (Union[int, geom.Point]): vertex to find closest other vertex from, 
        supplied either as an integer index or a Point geometry.
      return_id (bool, optional): if True returns the closest vertex id, if  
        False the point. Defaults to True.

    Returns:
      Vertex: the closest Vertex in self.points.
    """
    distances_to_points = sorted([(pt.distance(other.point), other) 
                                  for other in self.points.values()],
                                  key = lambda x: x[0])
    return distances_to_points[0][1]

  def _setup_point_and_tile_relations(self) -> None:
    for tile_id, tile in enumerate(self.tiles):
      tile.points = []
      new_points_on_boundary = [
        (i, v) for i, v in self.points.items()
        if v.point.distance(tile.shape) <= 2 * tiling_utils.RESOLUTION 
        and i not in tile.corners]
      for i, (c1, c2) in enumerate(zip(tile.corners, 
                                       tile.corners[1:] + tile.corners[:1])):
        ls = geom.LineString([self.points[c1].point, self.points[c2].point])
        inserts = [(j, pt.point) for j, pt in new_points_on_boundary
                  if pt.point.distance(ls) <= 2 * tiling_utils.RESOLUTION]
        x_along = sorted([(ls.line_locate_point(pt), j) for j, pt in inserts])
        to_insert = [j for d, j in x_along]
        all_points = [c1] + to_insert + [c2]
        tile.points.extend(all_points[:-1])
        for x1, x2 in zip(all_points[:-1], all_points[1:]):
          self.points[x1].add_tile(tile_id)
          self.points[x1].add_neighbour(x2)
          self.points[x2].add_neighbour(x1)
      tile.set_shape(geom.Polygon([self.points[i].point for i in tile.points]))

  def _setup_boundary_points(self) -> list[int]:
    """Sets up the points_on_boundary list. 
    """
    patch = tiling_utils.safe_union(
      gpd.GeoSeries([tile.shape for tile in self.tiles]), as_polygon = True)
    patch = tiling_utils.get_clean_polygon(patch)
    bounds = tiling_utils.get_corners(patch, repeat_first = False)
    for v in self.points.values():
      v.on_boundary = any(
        [v.point.distance(p) <= 2 * tiling_utils.RESOLUTION for p in bounds])

  def _assign_vertices(self) -> None:
    for v in self.points.values():
      v.is_tiling_vertex = \
        not v.on_boundary and len(v.tiles) > 2 \
        or v.on_boundary and len(v.neighbours) > 2

  def _setup_edges(self) -> None:
    self.edges = {}
    n_edges = 0
    for tile_id, tile in enumerate(self.tiles):
      vertices = [self.points[v] for v in tile.points]
      tile.edges = []
      if len(vertices) > 1:
        for v1, v2 in zip(vertices, vertices[1:] + vertices[:1]):
          if not (v1.on_boundary and v2.on_boundary):
            idx1 = tile.points.index(v1.ID)
            idx2 = tile.points.index(v2.ID)
            if idx1 < idx2:
              corners = tile.points[idx1:(idx2 + 1)]
            else:
              corners = tile.points[idx1:] + tile.points[:(idx2 + 1)]
            edges = [e.vertices for e in self.edges.values()]
            if not corners in edges:
              r_corners = list(reversed(corners))
              if r_corners in edges:
                edge_id = list(self.edges.keys())[edges.index(r_corners)]
                self.edges[edge_id].left_polygon = tile_id
                tile.edges.append(-edge_id - 1)
              else:
                self.edges[n_edges] = \
                  Edge([self.points[v] for v in corners], n_edges)
                self.edges[n_edges].right_polygon = tile_id
                tile.edges.append(n_edges)
                n_edges = n_edges + 1


  def _identify_distinct_tiles(self) -> None:
    self.example_tile_shapes = []
    offsets = []
    n_tiles = self.tile_unit.tiles.geometry.shape[0]
    for tile_id in range(n_tiles):
      this_tile = self.tiles[tile_id]
      s = Symmetries(this_tile.shape)
      transforms = [s.get_matching_transforms(other) 
                    for other in self.example_tile_shapes]
      matches = [not t is None for t in transforms]
      if any(matches):
        match = matches.index(True)
        this_tile.group = match
        offset = transforms[match]["offset"]
        offsets.append(offset)
      else:
        self.example_tile_shapes.append(this_tile.shape)
        offsets.append(0)
        this_tile.group = self.n_tile_groups
        self.n_tile_groups = self.n_tile_groups + 1
    # now copy assignment from the central TileUnit to everything else
    for i, tile in enumerate(self.tiles):
      tile.group = self.tiles[i % n_tiles].group
      tile.offset_corners(offsets[i % n_tiles])


  def _label_tiles(self) -> None:
    first_letter = 0
    for group in range(self.n_tile_groups):
      tiles = [t for t in self.tiles if t.group == group]
      s = Symmetries(self.example_tile_shapes[group])
      vlabels = list(s.get_unique_labels(offset = first_letter)["rotations"][0])
      elabels = self._get_edge_labels_from_vertex_labels(vlabels)
      for tile in tiles:
        tile.vertex_labels = vlabels
        tile.edge_labels = elabels
      first_letter = first_letter + len(set(vlabels))

  def _get_edge_labels_from_vertex_labels(self, 
                                          vlabels:list[str]) -> list[str]:
    edge_labels = [a + b for a, b in zip(vlabels, vlabels[1:] + vlabels[:1])]
    letter = LETTERS.index(min(list(vlabels)))
    elabels = {}
    for e in edge_labels:
      if not e in elabels:
        if e[::-1] in elabels:
          elabels[e] = elabels[e[::-1]]
        else:
          elabels[e] = letters[letter]
          letter = letter + 1
      # if not e in elabels:
      #   if e[::-1] in elabels:
      #     label = elabels[e[::-1]]
      #     elabels[e] = label + "-"
      #     elabels[e[::-1]] = label + "+"
      #   else:
      #     elabels[e] = letters[letter]
      #     letter = letter + 1
    return [elabels[l] for l in edge_labels]

  def _label_vertices(self) -> None:
    uniques = set()
    for v in self.points.values():
      label = ""
      for tile_id in v.tiles:
        tile = self.tiles[tile_id]
        label = label + \
          tile.vertex_labels[tile.points.index(v.ID)]
      v.label = min(list(label))
      uniques.add(v.label)
    for v in self.points.values():
      v.label = LETTERS[list(uniques).index(v.label)]

  def _label_edges(self) -> None:
    elabels = {}
    for tile in self.tiles:
      labels = tile.edge_labels
      for j, e in enumerate(tile.edges):
        elabels[e] = labels[j]
    for e, label in {i: lbl for i, lbl in elabels.items() if i >= 0}.items():
      if -e - 1 in elabels:
        self.edges[e].label = "".join(sorted([label, elabels[-e - 1]]))
      else:
        self.edges[e].label = label
    u_labels = sorted(set([e.label for e in self.edges.values()]))
    lookup = dict(zip(u_labels, letters[:len(u_labels)]))
    for e in self.edges.values():
      e.label = lookup[e.label]

  def _reorder_vertex_incident_tiles(self) -> None:
    for v in self.points.values():
      p0 = v.point
      incident_tiles = v.tiles
      cw_order = tiling_utils.order_of_pts_cw_around_centre(
        [self.tiles[t].centre for t in incident_tiles], p0)
      v.tiles = [incident_tiles[i] for i in cw_order]

  def _order_vertex_incident_edges(self) -> None:
    for i, v in self.points.items():
      incident_edges = [(i, e) for i, e in self.edges.items()
                        if e.vertices[0] == i or e.vertices[-1] == i]
      nearest_points = [e.vertices[1] if i == e.vertices[0] else e.vertices[-2]
                        for i, e in incident_edges]
      cw_order = tiling_utils.order_of_pts_cw_around_centre(
        [self.points[p].point for p in nearest_points], v.point)
      v.edges = [incident_edges[i][0] for i in cw_order]

  def _cyclic_sort_first(self, lst:Iterable) -> Iterable:
    cycle = lst + lst
    cycles = [cycle[i:] + cycle[:i] for i in range(len(cycle))]
    return sorted(cycles)[0][:len(lst)]

  def _trim_to_r1(self):
    r1_n = len(self.tile_unit.get_local_patch(include_0 = True).geometry)
    self.tiles = self.tiles[:r1_n]
    pts = set()
    for tile in self.tiles:
      pts = pts.union(tile.points)
    self.points = {k: v for k, v in self.points.items() if v.ID in pts}
    unique_labels = list(set([v.label for v in self.points.values()]))
    for v in self.points.values():
      v.label = LETTERS[unique_labels.index(v.label)]
    es = set()
    for tile in self.tiles:
      es = es.union(tile.edges)
    self.edges = {k: e for k, e in self.edges.items() if k in es}
    unique_labels = list(set([e.label for e in self.edges.values()]))
    for e in self.edges.values():
      e.label = letters[unique_labels.index(e.label)]
    self.dual_tiles = {k: v for k, v in self.dual_tiles.items() if k in pts}

  def get_tile_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries([t.shape for t in self.tiles])

  def get_tile_centre_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries([t.centre for t in self.tiles])

  def get_point_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries([v.point for v in self.points.values()])

  def get_vertex_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries([v.point for v in self.points.values()
                          if v.is_tiling_vertex])

  def get_edge_geoms(self, offset:float = 0.0) -> gpd.GeoSeries:
    return gpd.GeoSeries([e.line.parallel_offset(offset)
                          for e in self.edges.values()])

  def get_edge_labels(self, 
                       include_ids:bool, 
                       include_corners:bool) -> list[str]:
    labels = []
    for i, corners in self.edges.items():
      parts = []
      if include_ids: parts.append(f"{i}") 
      if include_corners:
        parts.append(u"\u00b7".join([str(c) for c in corners]))
      labels.append(":".join(parts))
    return labels

  def _generate_dual(self) -> gpd.GeoSeries:
    """Alternative approach to creating the dual tiie unit for a tile unit.

    Returns:
        gpd.GeoDataFrame: _description_
    """
    self.dual_tiles = {}
    for v in self.points.values():
      if v.is_tiling_vertex:
        if len(v.tiles) > 2:
          self.dual_tiles[v.ID] = \
            geom.Polygon([self.tiles[i].centre for i in v.tiles])


  def plot(self, show_original_tiles: bool = True,  
           show_tile_centres: bool = False,
           show_tile_vertex_labels: bool = False, 
           show_tile_edge_labels: bool = False, 
           show_vertex_ids: bool = False,
           show_vertex_labels: bool = True,
           show_edges: bool = True,
           offset_edges: bool = True,
           show_edge_labels:bool = False,
           show_dual_tiles: bool = False) -> pyplot.Axes:
    
    fig = pyplot.figure(figsize = (10, 10))
    ax = fig.add_axes(111)
    
    extent = self._patch.total_bounds
    dist = max([extent[2] - extent[0], extent[3] - extent[1]])
    
    if show_original_tiles:
      self.get_tile_geoms().plot(
        ax = ax, fc = "dodgerblue", ec = "#333366", alpha = 0.2, lw = 0.75)
    
    if show_tile_centres:
      for i, tile in enumerate(self.tiles):
        ax.annotate(i, xy = (tile.centre.x, tile.centre.y), color = "b", 
                    fontsize = 10, ha = "center", va = "center")
    
    if show_vertex_labels:
      for v in self.points.values():
          ax.annotate(v.ID if show_vertex_ids else v.label, 
                      xy = (v.point.x, v.point.y), color = "r", fontsize = 10,
                      ha = "center", va = "center",
                      bbox = dict(boxstyle="circle", lw=0, fc="#ffffff40"))
    
    if show_tile_vertex_labels:
      for tile in self.tiles:
        for posn, label in zip(tile.get_vertex_label_positions(), 
                               tile.vertex_labels):
          ax.annotate(label, xy = (posn.x, posn.y), fontsize = 8, 
                      ha = "center", va = "center")
    
    if show_tile_edge_labels:
      for tile in self.tiles:
        for posn, label in zip(tile.get_edge_label_positions(),
                               tile.edge_labels):
          ax.annotate(label, xy = (posn.x, posn.y),
                      ha = "center", va = "center", fontsize = 8)
    
    if show_edges:
      edges = self.get_edge_geoms(dist / 200 if offset_edges else 0)
      edges.plot(ax = ax, color = "forestgreen", ls = "dashed", lw = 1)
    
    if show_edge_labels:
      for e in self.edges.values():
        c = e.line.centroid
        ax.annotate(e.label, xy = (c.x, c.y), ha = "center", va = "center",
          fontsize = 10 if show_edge_labels else 7, color = "k")
    
    if show_dual_tiles:
      gpd.GeoSeries(self.dual_tiles).buffer(
        -dist / 400, join_style = 2, cap_style = 3).plot(
          ax = ax, fc = "red", alpha = 0.1)
    
    pyplot.axis("off")
    return ax


  def get_new_tile_unit(self, edge_transforms):
    edges_to_transform = list(edge_transforms.keys())
    tile_edges = self.tile_edges[:self.tile_unit.tiles.geometry.shape[0]]
    new_tiles = []

    for t in tile_edges:
      new_t = []
      for edge in t:
        if edge >= 0:
          pts = self.edges[edge]
          e_label = self.edge_labels[edge]
        else:
          pts = (self.edges[-edge - 1])[::-1]
          e_label = self.edge_labels[-edge - 1]
        # e_label = f"{self.point_labels[pts[0]]}{self.point_labels[pts[-1]]}"
        # e_label_r = e_label[::-1]
        if not e_label in edges_to_transform:
          new_t.extend([self.points[i] for i in pts[:-1]])
        else:
          if e_label in edge_transforms:
            transform = edge_transforms[e_label]
            transform_type = transform["type"]
            transform_to_apply = {k: v for k, v in transform.items()
                                  if k != "type"}
          else:
            continue
          if transform_type == "zigzag":
            p0 = self.points[pts[0]]
            p1 = self.points[pts[1]]
            ls = tiling_utils.zigzag_between_points(
              p0, p1, **transform_to_apply)
            new_t.extend([geom.Point(xy) for xy in ls.coords][:-1])
          elif transform_type == "stretch":
            ls = tiling_utils.stretch_between_points(
              p0, p1, **transform_to_apply)
            if len(new_t) == 0:
              new_t = [geom.Point(xy) for xy in ls.coords]
            else:
              new_t = new_t[:-1] + [geom.Point(xy) for xy in ls.coords]
      new_tiles.append(geom.Polygon(new_t))

    new_tile_unit = copy.deepcopy(self.tile_unit)
    new_tile_unit.tiles.geometry = gpd.GeoSeries(new_tiles)
    new_tile_unit.setup_regularised_prototile_from_tiles()
    return new_tile_unit