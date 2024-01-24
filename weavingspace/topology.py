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
class Topology:
  """Class to represent topology of a Tileable.
  """

  tile_unit: TileUnit = None
  _patch: gpd.GeoDataFrame = None

  points: dict[int, geom.Point] = None
  points_on_boundary: dict[int, int] = None
  point_tiles: dict[int, list[int]] = None
  point_edges: dict[int, int] = None
  point_labels: dict[int, str] = None

  is_vertex: dict[int, bool] = None

  edges: dict[int, tuple[int]] = None
  edge_lefts: dict[int, int] = None
  edge_rights: dict[int, int] = None
  edge_labels: dict[int, str] = None

  tiles: list[geom.Polygon] = None
  tile_ids: list[str] = None
  distinct_tiles: dict[int, list[int]] = None
  example_tiles: list[geom.Polygon] = None
  tile_group: list[int] = None
  tile_centres: list[geom.Point] = None
  tile_adjacents: dict[int, list[int]] = None
  tile_corners: tuple[list[int]] = None
  tile_points: list[list[int]] = None
  tile_vertex_labels: list[list[str]] = None
  tile_edge_labels: list[list[str]] = None
  
  dual_tiles: list[geom.Polygon] = None

  def __init__(self, unit: TileUnit):
    """Builds topology for a supplied TileUnit. WeaveUnits will also work
    but results may be a bit strange, especially if aspect < 1.

    Args:
        unit (TileUnit): the weavingspace.TileUnit whose topology we want.
    """
    self.tile_unit = unit # keep this for reference
    self._patch = self.tile_unit.get_local_patch(r = 2, include_0 = True)
    # save the original geometries so we have them on hand, but make sure
    # they are clean first!
    self.tiles = [tiling_utils.get_clean_polygon(g) 
                  for g in self._patch.geometry]
    self.tile_ids = self._patch.tile_id
    self.tile_centres = [tiling_utils.incentre(t) for t in self.tiles]
    # build basic topology of tiles, vertices, and edges
    self._setup_unique_points_and_remake_tiles()
    self._setup_point_and_tile_relations()
    self._setup_boundary_points()
    self._assign_vertices()
    self._setup_edges()
    self._reorder_vertex_incident_tiles()
    self._order_vertex_incident_edges()
    # identify groups of distinct tiles
    self._identify_distinct_tiles()
    self._label_tiles()
    self._label_vertices()
    self._label_edges()
    self._generate_dual()
    self._filter_back_to_r1()


  def _setup_unique_points_and_remake_tiles(self) -> None:
    """Sets up the list of unique point locations in the tiling.
    """
    self.points = {}
    self.tile_corners = []
    n = 0
    for i, tile in enumerate(self.tiles):
      tile_corners = []
      corners = tiling_utils.get_corners(tile, repeat_first = False)
      for corner in corners:
        # only add if not close to an existing point
        if not any ([corner.distance(point) <= 2 * tiling_utils.RESOLUTION 
                     for point in self.points.values()]):
          self.points[n] = corner
          tile_corners.append(n)
          n = n + 1
        else:
          tile_corners.append(self._get_closest_point(corner, return_id = True))
      # replace the tile with this version
      self.tiles[i] = geom.Polygon([self.points[i] for i in tile_corners])
      self.tile_corners.append(tile_corners)


  def _setup_point_and_tile_relations(self) -> None:
    self.point_tiles = defaultdict(set)
    self.point_neighbours = defaultdict(set)
    self.tile_points = []
    for tile_id, tile in enumerate(self.tiles):
      tile_vertices = []
      corners = self.tile_corners[tile_id]
      new_points_on_boundary = [
        (i, pt) for i, pt in self.points.items()
        if pt.distance(tile) <= 2 * tiling_utils.RESOLUTION 
        and i not in corners]
      for i, (c1, c2) in enumerate(zip(corners, corners[1:] + corners[:1])):
        ls = geom.LineString([self.points[c1], self.points[c2]])
        inserts = [(j, pt) for j, pt in new_points_on_boundary
                  if pt.distance(ls) <= 2 * tiling_utils.RESOLUTION]
        x_along = sorted([(ls.line_locate_point(pt), j) for j, pt in inserts])
        to_insert = [j for d, j in x_along]
        all_points = [c1] + to_insert + [c2]
        tile_vertices.extend(all_points[:-1])
        for x1, x2 in zip(all_points[:-1], all_points[1:]):
          self.point_tiles[x1].add(tile_id)
          self.point_neighbours[x1].add(x2)
          self.point_neighbours[x2].add(x1)
      self.tile_points.append(tile_vertices)
    self.tiles = [geom.Polygon(self.points[v] for v in vs)
                  for vs in self.tile_points]


  def _setup_boundary_points(self) -> list[int]:
    """Sets up the points_on_boundary list. 
    """
    boundary = [
      geom.Point(xy)
      for xy in gpd.GeoSeries(self.tiles).unary_union.exterior.coords][:-1]
    self.points_on_boundary = [
      i for i, pt in self.points.items()
      if any([pt.distance(b) <= 2 * tiling_utils.RESOLUTION for b in boundary])]


  def _assign_vertices(self) -> None:
    self.is_vertex = {}
    for i, pt in enumerate(self.points):
      self.is_vertex[i] = \
        i not in self.points_on_boundary and len(self.point_tiles[i]) > 2 \
        or i in self.points_on_boundary and len(self.point_neighbours[i]) > 2


  def _setup_edges(self) -> None:
    self.edges = {}
    self.edge_lefts = {}
    self.edge_rights = {}
    self.tile_edges = []
    n_edges = 0
    for tile_id, all_tile_points in enumerate(self.tile_points):
      vertices = [v for v in all_tile_points if self.is_vertex[v]]
      tile_edges = []
      if len(vertices) > 1:
        for v1, v2 in zip(vertices, vertices[1:] + vertices[:1]):
          if not (v1 in self.points_on_boundary and v2 in self.points_on_boundary):
            idx1 = all_tile_points.index(v1)
            idx2 = all_tile_points.index(v2)
            if idx1 < idx2:
              corners = all_tile_points[idx1:(idx2 + 1)]
            else:
              corners = all_tile_points[idx1:] + all_tile_points[:(idx2 + 1)]
            if not corners in self.edges.values():
              r_corners = list(reversed(corners))
              if r_corners in self.edges.values():
                id = list(self.edges.keys())[
                  list(self.edges.values()).index(r_corners)]
                self.edge_lefts[id] = tile_id
                tile_edges.append(-id - 1)
              else:
                self.edges[n_edges] = corners
                self.edge_rights[n_edges] = tile_id
                tile_edges.append(n_edges)
                n_edges = n_edges + 1
      self.tile_edges.append(tile_edges)


  def _identify_distinct_tiles(self) -> None:
    self.distinct_tiles = defaultdict(set)
    self.tile_group = []
    offsets = []
    self.example_tiles = []
    n_groups = 0
    n_tiles = self.tile_unit.tiles.geometry.shape[0]
    # r1_n = self.tile_unit.get_local_patch(r = 1, include_0 = True).shape[0]
    for tile_id in range(n_tiles):
      this_tile = self.tiles[tile_id]
      s = Symmetries(this_tile)
      transforms = [s.get_matching_transforms(other) 
                    for other in self.example_tiles]
      matches = [not t is None for t in transforms]
      if any(matches):
        match = matches.index(True)
        self.distinct_tiles[match].add(tile_id)
        self.tile_group.append(match)
        offset = transforms[match]["offset"]
        offsets.append(offset)
      else:
        self.example_tiles.append(this_tile)
        self.distinct_tiles[n_groups].add(tile_id)
        offsets.append(0)
        self.tile_group.append(n_groups)
        n_groups = n_groups + 1
    # now copy assignment from the central TileUnit to everything else
    for i in range(len(self.tiles)):
      self.distinct_tiles[self.tile_group[i % n_tiles]].add(i)
      self._offset_tile_corners(i, offsets[i % n_tiles])
    self.distinct_tiles = {k: list(v) for k, v in self.distinct_tiles.items()}


  def _offset_tile_corners(self, idx, offset):
    self.tiles[idx] = tiling_utils.offset_polygon_corners(
      self.tiles[idx], -offset)
    self.tile_corners[idx] = \
      self.tile_corners[idx][offset:] + self.tile_corners[idx][:offset]
    self.tile_points[idx] = \
      self.tile_points[idx][offset:] + self.tile_points[idx][:offset]
    self.tile_edges[idx] = \
      self.tile_edges[idx][offset:] + self.tile_edges[idx][:offset]


  def _label_tiles(self) -> None:
    tile_vertex_labels = {}
    tile_edge_labels = {}
    first_letter = 0
    # find distinct tile types
    for group, tiles in self.distinct_tiles.items():
      s = Symmetries(self.example_tiles[group])
      vlabels = s.get_unique_labels(offset = first_letter)["rotations"][0]
      elabels = self._get_edge_labels_from_vertex_labels(vlabels)
      for tile_id in tiles:
        tile_vertex_labels[tile_id] = vlabels
        tile_edge_labels[tile_id] = elabels
        # my_label = label
        # points = self.tile_points[tile_id]
        # corners = self.tile_corners[tile_id]
        # for i, point in enumerate(points):
        #   if not point in corners:
        #     my_label = my_label[:i] + LETTERS[-group - 1] + my_label[i:]
        # tile_labels[tile_id] = my_label
      first_letter = first_letter + len(set(vlabels))
    self.tile_vertex_labels = \
      [tile_vertex_labels[i] for i, t in enumerate(self.tiles)]
    self.tile_edge_labels = \
      [tile_edge_labels[i] for i, t in enumerate(self.tiles)]


  def _get_edge_labels_from_vertex_labels(self, vertex_labels:str) -> list[str]:
    vlabels = list(vertex_labels)
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
    point_labels = {}
    for pt_id, incident_tiles in self.point_tiles.items():
      label = ""
      for tile_id in incident_tiles:
        label = label + \
          self.tile_vertex_labels[tile_id][self.tile_points[tile_id].index(pt_id)]
      # label1 = self._cyclic_sort_first(label)
      # label2 = self._cyclic_sort_first(label[::-1])
      # point_labels[pt_id] = min(label1, label2)
      # point_labels[pt_id] = self._cyclic_sort_first(label)
      point_labels[pt_id] = min(list(label))
    uniques = list(set(point_labels.values()))
    self.point_labels = {k: LETTERS[uniques.index(v)] 
                         for k, v in point_labels.items()}


  def _label_edges(self) -> None:
    self.edge_labels = {}
    elabels = {}
    for i, es in enumerate(self.tile_edges):
      labels = self.tile_edge_labels[i]
      for j, e in enumerate(es):
        elabels[e] = labels[j]
    for e, label in elabels.items():
      if -e - 1 in elabels:
        self.edge_labels[e] = "".join(sorted([label, elabels[-e - 1]]))
      else:
        self.edge_labels[e] = label
    u_labels = sorted(set(self.edge_labels.values()))
    lookup = dict(zip(u_labels, letters[:len(u_labels)]))
    self.edge_labels = {e: lookup[l] for e, l in self.edge_labels.items()}


  def _reorder_vertex_incident_tiles(self) -> None:
    cw_ordered_v_incident_tiles = {}
    for v, incidents in self.point_tiles.items():
      p0 = self.points[v]
      incident_tiles = list(incidents)
      cw_order = tiling_utils.order_of_pts_cw_around_centre(
        [self.tile_centres[t] for t in incident_tiles], p0)
      cw_ordered_v_incident_tiles[v] = [incident_tiles[i] for i in cw_order]
    self.point_tiles = cw_ordered_v_incident_tiles


  def _order_vertex_incident_edges(self) -> None:
    self.point_edges = {}
    for i, p0 in self.points.items():
      incident_edges = [k for k, e in self.edges.items()
                        if e[0] == i or e[-1] == i]
      nearest_points = [self.edges[e][1] if i == self.edges[e][0] 
                        else self.edges[e][-2] for e in incident_edges]
      cw_order = tiling_utils.order_of_pts_cw_around_centre(
        [self.points[p] for p in nearest_points], p0)
      cw_ordered_edges = [incident_edges[i] for i in cw_order]
      self.point_edges[i] = cw_ordered_edges


  def remake_tile_from_edges(self, i:int) -> geom.Polygon:
    edges = self.tile_edges[i]
    vertices = []
    for e in edges:
      if e < 0:
        points = reversed([self.points[i] for i in self.edges[-e - 1][1:]])
        vertices.extend(points)
      else:
        points = [self.points[i] for i in self.edges[e][:-1]]
        vertices.extend(points)
    return geom.Polygon(vertices)


  def get_tile_vertices(self, i: int) -> list[int]:
    edges = self.tile_edges[i]
    vertices = []
    for e in edges:
      if e < 0:
        vertices.append(self.edges[-e - 1][-1])
      else:
        vertices.append(self.edges[e][0])
    return vertices


  def _cyclic_sort_first(self, lst:Iterable) -> Iterable:
    cycle = lst + lst
    cycles = [cycle[i:] + cycle[:i] for i in range(len(cycle))]
    return sorted(cycles)[0][:len(lst)]


  def _filter_back_to_r1(self):
    r1_n = len(self.tile_unit.get_local_patch(include_0 = True).geometry)
    self.tiles = self.tiles[:r1_n]
    self.tile_centres = self.tile_centres[:r1_n]
    self.tile_corners = self.tile_corners[:r1_n]
    self.tile_points = self.tile_points[:r1_n]
    self.tile_vertex_labels = self.tile_vertex_labels[:r1_n]
    self.tile_edge_labels = self.tile_edge_labels[:r1_n]
    self.tile_edges = self.tile_edges[:r1_n]
    self.tile_group = self.tile_group[:r1_n]
    distinct_tiles = {}
    for i, g in self.distinct_tiles.items():
      distinct_tiles[i] = [t for t in g if t < r1_n]
    self.distinct_tiles = distinct_tiles  

    pts = set()
    for t in self.tile_points:
      pts = pts.union(t)
    self.points = {k: v for k, v in self.points.items() if k in pts}
    self.is_vertex = {k: v for k, v in self.is_vertex.items() if k in pts}
    self.points_on_boundary = None
    self.point_neighbours = \
      {k: v for k, v in self.point_neighbours.items() if k in pts}
    for p, N in self.point_neighbours.items():
      self.point_neighbours[p] = [p for p in N if p in pts]
    self.point_tiles = {k: v for k, v in self.point_tiles.items() if k in pts}
    for p, T in self.point_tiles.items():
      self.point_tiles[p] = [t for t in T if t < r1_n]
    self.point_labels = {k: v for k, v in self.point_labels.items() if k in pts}
    self.point_edges = {k: v for k, v in self.point_edges.items() if k in pts}
    unique_labels = list(set(self.point_labels.values()))
    self.point_labels = {k: LETTERS[unique_labels.index(v)] 
                         for k, v in self.point_labels.items()}

    self.dual_tiles = {k: v for k, v in self.dual_tiles.items() if k in pts}

    self.edges = {k: vs for k, vs in self.edges.items()
                  if all([v in pts for v in vs])}
    for p, N in self.point_edges.items():
      self.point_edges[p] = [e for e in N if e in self.edges]
    self.edge_lefts = {k: v for k, v in self.edge_lefts.items() 
                       if k in self.edges and v < r1_n}
    self.edge_rights = {k: v for k, v in self.edge_rights.items() 
                        if k in self.edges and v < r1_n}
    self.edge_labels = {k: v for k, v in self.edge_labels.items() 
                        if k in self.edges}
    unique_labels = list(set(self.edge_labels.values()))
    self.edge_labels = {k: letters[unique_labels.index(v)] 
                         for k, v in self.edge_labels.items()}


  def _get_closest_point(self, point:Union[int, geom.Point], 
                         return_id:bool = True) -> Union[int, geom.Point]:
    """Finds the closest vertex to that supplied from the Topology's vertex 
    list.

    Args:
      v (Union[int, geom.Point]): vertex to find closest other vertex from, 
        supplied either as an integer index or a Point geometry.
      return_id (bool, optional): if True returns the closest vertex id, if  
        False the point. Defaults to True.

    Returns:
      Union[int, geom.Point]: the closest point in self.points as an index
        or Point geometry.
    """
    this_pt = self.points[point] if isinstance(point, int) else point
    distances_to_points = sorted([(this_pt.distance(other_pt), i) 
                                  for i, other_pt in self.points.items()])
    if return_id:
      return distances_to_points[0][1]
    else:
      return self.points[distances_to_points[0][1]]


  def get_geometries_induced_by_tiles(self, tile_ids, type = "vertices"):
    if type == "vertices":
      vs = set()
      for t in [self.tile_points[i] for i in tile_ids]:
        vs = vs.union(t)
      return {v: self.points[v] for v in vs}
    if type == "edges":
      es = set()
      for t in [self.tile_edges[i] for i in tile_ids]:
        es = es.union([e if e >= 0 else -e -1 for e in t])
      return {e: geom.LineString([self.points[v] for v in self.edges[e]]) 
              for e in es}


  def get_tile_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries(self.tiles)


  def get_tile_centre_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries(self.tile_centres)
  
  
  def get_point_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries(self.points)


  def get_vertex_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries([p for i, p in self.points.items()
                          if self.is_vertex[i]])


  def get_edge_geoms(self, offset:float = 0.0) -> gpd.GeoSeries:
    geoms = gpd.GeoSeries([
      geom.LineString([self.points[vs[0]], 
                       self.points[vs[-1]]]).parallel_offset(offset) 
      for vs in self.edges.values()])
    return geoms


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
    for v, is_vertex in self.is_vertex.items():
      if is_vertex:
        points = self.point_tiles[v]
        if len(points) > 2:
          self.dual_tiles[v] = \
            geom.Polygon([self.tile_centres[i] for i in points])


  def plot(self, show_original_tiles: bool = True,  
           show_tile_centres: bool = False,
           show_tile_vertex_labels: bool = False, 
           show_tile_edge_labels: bool = False, 
           show_vertex_ids: bool = False,
           show_vertex_labels: bool = True,
           show_edges: bool = True,
           offset_edges: bool = True,
           show_edge_labels:bool = False,
           show_edge_ids: bool = False,
           show_edge_corners: bool = False,
           show_dual_tiles: bool = False) -> pyplot.Axes:
  
    fig = pyplot.figure(figsize = (10, 10))
    ax = fig.add_axes(111)
    
    extent = self._patch.total_bounds
    dist = max([extent[2] - extent[0], extent[3] - extent[1]])

    if show_original_tiles:
      self.get_tile_geoms().plot(
        ax = ax, fc = "dodgerblue", ec = "#333366", alpha = 0.2, lw = 0.75)

    if show_tile_centres == "id":
      # for i, pt in enumerate(self.get_tile_centre_geoms()):
      for i, pt in zip(self.tile_vertex_labels, self.get_tile_centre_geoms()):
        ax.annotate(i, xy = (pt.x, pt.y), color = "b", fontsize = 9,
                    ha = "center", va = "center")
    elif show_tile_centres:
      for i, pt in enumerate(self.get_tile_centre_geoms()):
        ax.annotate(i, xy = (pt.x, pt.y), color = "b", fontsize = 9,
                    ha = "center", va = "center")

    if show_vertex_labels:
      if show_vertex_ids:
        for i, pt in self.points.items():
          ax.annotate(i, xy = (pt.x, pt.y), color = "r", fontsize = 9,
                      ha = "center", va = "center")
      else:
        for i, lbl in self.point_labels.items():
          pt = self.points[i]
          ax.annotate(lbl, xy = (pt.x, pt.y), color = "r", fontsize = 9,
                      ha = "center", va = "center",
                      bbox = dict(boxstyle="circle", lw=0, fc="#ffffff40"))

    if show_tile_vertex_labels:
      for t, c, l in zip(self.tile_points, self.tile_centres, self.tile_vertex_labels):
        for corner, label in zip(t, l):
          posn = geom.LineString(
            [self.points[corner], c]).line_interpolate_point(0.2, True)
          ax.annotate(label, xy = (posn.x, posn.y), fontsize = 7, 
                      ha = "center", va = "center")

    if show_tile_edge_labels:
      for i, c in enumerate(self.tile_centres):
        a = self.tiles[i].area
        d = (a ** 0.5) / 8
        pts = [self.points[v] for v in self.tile_points[i]]
        tile_edges = [geom.LineString([p0, p1])
                      for p0, p1 in zip(pts, pts[1:] + pts[:1])]
        midpoints = [l.centroid for l in tile_edges]
        labelpts = [geom.LineString([m, c]).line_interpolate_point(d) 
                    for m in midpoints]
        for label, labelpoint in zip(self.tile_edge_labels[i], labelpts):
          ax.annotate(label, xy = (labelpoint.x, labelpoint.y),
                      ha = "center", va = "center", fontsize = 7)
    
    # if show_vertices:
    #   self.get_vertex_geoms().plot(
    #     ax = ax, ec = "red", lw = 0.5, fc = "#00000000",
    #     marker = "o", markersize = 100)
    
    if show_edges:
      edges = self.get_edge_geoms(dist / 200 if offset_edges else 0)
      edges.plot(ax = ax, color = "forestgreen", ls = "dashed", lw = 1)
      if show_edge_labels:
        for i, e in self.edges.items():
          c = edges[i].centroid
          ax.annotate(self.edge_labels[i], xy = (c.x, c.y),
                      ha = "center", va = "center", fontsize = 9)
      elif show_edge_ids or show_edge_corners:
        labels = self.get_edge_labels(show_edge_ids, show_edge_corners)
        for label, ls in zip(labels, edges): 
          ax.annotate(label, xy = (ls.centroid.x, ls.centroid.y),
                      color = "k", ha = "center", va = "top", fontsize = 7)
    
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