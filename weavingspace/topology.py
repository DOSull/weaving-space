#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
from collections import defaultdict
from typing import Union
from typing import Iterable
import string

import geopandas as gpd
import shapely.geometry as geom
# import shapely.affinity as affine
import shapely.ops

import matplotlib.pyplot as pyplot

from weavingspace import TileUnit
from weavingspace import Symmetries
import weavingspace.tiling_utils as tiling_utils

LETTERS = string.ascii_letters.upper()

@dataclass
class Topology:
  """Class to represent topology of a Tileable.
  """

  tile_unit: TileUnit = None
  _patch: gpd.GeoDataFrame = None

  points: dict[int, geom.Point] = None
  points_on_boundary: dict[int, int] = None
  point_tiles: dict[int, list[int]] = None
  point_neighbours: dict[int, list[int]] = None
  point_labels: dict[int, str] = None

  is_vertex: dict[int, bool] = None

  edges: list[tuple[int]] = None
  edge_lefts: dict[int, int] = None
  edge_rights: dict[int, int] = None

  tiles: list[geom.Polygon] = None
  distinct_tiles: dict[int, list[int]] = None
  example_tiles: list[geom.Polygon] = None
  tile_group: list[int] = None
  tile_centres: list[geom.Point] = None
  tile_corners: tuple[list[int]] = None
  tile_points: list[list[int]] = None
  tile_labels: list[list[str]] = None

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
    self.tile_centres = [tiling_utils.incentre(t) for t in self.tiles]
    self._setup_unique_points_and_rebuild_tiles()
    self._identify_distinct_tiles()
    self._label_tile_corners()
    self._setup_boundary_points()
    self._setup_point_and_tile_relations()
    self._assign_vertices()
    self._reorder_vertex_incident_tiles()
    self._label_vertices()
    self._setup_edges()
    self._filter_back_to_r1()


  def _setup_unique_points_and_rebuild_tiles(self) -> None:
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
  
  
  def _identify_distinct_tiles(self) -> None:
    self.distinct_tiles = defaultdict(list)
    self.tile_group = []
    self.example_tiles = []
    n_groups = 0
    for tile_id, this_tile in enumerate(self.tiles):
      this_tile = self.tiles[tile_id]
      s = Symmetries(this_tile)
      transforms = [s.get_matching_transforms(other) 
                    for other in self.example_tiles]
      matches = [not t is None for t in transforms]
      if any(matches):
        match = matches.index(True)
        self.distinct_tiles[match].append(tile_id)
        self.tile_group.append(match)
        this_tile = s._shift_polygon_corners(this_tile, 
                                             -transforms[match]["offset"])
      else:
        self.example_tiles.append(this_tile)
        self.distinct_tiles[n_groups].append(tile_id)
        self.tile_group.append(n_groups)
        n_groups = n_groups + 1
    self.distinct_tiles = dict(self.distinct_tiles)

      
  def _label_tile_corners(self) -> None:
    tile_labels = {}
    first_letter = 0
    # find distinct tile types
    for group, tiles in self.distinct_tiles.items():
      labels = Symmetries(self.example_tiles[group]).get_unique_labels(
        offset = first_letter)["rotations"][0]
      for tile_id in tiles:
        tile_labels[tile_id] = labels
      first_letter = first_letter + len(set(labels))
    self.tile_labels = [tile_labels[i] for i, t in enumerate(self.tiles)]


  def _label_vertices(self) -> None:
    point_labels = {}
    for pt_id, incident_tiles in self.point_tiles.items():
      label = ""
      for tile_id in incident_tiles:
        label = label + \
          self.tile_labels[tile_id][self.tile_points[tile_id].index(pt_id)]
      point_labels[pt_id] = self._cyclic_sort_first(label)
    uniques = list(set(point_labels.values()))
    self.point_labels = {k: LETTERS[uniques.index(v)] 
                         for k, v in point_labels.items()}


  def _setup_boundary_points(self) -> list[int]:
    """Sets up the points_on_boundary list. 
    """
    boundary = [
      geom.Point(xy)
      for xy in gpd.GeoSeries(self.tiles).unary_union.exterior.coords][:-1]
    self.points_on_boundary = [
      i for i, pt in self.points.items()
      if any([pt.distance(b) <= 2 * tiling_utils.RESOLUTION for b in boundary])]


  def _setup_point_and_tile_relations(self) -> None:
    self.point_tiles = defaultdict(set)
    self.point_neighbours = defaultdict(set)
    tile_points = {}
    for group, tiles in self.distinct_tiles.items():
      for tile_id in tiles:
        tile = self.tiles[tile_id]
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
          label = self.tile_labels[tile_id]
          self.tile_labels[tile_id] = \
            label[:i+1] + LETTERS[-group - 1] * len(to_insert) + label[i+1:]
        tile_points[tile_id] = tile_vertices
    self.tile_points = []
    for i in range(len(tile_points)):
      self.tile_points.append(tile_points[i])
      
    
  def _assign_vertices(self) -> None:
    self.is_vertex = {}
    for i, pt in enumerate(self.points):
      self.is_vertex[i] = \
        i not in self.points_on_boundary and len(self.point_tiles[i]) > 2 \
        or i in self.points_on_boundary and len(self.point_neighbours[i]) > 2
    
    
  def _reorder_vertex_incident_tiles(self) -> None:
    cw_ordered_v_incident_tiles = {}
    for v, incidents in self.point_tiles.items():
      p0 = self.points[v]
      incident_tiles = list(incidents)
      cw_order = tiling_utils.order_of_pts_cw_around_centre(
        [self.tile_centres[t] for t in incident_tiles], p0)
      cw_ordered_v_incident_tiles[v] = [incident_tiles[i] for i in cw_order]
    self.point_tiles = cw_ordered_v_incident_tiles


  def _setup_edges(self) -> None:
    self.edges = []
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
            if not corners in self.edges:
              r_corners = list(reversed(corners))
              if r_corners in self.edges:
                id = self.edges.index(r_corners)
                self.edge_lefts[id] = tile_id
                tile_edges.append(-id - 1)
              else:
                self.edges.append(corners)
                self.edge_rights[n_edges] = tile_id
                tile_edges.append(n_edges)
                n_edges = n_edges + 1
      self.tile_edges.append(tile_edges)
      
  
  def get_tile(self, i:int) -> geom.Polygon:
    edges = self.tile_edges[i]
    vertices = []
    for e in edges:
      if e < 0:
        vertices.append(self.points[self.edges[-e - 1][-1]])
      else:
        vertices.append(self.points[self.edges[e][0]])
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
    self.tile_labels = self.tile_labels[:r1_n]

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
    unique_labels = list(set(self.point_labels.values()))
    self.point_labels = {k: LETTERS[unique_labels.index(v)] 
                         for k, v in self.point_labels.items()}

    self.edges = [e for e in self.edges if all([v in pts for v in e])]
    self.edge_lefts = {k: v for k, v in self.edge_lefts.items() 
                       if k in self.edges and v < r1_n}
    self.edge_rights = {k: v for k, v in self.edge_rights.items() 
                        if k in self.edges and v < r1_n}


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
      for vs in self.edges])
    return geoms
    
  
  def get_edge_labels(self, 
                       include_ids:bool, 
                       include_corners:bool) -> list[str]:
    labels = []
    for i, corners in enumerate(self.edges):
      parts = []
      if include_ids: parts.append(f"{i}") 
      if include_corners:
        parts.append(u"\u00b7".join([str(c) for c in corners]))
      labels.append(":".join(parts))
    return labels


  def get_dual_tile_geoms(self) -> gpd.GeoSeries:
    """Alternative approach to creating the dual tiie unit for a tile unit.

    Returns:
        gpd.GeoDataFrame: _description_
    """
    dual_faces = []
    for v, is_vertex in self.is_vertex.items():
      if is_vertex:
        points = self.point_tiles[v]
        if len(points) > 2:
          dual_faces.append(
            geom.Polygon([self.tile_centres[i] for i in points]))
    # TODO: label and select which faces to include so that they are tileable
    return gpd.GeoSeries(dual_faces)


  def plot(self, show_original_tiles: bool = True,  
           show_tile_centres: bool = True,
           show_tile_vertex_labels: bool = False, 
          #  show_vertices: bool = False,
           show_vertex_ids: bool = False,
           show_vertex_labels: bool = False,
           show_edges: bool = True,
           offset_edges: bool = True,
           show_edge_ids: bool = False,
           show_edge_corners: bool = False,
           show_dual_tiles: bool = False) -> pyplot.Axes:
  
    fig = pyplot.figure(figsize = (10, 10))
    ax = fig.add_axes(111)

    if show_original_tiles:
      self.get_tile_geoms().plot(
        ax = ax, fc = "dodgerblue", ec = "#333399", alpha = 0.25, lw = 0.5)

    if show_tile_centres:
      for i, pt in enumerate(self.get_tile_centre_geoms()):
        ax.annotate(i, xy = (pt.x, pt.y), color = "b", fontsize = 9,
                    ha = "center", va = "center")

    if show_tile_centres == "id":
      # for i, pt in enumerate(self.get_tile_centre_geoms()):
      for i, pt in zip(self.tile_labels, self.get_tile_centre_geoms()):
        ax.annotate(i, xy = (pt.x, pt.y), color = "b", fontsize = 9,
                    ha = "center", va = "center")

    if show_vertex_labels:
      # self.get_point_geoms().plot(ax = ax, color = "r", markersize = 5)
      if show_vertex_ids:
        for i, pt in enumerate(self.get_point_geoms()):
          ax.annotate(i, xy = (pt.x, pt.y), color = "r", fontsize = 9,
                      ha = "center", va = "center")
      else:
        for i, lbl in self.point_labels.items():
          pt = self.points[i]
          ax.annotate(lbl, xy = (pt.x, pt.y), color = "r", fontsize = 9,
                      ha = "center", va = "center")

    if show_tile_vertex_labels:
      for t, c, l in zip(self.tile_points, self.tile_centres, self.tile_labels):
        for corner, label in zip(t, l):
          posn = geom.LineString(
            [self.points[corner], c]).line_interpolate_point(0.2, True)
          ax.annotate(label, xy = (posn.x, posn.y), fontsize = 7, 
                      ha = "center", va = "center")
    
    # if show_vertices:
    #   self.get_vertex_geoms().plot(
    #     ax = ax, ec = "red", lw = 0.5, fc = "#00000000",
    #     marker = "o", markersize = 100)
    
    if show_edges:
      edges = self.get_edge_geoms(
        self.tile_unit.spacing / 40 if offset_edges else 0)
      edges.plot(ax = ax, color = "forestgreen", ls = "dashed", lw = 1)
      if show_edge_ids or show_edge_corners:
        labels = self.get_edge_labels(show_edge_ids, show_edge_corners)
        for label, ls in zip(labels, edges): 
          ax.annotate(label, xy = (ls.centroid.x, ls.centroid.y),
                      color = "k", ha = "center", va = "top", fontsize = 7)
    
    if show_dual_tiles:
      self.get_dual_tile_geoms().buffer(
        -self.tile_unit.spacing / 100, join_style = 2, cap_style = 3).plot(
          ax = ax, fc = "red", alpha = 0.15)
    
    pyplot.axis("off")
    return ax
