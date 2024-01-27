#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
from collections import defaultdict
from typing import Union
from typing import Iterable
import copy
import pickle
import string

import geopandas as gpd
import shapely.geometry as geom
import matplotlib.pyplot as pyplot

from weavingspace import TileUnit
from weavingspace import Symmetries
import weavingspace.tiling_utils as tiling_utils

from time import perf_counter

LETTERS = string.ascii_letters.upper()
letters = string.ascii_letters.lower()


@dataclass
class Tile(object):
  ID: int = None
  corners: list["Vertex"] = None  # all the corners incl tiling vertices
  edges: list["Edge"] = None
  edges_cw: list[bool] = None
  vertex_labels: list[str] = None
  edge_labels: list[str] = None
  shape: geom.Polygon = None
  centre: geom.Point = None
  group: int = None
  

  def __init__(self, ID:int):
    # self.set_shape(tiling_utils.get_clean_polygon(shape))
    self.ID = ID
    self.corners = []
    self.edges = []
    self.edges_cw = []
    self.vertex_labels = []
    self.edge_labels = []


  def __repr__(self) -> str:
    return f"Tile {self.get_corner_IDS()}"


  def get_corner_IDS(self) -> list[int]:
    return [c.ID for c in self.corners]


  def get_edge_IDS(self) -> list[int]:
    return [e.ID for e in self.edges]


  def set_shape(self, shape):
    self.shape = shape
    self.centre = tiling_utils.incentre(
      tiling_utils.get_clean_polygon(self.shape))
  

  def get_current_shape(self):
    if self.shape is None:
      # this is used during initialisation of the Topology
      return geom.Polygon([c.point for c in self.corners])
    else:
      return self.shape


  def set_shape_from_corners(self):
    self.set_shape(
      geom.Polygon([v.point for v in self.corners]))


  def set_corners_from_edges(self, update_shape:bool = True):
    self.corners = []
    for e, cw in zip(self.edges, self.edges_cw):
      if cw:
        self.corners.extend(e.corners[:-1])
      else:
        self.corners.extend(reversed(e.corners[1:]))
    if update_shape:
      self.set_shape_from_corners()


  def insert_vertex_at(self, v:"Vertex", i:int, update_shape:bool = False):
    # This is not a generic vertex insertion method: it is only for use during 
    # Topology initialisation, and does not guarantee correct maintenance of 
    # all tile, edge and vertex relations in the general case!
    print(f"corners={self.get_corner_IDS()},v={v.ID},{i=}")
    self.corners = self.corners[:i] + [v] + self.corners[i:]
    print(f"edges={self.get_edge_IDS()}\ncw={self.edges_cw}")
    old_edge = self.edges[i - 1]
    old_edge_ID = old_edge.ID
    print(f"{old_edge=}")
    new_edges = old_edge.insert_vertex(v, self.corners[i - 1])
    if i >= len(self.edges_cw) or self.edges_cw[i]:
      self.edges = self.edges[:(i-1)] + new_edges + self.edges[i:]
      self.edges_cw = self.edges_cw[:(i-1)] + [True] * 2 + self.edges_cw[i:]
    else:
      self.edges = self.edges[:(i-1)] + new_edges + self.edges[i:]
      self.edges_cw = self.edges_cw[:(i-1)] + [True, True] + self.edges_cw[i:]
    print(f"corners={self.get_corner_IDS()}")
    print(f"edges={self.get_edge_IDS()}\ncw={self.edges_cw}\n")
    if update_shape:
      self.set_shape_from_corners()
    return old_edge_ID, new_edges


  def offset_corners(self, offset):
    self.shape = tiling_utils.offset_polygon_corners(self.shape, offset)
    self.corners = self.corners[offset:] + self.corners[:offset]
    self.edges = self.edges[offset:] + self.edges[:offset]
    self.edges_cw = self.edges_cw[offset:] + self.edges_cw[:offset]


  def get_edge_label(self, e):
    return self.edge_labels[self.get_edge_IDS().index(e.ID)]


  def get_corner_label(self, v):
    return self.edge_labels[self.get_corner_IDS().index(v.ID)]


  def get_vertex_label_positions(self) -> list[geom.Point]:
    # this convoluted method required because various buffer and parallel
    # offset options seem to return clean shapes with points removed!
    d = (self.shape.area ** 0.5) / 8
    corners = [c.point for c in self.corners]
    c = self.centre
    return [geom.LineString([p, c]).line_interpolate_point(d) for p in corners]
    

  def get_edge_label_positions(self) -> list[geom.Point]:
    d = (self.shape.area ** 0.5) / 8
    sides = [e.get_geometry(cw) for e, cw in zip(self.edges, self.edges_cw)]
    c = self.centre
    return [geom.LineString(
      [s.centroid, c]).line_interpolate_point(d) for s in sides]


@dataclass
class Vertex:
  point: geom.Point = None
  ID: int = None
  tiles: list["Tile"] = None
  neighbours: list[int] = None
  label: str = None
  on_boundary: bool = False
  is_tiling_vertex: bool = True   # TODO: fix this...
  

  def __init__(self, point:geom.Point, ID:int):
    self.point = point
    self.ID = ID
    self.tiles = []
    self.neighbours = []


  def __repr__(self) -> str:
    return f"Vertex ID {self.ID} at {self.point}"


  def add_tile(self, tile):
    if not tile in self.tiles:
      self.tiles.append(tile)
      

  def add_neighbour(self, vertex_id):
    if not vertex_id in self.neighbours:
      self.neighbours.append(vertex_id)


@dataclass
class Edge:
  ID: tuple[int] = None
  vertices: list["Vertex"] = None  # only (end) points which are tiling vertices
  corners: list["Vertex"] = None   # all the points
  left_tile: "Tile" = None
  right_tile: "Tile" = None
  label: str = None
  

  def __init__(self, vs:list[Vertex]):
    self.corners = vs
    self.vertices = [v for v in self.corners if v.is_tiling_vertex]
    self.ID = tuple(v.ID for v in self.vertices)


  def __repr__(self) -> str:
    return f"Edge ID {self.ID} {[c.ID for c in self.corners]}"


  def get_corner_IDs(self) -> list[int]:
    return [c.ID for c in self.corners]


  def get_vertex_IDs(self) -> list[int]:
    return [v.ID for v in self.vertices]


  def insert_vertex(self, v:"Vertex", predecessor:"Vertex") -> list["Edge"]:
    i = self.corners.index(predecessor)
    new_edge = Edge([v] + self.corners[(i+1):])
    new_edge.left_tile = self.left_tile
    new_edge.right_tile = self.right_tile
    self.corners = self.corners[:(i+1)] + [v]
    self.vertices = [self.vertices[0], v]
    self.ID = tuple(v.ID for v in self.vertices)
    return [self, new_edge]


  def get_geometry(self, forward = True):
    if forward:
      return geom.LineString([v.point for v in self.corners])
    else:
      return geom.LineString([v.point for v in self.corners[::-1]])


  def get_topology(self, forward = True):
    if forward:
      return geom.LineString([v.point for v in self.vertices])
    else:
      return geom.LineString([v.point for v in self.vertices[::-1]])


class Topology:
  """Class to represent topology of a Tileable.
  """
  tile_unit: TileUnit
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
    self._initialise_vertices()
    self._setup_edges()
    self._copy_base_tiles_to_patch()
    self._generate_dual()
    # self._trim_to_r1()
    self._label_elements()


  def __repr__(self) -> str:
    return f"Topology"


  def _initialise_vertices(self):
    self._initialise_points_into_tiles()
    self._setup_vertex_tile_relations()
    # self._setup_boundary_points()
    # self._assign_vertices()
    self._reorder_vertex_incident_tiles()


  def _label_elements(self):
    self._identify_distinct_tiles()
    self._label_tiles()
    self._label_vertices()
    self._label_edges()


  def _initialise_points_into_tiles(self) -> None:
    """Sets up the list of unique point locations in the tiling.
    """
    shapes = self.tile_unit.get_local_patch(r = 1, include_0 = True).geometry
    shapes = [tiling_utils.get_clean_polygon(s) for s in shapes]
    self.tiles = []
    n_tiles = 0
    self.points = {}
    n_points = 0
    for shape in shapes:
      tile = Tile(n_tiles)
      self.tiles.append(tile)
      n_tiles = n_tiles + 1
      tile.corners = []
      corners = tiling_utils.get_corners(shape, repeat_first = False)
      for c in corners:
        prev_vertex = None
        for p in reversed(self.points.values()):
          if c.distance(p.point) <= 2 * tiling_utils.RESOLUTION:
            prev_vertex = p
            break
        if prev_vertex is None:
          v = self.add_vertex(c)
          tile.corners.append(v)
          n_points = n_points + 1
        else:
          tile.corners.append(prev_vertex)


  def _setup_vertex_tile_relations(self) -> None:
    for tile in self.tiles:
      corners = []
      current_tile_corner_IDs = tile.get_corner_IDS()
      # unfortunately (for performance) we need the current shape (although
      # it has not yet been set) to check for vertices incident on the tile
      shape = tile.get_current_shape()
      new_points_on_boundary = [
        v for v in self.points.values()
        if v.ID not in current_tile_corner_IDs and
        v.point.distance(shape) <= 2 * tiling_utils.RESOLUTION]
      for c1, c2 in zip(tile.corners, tile.corners[1:] + tile.corners[:1]):
        if len(new_points_on_boundary) > 0:
          ls = geom.LineString([c1.point, c2.point])
          inserts = [v for v in new_points_on_boundary 
                     if v.point.distance(ls) <= 2 * tiling_utils.RESOLUTION]
          d_along = sorted([(ls.line_locate_point(v.point), v)
                            for v in inserts], key = lambda x: x[0])
          to_insert = [v for d, v in d_along]
          all_points = [c1] + to_insert + [c2]
        else:
          all_points = [c1, c2]
        corners.extend(all_points[:-1])
        for x1, x2 in zip(all_points[:-1], all_points[1:]):
          x1.add_tile(tile)
          x1.add_neighbour(x2.ID)
          x2.add_neighbour(x1.ID)
      tile.corners = corners
      tile.set_shape_from_corners()


  def _copy_base_tiles_to_patch(self):
    n_base = self.tile_unit.tiles.shape[0]
    n_r1 = len(self.tiles)
    for base_tile in self.tiles[:n_base]:
      for other_tile in self.tiles[base_tile.ID:n_r1:n_base]:
        self._map_reference_tile_onto_tile(base_tile, other_tile)


  def _map_reference_tile_onto_tile(self, ti:Tile, tj:Tile):
    to_add = len(ti.corners) - len(tj.corners)
    if to_add > 0:
      dxy = (tj.centre.x - ti.centre.x, tj.centre.y - ti.centre.y)
      print(f"{ti.ID=},{tj.ID=},{dxy=}")
      for i, tic in enumerate([c.point for c in ti.corners]):
        tjc = tj.corners[i % len(tj.corners)].point
        if abs((tjc.x - tic.x) - dxy[0]) > 10 * tiling_utils.RESOLUTION or \
           abs((tjc.y - tic.y) - dxy[1]) > 10 * tiling_utils.RESOLUTION:
          # add vertex to tj
          v = self.add_vertex(geom.Point(tic.x + dxy[0], tic.y + dxy[1]))
          v.is_tiling_vertex = True
          old_edge, new_edges = tj.insert_vertex_at(v, i)
          del self.edges[old_edge]
          for e in new_edges:
            self.edges[e.ID] = e
          to_add = to_add - 1
          if to_add == 0:
            break


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
      # v.is_tiling_vertex = \
      #   not v.on_boundary and len(v.tiles) > 2 \
      #   or v.on_boundary and len(v.neighbours) > 2
      v.is_tiling_vertex = True


  def _setup_edges(self) -> None:
    self.edges = {}
    n_edges = 0
    for tile in self.tiles:
      tile.edges = []
      tile.edges_cw = []
      vertices = [v for v in tile.corners if v.is_tiling_vertex]
      # note that through here finding ints in lists is much faster than
      # finding Vertex objects, hence we assemble and use lists of Vertex IDs
      if len(vertices) > 1:
        for v1, v2 in zip(vertices, vertices[1:] + vertices[:1]):
          if not (v1.on_boundary and v2.on_boundary):
            corner_IDs = tile.get_corner_IDS()
            idx1 = corner_IDs.index(v1.ID)
            idx2 = corner_IDs.index(v2.ID)
            if idx1 < idx2:
              corners = [c for c in corner_IDs[idx1:(idx2 + 1)]]
            else:
              corners = [c for c in corner_IDs[idx1:] + corner_IDs[:(idx2 + 1)]]
            ID = (corners[0], corners[-1])
            if not ID in self.edges:
              r_ID = ID[::-1]
              if r_ID in self.edges:
                e = self.edges[r_ID]
                e.left_tile = tile
                tile.edges.append(e)
                tile.edges_cw.append(False)
              else:
                e = self.add_edge([self.points[c] for c in corners])
                e.right_tile = tile
                tile.edges.append(e)
                tile.edges_cw.append(True)


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
    n_tiles = self.tile_unit.tiles.shape[0]
    for group in range(self.n_tile_groups):
      tiles = [t for t in self.tiles[:n_tiles] if t.group == group]
      s = Symmetries(self.example_tile_shapes[group])
      vlabels = list(s.get_unique_labels(offset = first_letter)["rotations"][0])
      elabels = self._get_edge_labels_from_vertex_labels(vlabels)
      for tile in tiles:
        tile.vertex_labels = vlabels
        tile.edge_labels = elabels
      first_letter = first_letter + len(set(vlabels))
    for tile in self.tiles[n_tiles:]:
      tile.vertex_labels = self.tiles[tile.ID % n_tiles].vertex_labels
      tile.edge_labels = self.tiles[tile.ID % n_tiles].edge_labels


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
    n_tiles = self.tile_unit.tiles.shape[0]
    # first label vertices in the core tiles
    vs = set()
    for t in self.tiles[:n_tiles]:
      vs = vs.union(t.get_corner_IDS())
    for vi in vs:
      v = self.points[vi]
      label = ""
      for tile in v.tiles:
        label = label + \
          tile.vertex_labels[tile.corners.index(v)]
      v.label = "".join(self._cyclic_sort_first(list(label)))
      uniques.add(v.label)
    for vi in vs:
      v = self.points[vi]
      v.label = LETTERS[list(uniques).index(v.label)]
    # now copy to corresponding vertices in the rest of tiling
    for ti, t in enumerate(self.tiles):
      if ti >= n_tiles:
        t0 = self.tiles[ti % n_tiles]
        for v0, v in zip(t0.corners, t.corners):
          if self.points[v.ID].label is None:
            self.points[v.ID].label = self.points[v0.ID].label


  def _label_edges(self) -> None:
    uniques = set()
    labelled = set()
    # first label the core edges from the central tile unit
    n_tiles = self.tile_unit.tiles.shape[0]
    for t in self.tiles[:n_tiles]:
      for e in t.edges:
        rt, lt = e.right_tile, e.left_tile
        e.label = "".join(sorted([rt.get_edge_label(e),
                                  lt.get_edge_label(e)]))
        uniques.add(e.label)
        labelled.add(e.ID)
    for ei in labelled:
      e = self.edges[ei]
      e.label = letters[list(uniques).index(e.label)]
    # now copy to corresponding edges in the rest of tiling
    for ti, t in enumerate(self.tiles):
      if ti >= n_tiles:
        t0 = self.tiles[ti % n_tiles]
        for e0, e in zip(t0.edges, t.edges):
          if e.label is None:
            e.label = e0.label


  def _reorder_vertex_incident_tiles(self) -> None:
    for v in self.points.values():
      p0 = v.point
      incident_tiles = v.tiles
      cw_order = tiling_utils.order_of_pts_cw_around_centre(
        [t.centre for t in incident_tiles], p0)
      v.tiles = [incident_tiles[i] for i in cw_order]


  def _cyclic_sort_first(self, lst:Iterable) -> Iterable:
    cycle = lst + lst
    cycles = [cycle[i:] + cycle[:i] for i in range(len(cycle))]
    return sorted(cycles)[0][:len(lst)]


  def add_vertex(self, pt:geom.Point) -> "Vertex":
    n = 0 if len(self.points) == 0 else max(self.points.keys()) + 1
    v = Vertex(pt, n)
    self.points[n] = v
    return v

  def add_edge(self, vs:list["Vertex"]) -> "Edge":
    e = Edge(vs)
    self.edges[e.ID] = e
    return e


  # def _trim_to_r1(self):
  #   # trim back tiles to the local patch radius 1
  #   r1_n = len(self.tile_unit.get_local_patch(include_0 = True).geometry)
  #   self.tiles = self.tiles[:r1_n]
  #   # compile lists of points and edges to retain
  #   pts, es = set(), set()
  #   for tile in self.tiles:
  #     pts = pts.union(tile.get_corner_IDS())
  #     es = es.union(tile.get_edge_IDS())
  #   # trim points back
  #   self.points = {k: v for k, v in self.points.items() if v.ID in pts}
  #   for v in self.points.values():
  #     v.tiles = [t for t in v.tiles if t.ID < r1_n]
  #   # trim edges back
  #   self.edges = {k: e for k, e in self.edges.items() if k in es}
  #   for e in self.edges.values():
  #     if not e.right_tile.ID < r1_n:
  #       e.right_tile = None
  #     if not e.left_tile.ID < r1_n:
  #       e.left_tile = None
  #   self.dual_tiles = {k: v for k, v in self.dual_tiles.items() if k in pts}


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
    return gpd.GeoSeries([e.get_topology().parallel_offset(offset)
                          for e in self.edges.values()])


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
            geom.Polygon([t.centre for t in v.tiles])


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
    
    extent = gpd.GeoSeries([t.shape for t in self.tiles]).total_bounds
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
                      ha = "center", va = "center", color = "b")
    
    if show_tile_edge_labels:
      for tile in self.tiles:
        for posn, label in zip(tile.get_edge_label_positions(),
                               tile.edge_labels):
          ax.annotate(label, xy = (posn.x, posn.y), color = "b",
                      ha = "center", va = "center", fontsize = 8)
    
    if show_edge_labels:
      if show_edges:
        edges = self.get_edge_geoms(dist / 200 if offset_edges else 0)
        edges.plot(ax = ax, color = "forestgreen", ls = "dashed", lw = 1)
      else:
        edges = [e.get_geometry() for e in self.edges.values()]
      for l, e in zip(edges, self.edges.values()):
        c = l.centroid
        ax.annotate(e.label, xy = (c.x, c.y), ha = "center", va = "center",
          fontsize = 10 if show_edge_labels else 7, color = "k")
      
    if show_dual_tiles:
      gpd.GeoSeries(self.dual_tiles).buffer(
        -dist / 400, join_style = 2, cap_style = 3).plot(
          ax = ax, fc = "red", alpha = 0.1)
    
    pyplot.axis("off")
    return ax


  def transform_edges(self, new_topology:bool, apply_to_tiles:bool,
                      label:str, type:str, **kwargs) -> "Topology":
    topo = pickle.loads(pickle.dumps(self)) if new_topology else self
    for e in topo.edges.values():
      if e.label in label:
        if type == "zigzag":
          topo.zigzag_edge(e, **kwargs)
    if apply_to_tiles:
      for t in topo.tiles:
        t.set_shape_from_corners()
    n = topo.tile_unit.tiles.shape[0]
    topo.tile_unit.tiles.geometry = gpd.GeoSeries(
      [tiling_utils.get_clean_polygon(topo.tiles[i].shape) for i in range(n)])
    topo.tile_unit.setup_regularised_prototile_from_tiles()
    return topo


  def zigzag_edge(self, edge:"Edge", 
                  n:int = 2, h:float = 0.5, smoothness:int = 0):
    p0 = edge.vertices[0].point
    p1 = edge.vertices[1].point
    ls = tiling_utils.zigzag_between_points(
      p0, p1, n, h, smoothness)
    # remove current corners
    self.points = {k: v for k, v in self.points.items()
                   if not k in edge.get_corner_IDs()[1:-1]}
    new_corners = []
    for xy in ls.coords:
      v = self.add_vertex(geom.Point(xy))
      new_corners.append(v)
    edge.corners = edge.vertices[:1] + new_corners + edge.vertices[-1:]
    if not edge.right_tile is None:
      edge.right_tile.set_corners_from_edges(False)
    if not edge.left_tile is None:
      edge.left_tile.set_corners_from_edges(False)
