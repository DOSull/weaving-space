#!/usr/bin/env python
# coding: utf-8

"""Together the `Topology`, `Tile`, `Vertex`, `Edge` and `weavingspace.symmetry.
Symmetry` and `weavingspace.symmetry.Symmetries` classes enable extraction of 
the topological structure of periodic `weavingspace.tileable.Tileable` objects 
so that modification of equivalent tiles can be carried out while retaining 
tileability. It is important to note that these are not fully generalised 
classes and methods, that is, the Topology object that is supported is not a 
permanent 'backing' data structure for our Tileable objects. While it might 
become that in time, as at Feb 2024 it is not such a data structure. Instead 
usage is

    tile = TileUnit(...)
    topology = Topology(tile)
    topology = topology.transform_*(...)
    new_tile = topology.tile_unit

Topology plot function is necessary for a user to be able to see what they are
doing, because how edges and vertices in a tiling are labelled under tile
equivalences is an essential step in the process.

Note also that these classes do not accurately represent the distinctions made
in the mathematical literature between tiling vertices and tile corners, or 
between tiling edges and tile sides.
"""
from collections import defaultdict
from typing import Union
from typing import Iterable
import pickle
import string

import numpy as np
import geopandas as gpd
import shapely.geometry as geom
import shapely.affinity as affine
import matplotlib.pyplot as pyplot

from weavingspace import Tileable
from weavingspace import Symmetries
import weavingspace.tiling_utils as tiling_utils

from time import perf_counter

ALPHABET = string.ascii_letters.upper()
alphabet = string.ascii_letters.lower()

class Tile(object):
  """Class to capture and manipulate essential features of polygons in a tiling.
  """
  ID: int
  """integer ID number which indexes the Tile in the containing Topology tiles 
  list."""
  corners: list["Vertex"]
  """list of Vertex objects. This includes all corners of the original polygon 
  and any tiling vertices induced by (for example) a the corner of an adjacent 
  tile lying halfway along an edge of the original polygon on which this tile 
  is based. Vertex objects are stored in strictly clockwise sequence."""
  edges: list["Edge"]
  """list of Edge objects that together compose the tile boundary."""
  edges_CW: list[bool]
  """list of Edge direction. Edges are stored only once in a Topology so some 
  edges are in clockwise order and others  are in counter-clockwise order. 
  These boolean flags are True if the corresponding Edge is clockwise, False if 
  counter-clockwise."""
  label: str
  """tile_id label from the tileable source"""
  vertex_labels: list[str]
  """list of (upper case) letter labels of the tile corners (i.e. all corners, 
  not only tiling vertices)."""
  edge_labels: list[str]
  """list of (lower case) letter labels of the tile edges (tiling edges, not 
  tile sides)."""
  shape: geom.Polygon = None
  """the tile geometry (which may include some redundant points along sides 
  where neighbouring tiles induce a tiling vertex). So for example a rectangle 
  might have additional points along its sides:
      
        +---+-------+
        |   |   2   |
        | 1 A---B---E---+
        |   |   |   4   |
        +---C 3 D-------+
            |   |
            +---+
      
  In the above Tile 1 has additional point A, 2 has B and 3 has C and D induced 
  by the corners of neighbouring tiles."""
  centre: geom.Point = None
  """a point centre for the Tile (determined by weavingspace.tiling_utils.
  incentre)."""
  group: int = None
  """the tile shape group of this tile in its containing Topology."""
  equivalence_class: int = None
  """the tile equivalence class of this tile its containing Topology"""
  
  def __init__(self, ID:int):
    """Class constructor.

    Args:
      ID (int): Tile ID which indexes it in the containing Topology tiles list.
    """
    # self.set_shape(tiling_utils.get_clean_polygon(shape))
    self.ID = ID
    self.corners = []
    self.edges = []
    self.edges_CW = []
    self.vertex_labels = []
    self.edge_labels = []

  def __str__(self) -> str:
    """Returns string representation of the Tile.

    Returns:
      str: string including Tile ID, list of corner vertex IDs and list of
        edge IDs.
    """
    return f"Tile {self.ID} Corners: {self.get_corner_IDs()} Edges: {self.get_edge_IDs()}"
  
  def __repr__(self) -> str:
    return str(self)

  def get_corner_IDs(self) -> list[int]:
    """Convenience method to return list of corner IDs (not Vertex objects).

    Returns:
      list[int]: list of integer IDs of tile corners.
    """
    return [c.ID for c in self.corners]

  def get_edge_IDs(self) -> list[tuple[int]]:
    """Convenience method to return list of edge IDs (not Edge objects).

    Returns:
      list[tuple[int]]: list of 2-element tuples of the start and end Vertex IDs
        of each edge.
    """
    return [e.ID for e in self.edges]

  def set_shape_from_corners(self):
    """Sets the shape attribute based on the current list of corners, and sets 
    the associated tile centre.
    """
    self.shape = geom.Polygon([c.point for c in self.corners])
    self.centre = tiling_utils.incentre(
      tiling_utils.get_clean_polygon(self.shape))

  def set_corners_from_edges(self, update_shape:bool = True):
    """Sets the corners attribute from the edges attribute. Typically called 
    after modification of topology edges. Optionally the shape attribute is NOT
    updated, which may save time when multiple changes to the edges of a tile
    are in process (i.e., only update the shape after all changes are complete).

    Args:
      update_shape (bool, optional): if True the shape attribute will be 
        updated, otherwise not. Defaults to True.
    """
    self.corners = []
    for e, cw in zip(self.edges, self.edges_CW):
      if cw: # clockwise to extend by all but the first corner
        self.corners.extend(e.corners[:-1])
      else: # counter-clockwise so extend in reverse
        self.corners.extend(e.corners[1:][::-1])
    if update_shape:
      self.set_shape_from_corners()

  def set_edge_directions(self):
    """Sets up the edges_CW attribute by inspection of the edges list.
    It is (frankly!) hard to keep track of the correct sequence of CW/CCW order
    of edges as new ones are created or old ones merged. This method inspects
    the 'tail-head' relations between consecutive edges to set these flags 
    correctly.

    The test is simply to check if the 'tail' Vertex ID in each edge appears
    in the ID tuple of the following edge, i.e. if successive edge 
    IDs are (0, 1) (2, 1) or (0, 1) (1, 2), then edge (0, 1) is in clockwise
    direction, but if we have (0, 1) (2, 3) then it is not.
    """
    edge_IDs = self.get_edge_IDs()
    self.edges_CW = [e1[-1] in e2 for e1, e2 in 
                     zip(edge_IDs, edge_IDs[1:] + edge_IDs[:1])]

  def insert_vertex_at(self, v:"Vertex", i:int, 
                       update_shape:bool = False) -> tuple:
    """Method to insert the supplied Vertex into tile at index position i, 
    optionally updating the shape attribute. Both corners and edges attributes
    are updated, and the old edge IDs for removal and the new edge itself are
    returned to the calling context (the containing Topology) for update of its
    edges collection.

    This is NOT a generic vertex insertion method: it is only for use during 
    Topology initialisation, and does not guarantee correct maintenance of 
    all tile, edge and vertex relations in the general case---at any rate it
    has not been tested for this!

    Args:
        v (Vertex): the Vertex to insert.
        i (int): index position in current corners after which to insert
          supplied Vertex.
        update_shape (bool, optional): if True shape attribute is updated. 
          Defaults to False.

    Returns:
      tuple: the (tuple) ID of the old edge which should be deleted, and 
        the new Edges arising from insertion of this Vertex.
    """
    self.corners = self.corners[:i] + [v] + self.corners[i:]
    old_edge = self.edges[i - 1]
    # store current ID of the affected edge for return to calling context
    old_edge_ID = old_edge.ID
    new_edges = old_edge.insert_vertex(v, self.corners[i - 1])
    self.edges = self.edges[:(i-1)] + new_edges + self.edges[i:]
    self.set_edge_directions()
    if update_shape:
      self.set_shape_from_corners()
    return old_edge_ID, new_edges

  def merge_edges_at_vertex(self, v:"Vertex") -> tuple:
    """Method to merge the edges that meet at the supplied Vertex. It is 
    assumed that only two tiles are impacted this one, and its neighbour across
    the Edge on which v lies. Both are updated. For this reason the work is
    delegated to get_updated_edges_from_merge which is run on both affected
    tiles, but only determines the edges to remove and the new edge to be added
    once. See that method for details.
    
    Args:
      v (Vertex): Vertex at which to merge Edges. This should currently be an
        end
    
    Returns:
      tuple: 2 item list of the edge IDs to be removed and a new Edge object to
        be added by the calling context (i.e. the containing Topology).
    """
    to_remove, new_edge = self.get_updated_edges_from_merge(v)
    if len(v.tiles) > 1:
      v.tiles[1].get_updated_edges_from_merge(v, new_edge)
    return to_remove, new_edge

  def get_updated_edges_from_merge(self, v:"Vertex", new_edge:"Edge" = None):
    """Updates the edges and edges_CW attributes based on insertion of 
    the supplied Vertex. If new_edge is supplied then the neighbour tile at
    v has already created the needed new Edge and this Edge is the one that
    will be 'slotted in' at the appropriate spot in the edges list.

    The edges_CW is also updated to maintain correct directions of the edges.
    The corners attribute is unaffected by these changes.

    Args:
      v (Vertex): existing Vertex at which to carry out the merge.
      new_edge (Edge, optional): if another Tile has already carried out this 
        merge this should be the resulting new Edge for insertion into this 
        Tile. Defaults to None (when the new Edge will be constructed).
  
    Returns:
      Union[None, tuple]: either None (if a new edge was supplied) or a tuple 
        of the two edge IDs to be removed and the new edge added for return to
        the calling context (i.e. the containing Topology).
    """
    # get the two edge list index positions in which Vertex v is found
    i, j = self.get_edge_IDs_including_vertex(v)
    if new_edge is None: # then we must make a new one
      # also record existing edge IDs to be removed
      to_remove = [self.edges[i].ID, self.edges[j].ID]
      new_edge = self.get_merged_edge(i, j)
      return_edge_updates = True
    else:
      return_edge_updates = False
    if abs(i - j) != 1:
      # edge indices 'wrap' around from end of edge list to start so drop 
      # first and last current edges and stick new one on at the end
      self.edges = self.edges[1:-1] + [new_edge] 
    else:
      # insert new edge into list in place of the two old ones
      self.edges = self.edges[:i] + [new_edge] + self.edges[j+1:]
    # update the edge directions
    self.set_edge_directions()
    if return_edge_updates:
      return to_remove, new_edge
    else:
      return None

  def get_edge_IDs_including_vertex(self, v:"Vertex") -> tuple[int]:
    """Gets the (two) index positions of the edges that include supplied Vertex.

    Args:
        v (Vertex): Vertex of interest.

    Returns:
      tuple[int]: two index positions of Edges in edges list that contain v.
    """
    return (i for i, e in enumerate(self.edges) if v.ID in e.ID)

  def get_merged_edge(self, i:int, j:int) -> "Edge":
    """Returns the new edge resulting from merging the two existing edges at
    index positions i and j in the edges list. For example, if the current list 
    of edge IDs was

        (0 1 2) (4 2) (4 5) (5 0)

    and the merge Vertex was #2, the resulting new edge is constructed from 
    vertices (0 1 2 4). 

    Returns:
      Edge: the requested new Edge.
    """
    # if i and j are not consecutive, then j is predecessor edge
    if abs(i - j) != 1:
      i, j = j, i
    # get edges and their directions
    ei, ej = self.edges[i], self.edges[j]
    CWi, CWj = self.edges_CW[i], self.edges_CW[j]
    # DON'T MESS WITH THIS!!!
    # for predecessors (the head) we want everything including the Vertex 
    # where the merge is occurring; for successors (the tail) we want all but 
    # the first Vertex (which is the one where the merge is occurring).
    # In both cases contingent on whether existing Edges are CW or CCW we
    # may need to flip the Vertex sequence to ensure that the merge Vertex
    # is in the middle of the new edge that will be created
    head = ei.corners if CWi else ei.corners[::-1]
    tail = ej.corners[1:] if CWj else ej.corners[::-1][1:]
    v_sequence = (head if CWi else head[::-1]) + (tail if CWj else tail[::-1])
    return Edge(v_sequence)

  def offset_corners(self, offset:int):
    """Shifts shape, corners, edges, and edges_CW by an offset amount. This is
    used to align tiles that are similar, which is required for correct 
    transfer of 'base' tile labelling on to 'radius 1' tiles during Topology 
    construction.

    Args:
      offset int: the number of positions to shift the lists.
    """
    self.shape = tiling_utils.offset_polygon_corners(self.shape, offset)
    self.corners = self.corners[offset:] + self.corners[:offset]
    self.edges = self.edges[offset:] + self.edges[:offset]
    self.edges_CW = self.edges_CW[offset:] + self.edges_CW[:offset]

  def get_edge_label(self, e:"Edge") -> str:
    """Returns edge label of specified edge.

    Args:
        e (Edge): Edge whose label is required.

    Returns:
        str: requested edge label.
    """
    return self.edge_labels[self.get_edge_IDs().index(e.ID)]

  def get_corner_label(self, v:"Vertex") -> str:
    """Returns corner label of specified corner.

    Args:
        v (Vertex): corner whose label is required.

    Returns:
        str: requested corner label.
    """
    return self.edge_labels[self.get_corner_IDs().index(v.ID)]

  def get_vertex_label_positions(self) -> list[geom.Point]:
    """Returns a viable location at which to position corner labels inside
    tile shape. The method is convoluted because a negative buffer may remove
    colinear corners resulting in fewer positions than we have corners in the
    tile shape!

    Returns:
        list[geom.Point]: list of locations.
    """
    d = (self.shape.area ** 0.5) / 8
    c = self.centre
    corners = [c.point for c in self.corners]
    return [geom.LineString([p, c]).line_interpolate_point(d) for p in corners]
    
  def get_edge_label_positions(self) -> list[geom.Point]:
    """Returns a reasonable location at which to position edge labels inside
    tile shape.

    Returns:
        list[geom.Point]: list of locations
    """
    d = (self.shape.area ** 0.5) / 8
    c = self.centre
    # note direction is important as edge might not be a simple line segment
    sides = [e.get_geometry(CW) for e, CW in zip(self.edges, self.edges_CW)]
    return [geom.LineString(
      [s.centroid, c]).line_interpolate_point(d) for s in sides]

  def angle_at(self, v:"Vertex") -> float:
    """Returns interior angle at the specified corner (in degrees).

    Args:
        v (Vertex): corner where angle is requested.

    Returns:
        float: angle at corner v in degrees.
    """
    i = self.corners.index(v)
    n = len(self.corners)
    return tiling_utils.get_inner_angle(self.corners[i-1].point,
                                        self.corners[i].point,
                                        self.corners[(i + 1) % n].point)


class Vertex:
  """Class to store attributes of a vertex in a tiling."""
  point: geom.Point
  """point (geom.Point): point location of the vertex."""
  ID: int
  """integer (mostly but not necessarily in sequence) of vertex keyed into the 
  points dictionary of the containing Topology."""
  tiles: list["Tile"]
  """list of Tiles incident on this vertex."""
  neighbours: list[int]
  """list of the immediately adjacent other corner IDs. Only required to 
  determine if a point is a tiling vertex (when it will have) three or more 
  neighbours, so only IDs are stored."""
  label: str = None
  """the (upper case letter) label of the vertex under the symmetries of the 
  tiling."""
  is_tiling_vertex: bool = True
  """is_tiling_vertex (bool): True if this is a tiling vertex, rather than a 
  tile corner. E.g., A below is a corner, not a tiling vertex. B is a tiling 
  vertex:
      
      +-------+
      | 1     |
      |   A---B---+
      |   | 2     |
      +---C   +---+
          |   |
          +---+"""
  
  def __init__(self, point:geom.Point, ID:int):
    """Class constructor.

    Args:
      point (geom.Point): point location of the vertex.
      ID (int): a unique integer ID (which will be its key in the containing
        Topology points dictionary).
    """
    self.point = point
    self.ID = ID
    self.tiles = []
    self.neighbours = []

  def __str__(self) -> str:
    """Returns string representation of Vertex.
    
    Returns:
        str: string including ID, point and list of incident Tiles.
    """
    return f"Vertex {self.ID} at {self.point} Tiles: {self.get_tile_IDs()}"

  def __repr__(self) -> str:
    return str(self)

  def get_tile_IDs(self) -> list[int]:
    """Convenience method to return list of Tile IDs (not the Tiles themselves).

    Returns:
        list[int]: list of Tile IDs
    """
    return [t.ID for t in self.tiles]

  def add_tile(self, tile:Tile):
    """Adds supplied Tile to the tiles list if it is not already present.

    Args:
        tile (Tile): Tile to add.
    """
    if not tile in self.tiles:
      self.tiles.append(tile)

  def add_neighbour(self, vertex_id:int):
    """Adds supplied ID to the neighbours list if it is not already present

    Args:
      vertex_id (int): ID to add to the neighbours list.
    """
    if not vertex_id in self.neighbours:
      self.neighbours.append(vertex_id)

  def clockwise_order_incident_tiles(self):
    """Reorders the tiles list clockwise (this is for dual tiling construction)
    """
    cw_order = tiling_utils.order_of_pts_cw_around_centre(
      [t.centre for t in self.tiles], self.point)
    self.tiles = [self.tiles[i] for i in cw_order]

  def is_interior(self) -> bool:
    """Tests if vertex is completely enclosed by its incident Tiles. Based on
    summing the interior angles of the incident tiles at this vertex.

    Returns:
        bool: True if vertex is completely enclosed by incident Tiles.
    """
    return abs(360 - sum([t.angle_at(self) for t in self.tiles])) \
                     < tiling_utils.RESOLUTION


class Edge:
  """Class to represent edges in a tiling (not tile sides) per the definitions
  in Grünbaum and Shephard.
  """
  ID: tuple[int]
  """IDs of the vertices at ends of the edge. Used as key in the containing 
  Topology's edges dictionary."""
  vertices: list["Vertex"]
  """two item list of the end vertices."""
  corners: list["Vertex"]
  """list of all the vertices in the edge (including its end vertices). In a 
  'normal' edge to edge tiling corners and vertices will be identical."""
  right_tile: "Tile" = None
  """the tile to the right of the edge traversed from its first to its last 
  vertex. Given clockwise winding default, all edges will have a right_tile."""
  left_tile: "Tile" = None
  """the tile to the left of the edge traversed from its first to its last 
  vertex. Exterior edges of the tiles in a Topology will not have a left_tile.
  """
  label: str = None
  """the (lower case letter) label of the edge under the symmetries of the 
  tiling."""
  
  def __init__(self, corners:list[Vertex]):
    """Class constructor. Initialises the corners and vertices lists and sets ID
    to (vertices[0].ID, vertices[1].ID). The vertices list is all the corners
    with is_tiling_vertex property True -- Note that during initialisation the
    default of this property is True until after the relations between tiles and
    vertices have been determined.

    Args:
      corners (list[Vertex]): list of all corners along the edge.
    """
    self.corners = corners
    self.vertices = [v for v in self.corners if v.is_tiling_vertex]
    self.ID = tuple(v.ID for v in self.vertices)

  def __str__(self) -> str:
    """Returns a string representation of  the Edge.

    Returns:
        str: include ID and a list of corner vertex IDs.
    """
    return f"Edge {self.ID} Corners: {[c.ID for c in self.corners]}"

  def __repr__(self) -> str:
    return str(self)

  def get_corner_IDs(self) -> list[int]:
    """Convenience method to get the IDs of edge corners.

    Returns:
        list[int]: IDs of all corners.
    """
    return [c.ID for c in self.corners]

  def get_vertex_IDs(self) -> list[int]:
    """Convenience method to get IDs of edge vertices.

    Returns:
        list[int]: list of IDs of the vertices.
    """
    return [v.ID for v in self.vertices]

  def insert_vertex(self, v:"Vertex", predecessor:"Vertex") -> list["Edge"]:
    """Inserts a vertex along this edge after the specified predecessor 
    Vertex and returns this edge modified and a new edge. 
    
    If the initial edge was (say) (0 1 2 5) and the predecessor was set to 1 
    the returned edges would be (0 1 v) and (v 2 5).
    """
    i = self.corners.index(predecessor)
    new_edge = Edge([v] + self.corners[(i+1):])
    if not self.right_tile is None:
      new_edge.right_tile = self.right_tile
    if not self.left_tile is None:
      new_edge.left_tile = self.left_tile
    self.corners = self.corners[:(i+1)] + [v]
    self.vertices = [self.vertices[0], v]
    self.ID = tuple(v.ID for v in self.vertices)
    return [self, new_edge]

  def get_geometry(self, forward = True) -> geom.LineString:
    """Returns a geom.LineString representing the geometry (including all
    corners) of this Edge, optionally starting at either end.

    Args:
      forward (bool, optional): if True the returned LineString starts at
        corners[0], else at corners[-1]. Defaults to True.

    Returns:
      geom.LineString: the required LineString.
    """
    if forward:
      return geom.LineString([v.point for v in self.corners])
    else:
      return geom.LineString([v.point for v in self.corners[::-1]])

  def get_topology(self, forward = True) -> geom.LineString:
    """Returns a LineString connecting the first and last corners (i.e. the
    vertices) of this tilin edge, optionally starting from either end.

    Args:
      forward (bool, optional): if True LineString starts at vertices[0], 
        else at vertices[1]. Defaults to True.

    Returns:
        geom.LineString: the required LineString.
    """
    if forward:
      return geom.LineString([v.point for v in self.vertices])
    else:
      return geom.LineString([v.point for v in self.vertices[::-1]])


class Topology:
  """Class to represent topology of a Tileable object.
  
  NOTE: It is important that get_local_patch return the tileable elements and 
  the translated copies in consistent sequence, i.e. if there are (say) four 
  tiles in the unit, the local patch should be 1 2 3 4 1 2 3 4 1 2 3 4 ... and
  so on. This is because self.tiles[i % n_tiles] is frequently used to reference
  the base unit Tile which corresponds to self.tiles[i].
  """
  tileable: Tileable
  """the Tileable on which the topology will be based."""
  tiles: list[Tile]
  """list of the Tiles in the topology. We use polygons returned by the tileable.get_local_patch method for these. That is the base tiles and 8 adjacent copies (for a rectangular tiling), or 6 adjacent copies (for a hexagonal tiling)."""
  points: dict[int, Vertex]
  """dictionary of all points (vertices and corners) in the tiling, keyed by Vertex ID."""
  edges: dict[int, Edge]
  """dictionary of the tiling edges, keyed by Edge ID."""
  unique_tile_shapes: list[geom.Polygon]
  """a 'reference' tile shape one per tile equivalence class."""
  dual_tiles: list[geom.Polygon]
  """list of geom.Polygons from which a dual tiling might be constructed."""
  n_tiles: int = 0
  """number of tiles in the base Tileable (retained for convenience)."""
  n_patch: int = 0
  """number of tiles in the first order patch of the base Tileable (retained 
  for convenience)."""
  n_tile_groups: int = 0
  """number of tile equivalence classes (as determined by this code base... NOT mathematical theory!)."""
  tile_groups: list[list["Tile"]]
  """list of lists of the tiles distinguished by shape and optionally tile_id"""
  tile_matching_transforms: list[tuple[float]]
  """shapely transform tuples that map tiles onto other tiles"""
  tile_equivalence_classes: list[tuple[int]]
  """list of lists of tile IDs in each equivalence class"""

  def __init__(self, unit: Tileable, ignore_tile_ids:bool = True):
    """Class constructor.

    Args:
      unit (Tileable): the Tileable whose topology is required.
    """
    # Note that the order of these setup steps is fragile and sometimes
    # not obviously so. For example vertex labelling will not operate 
    # correctly if the dual tile generation has not been carried out as
    # this operation reorders tiles around each vertex clockwise.
    self.tileable = unit # keep this for reference
    self.n_tiles = self.tileable.tiles.shape[0]
    self._initialise_points_into_tiles()
    self._setup_vertex_tile_relations()
    self._setup_edges()
    self._copy_base_tiles_to_patch()
    self._identify_distinct_tiles(ignore_tile_ids)
    self._find_tile_equivalence_classes()
    self._label_tiles()
    self._label_vertices()
    self._label_edges()
    self.generate_dual()

  def __str__(self) -> str:
    """Returns string representation of this Topology.
    
    Returns:
      str: a message that recommends examining the tiles, points and edges 
        attributes.

    """
    return f"Topology of Tileable with {self.n_tiles} tiles.\nExamine .tiles, .points and .edges for more details."

  def __repr__(self) -> str:
    return str(self)

  def _initialise_points_into_tiles(self, debug:bool = False) -> None:
    """Sets up dictionary of unique point locations in the tiling and assigns
    them to Tiles.
    
    Args:
      debug (bool): if True prints useful debugging information.
    """
    shapes = self.tileable.get_local_patch(r = 1, include_0 = True).geometry
    shapes = [tiling_utils.get_clean_polygon(s) for s in shapes]
    self.n_patch = len(shapes)
    labels = list(self.tileable.tiles.tile_id) * (len(shapes) // self.n_tiles)
    self.tiles = []
    self.points = {}
    for (i, shape), label in zip(enumerate(shapes), labels):
      tile = Tile(i)
      tile.label = label
      self.tiles.append(tile)
      tile.corners = []
      corners = tiling_utils.get_corners(shape, repeat_first = False)
      for c in corners:
        prev_vertex = None
        for p in reversed(self.points.values()):
          if c.distance(p.point) <= 2 * tiling_utils.RESOLUTION:
            # found an already existing vertex, so add to the tile and break
            prev_vertex = p
            tile.corners.append(prev_vertex)
            break
        if prev_vertex is None:
          # new vertex, add it to the topology dictionary and to the tile
          v = self.add_vertex(c)
          tile.corners.append(v)
          if debug: print(f"Added new Vertex {v} to Tile {i}")

  def _setup_vertex_tile_relations(self, debug:bool = False):
    """Determines relations between vertices and tiles. In particular vertices
    along tile edges that are not yet included in their list of vertices are
    added. Meanwhile vertex lists of incident tiles and neighbours are set up.
    
    Args:
      debug (bool): if True prints debugging information.
    """
    # we do this for all tiles in the radius-1 local patch
    for tile in self.tiles:
      if debug: print(f"Checking for vertices incident on Tile {tile.ID}")
      corners = []
      # performance is much improved using the vertex IDs
      initial_corner_IDs = tile.get_corner_IDs()
      # we need the current shape (not yet set) to check for incident vertices 
      shape = geom.Polygon([c.point for c in tile.corners])
      # get points incident on tile boundary, not already in the tile corners
      new_points = [v for v in self.points.values()
                    if v.ID not in initial_corner_IDs and
                    v.point.distance(shape) <= 2 * tiling_utils.RESOLUTION]
      # iterate over sides of the tile to see which the vertex is incident on
      for c1, c2 in zip(tile.corners, tile.corners[1:] + tile.corners[:1]):
        all_points = [c1, c2]
        if len(new_points) > 0:
          if debug: print(f"{[v.ID for v in new_points]} incident on tile")
          ls = geom.LineString([c1.point, c2.point])
          to_insert = [v for v in new_points 
                       if v.point.distance(ls) <= 2 * tiling_utils.RESOLUTION]
          if len(to_insert) > 0:
            # sort by distance along the side
            d_along = sorted([(ls.line_locate_point(v.point), v)
                              for v in to_insert], key = lambda x: x[0])
            to_insert = [v for d, v in d_along]
            all_points = all_points[:1] + to_insert + all_points[1:]
        corners.extend(all_points[:-1])
        for x1, x2 in zip(all_points[:-1], all_points[1:]):
          # x2 will add the tile and neigbour when we get to the next side
          x1.add_tile(tile)
          x1.add_neighbour(x2.ID)
      tile.corners = corners
      tile.set_shape_from_corners()

  def _setup_edges(self, debug:bool = False):
    """Sets up the tiling edges. 
    
    First vertices in the base tiles are classified
    as tiling vertices or not -- only these can be classified reliably (e.g 
    vertices on the perimeter are tricky).. Up to here all vertices have been 
    considered tiling vertices.
    
    Second edges are created by traversing tile corner lists. Edges are stored
    once only by checking for edges in the reverse direction already in the 
    edges dictionary. Edge right and left tiles are initialised.
    
    Third tile edge direction lists are initialised.
    """
    # classify vertices in the base tiles
    for tile in self.tiles[:self.n_tiles]:
      for v in tile.corners:
        v.is_tiling_vertex = len(v.neighbours) > 2
    if debug: print(f"Classified base tile vertices")
    self.edges = {}
    for tile in self.tiles:
      if debug: print(f"Adding edges from Tile {tile.ID}")
      tile.edges = []
      vertices = [v for v in tile.corners if v.is_tiling_vertex]
      # note that through here finding ints in lists is much faster than
      # finding Vertex objects, hence we use lists of IDs not Vertex objects
      if len(vertices) > 1:
        for v1, v2 in zip(vertices, vertices[1:] + vertices[:1]):
          corner_IDs = tile.get_corner_IDs()
          idx1 = corner_IDs.index(v1.ID)
          idx2 = corner_IDs.index(v2.ID)
          if idx1 < idx2:
            corners = [c for c in corner_IDs[idx1:(idx2 + 1)]]
          else:
            corners = [c for c in corner_IDs[idx1:] + corner_IDs[:(idx2 + 1)]]
          ID = (corners[0], corners[-1])
          if not ID in self.edges:
            # check that reverse direction edge is not present first
            r_ID = ID[::-1]
            if r_ID in self.edges:
              # if it is, then add it and set left_tile
              if debug: print(f"reverse edge {r_ID} found")
              e = self.edges[r_ID]
              e.left_tile = tile
              tile.edges.append(e)
            else:
              # we've found a new edge so make and add it
              if debug: print(f"adding new edge {corners}")
              e = self.add_edge([self.points[c] for c in corners])
              e.right_tile = tile
              tile.edges.append(e)
      # finally initialise the edge direction information in the tile
      tile.set_edge_directions()

  def _copy_base_tiles_to_patch(self):
    """Copies attributes of base tiles to their corresponding tiles in 
    the radius-1 patch. This requires:
     
    First, inserting any additional corners in the base tiles not found in the 
    radius-1 tiles. 
    
    Second, any vertices in the base tiles that are NOT tiling vertices are 
    applied to radius-1 tiles leading to merging of some edges into.
    """
    # the number of tiles in the base + radius-1
    n_r1 = len(self.tiles)
    # first add any missing vertices to the non-base tiles
    # we add all the missing vertices before doing any merges
    for base in self.tiles[:self.n_tiles]:
      for other in self.tiles[base.ID:n_r1:self.n_tiles]:
        self._match_reference_tile_vertices(base, other)
    # then merge any edges that meet at a corner
    for base in self.tiles[:self.n_tiles]:
      for other in self.tiles[base.ID:n_r1:self.n_tiles]:
        self._match_reference_tile_corners(base, other)

  def _match_reference_tile_vertices(self, tile1:Tile, tile2:Tile):
    """Adds vertices to tile2 so that it matches tile1 adjusting edges as
    required. This assumes the tiles are the same shape, but that tile2 may
    be missing some tiling vertices along some edges.

    Args:
      tile1 (Tile): reference tile.
      tile2 (Tile): tile to change.
    """
    to_add = len(tile1.corners) - len(tile2.corners)
    while to_add > 0:
      # find the reference x-y offset 
      dxy = (tile2.centre.x - tile1.centre.x, tile2.centre.y - tile1.centre.y)
      for i, t1c in enumerate([c.point for c in tile1.corners]):
        t2c = tile2.corners[i % len(tile2.corners)].point
        if abs((t2c.x - t1c.x) - dxy[0]) > 10 * tiling_utils.RESOLUTION or \
           abs((t2c.y - t1c.y) - dxy[1]) > 10 * tiling_utils.RESOLUTION:
          # add vertex to t2 by copying the t1 vertex appropriately offset
          # note that this might alter the length ot t2.corners
          v = self.add_vertex(geom.Point(t1c.x + dxy[0], t1c.y + dxy[1]))
          v.is_tiling_vertex = True
          old_edge, new_edges = tile2.insert_vertex_at(v, i)
          del self.edges[old_edge]
          for e in new_edges:
            self.edges[e.ID] = e
          to_add = to_add - 1

  def _match_reference_tile_corners(self, tile1:Tile, tile2:Tile):
    """Finds vertices that are corners in tile2 but vertices in tile1 and 
    updates tile2 to match - merging edges as required.

    Args:
        tile1 (Tile): reference tile.
        tile2 (Tile): tile to make match.
    """
    vs_to_change = [vj for vi, vj in zip(tile1.corners, tile2.corners)
                    if not vi.is_tiling_vertex and vj.is_tiling_vertex]
    if len(vs_to_change) > 0:
      for v in vs_to_change:
        v.is_tiling_vertex = False
        # it's a corner not an edge so will have no more than 2 v.tiles
        old_edges, new_edge = v.tiles[0].merge_edges_at_vertex(v)
        for e in old_edges:
          del self.edges[e]
        self.edges[new_edge.ID] = new_edge

  def _identify_distinct_tiles(self, ignore_tile_id_labels:bool = True):
    """Determines unique tiles based on their symmetries and shapes. At the
    same time assembles a list of the affine transforms under which matches
    occurs since these are potential symmetries of the tiling.
    
    Args:
      ignore_tile_id_labels (bool): if True only the shape of tiles matters;
        if False the tile_id label is also considered. Defaults to True.
    """
    self.unique_tile_shapes = []
    offsets = []
    self.tile_matching_transforms = {}
    n_symmetries = 0
    for tile in self.tiles[:self.n_tiles]:
      s = Symmetries(tile.shape)
      transforms = [s.get_matching_transforms(other.shape) 
                    for other in self.unique_tile_shapes]
      shape_matches = [not t is None for t in transforms]
      if ignore_tile_id_labels:
        label_matches = [True] * len(self.unique_tile_shapes)
      else:
        label_matches = [tile.label == other.label
                        for other in self.unique_tile_shapes]
      if any([x and y for x, y in zip(shape_matches, label_matches)]):
        match = shape_matches.index(True)
        offset = transforms[match]["offset"]
        offsets.append(offset)
        tile.group = match
        tile.offset_corners(offset)
        for s in transforms[match]["symmetries"]:
          self.tile_matching_transforms[n_symmetries] = s
          n_symmetries = n_symmetries + 1
      else:
        self.unique_tile_shapes.append(tile)
        offsets.append(0)
        tile.group = self.n_tile_groups
        self.n_tile_groups = self.n_tile_groups + 1
    # now copy assignments from the central TileUnit to everything else
    # and rotate tiles if required to ensure their corners match
    for tile in self.tiles[self.n_tiles:]:
      tile.group = self.tiles[tile.ID % self.n_tiles].group
      tile.offset_corners(offsets[tile.ID % self.n_tiles])
    # also add into the list of tile matching transforms the symmetries
    # of the prototile, which might also be symmetries of the tiling
    for s in Symmetries(self.tileable.prototile.geometry[0]).symmetries:
      self.tile_matching_transforms[n_symmetries] = s.transform
      n_symmetries = n_symmetries + 1
    # filter back the tile matching transforms to a unique set
    self._remove_duplicate_symmetries()
    self.tile_groups = [[t.ID for t in self.tiles if t.group == group]
                        for group in range(self.n_tile_groups)]

  def _remove_duplicate_symmetries(self):
    """Filters the supplied list of shapely affine transforms so that no 
    duplicates are retained.
    """
    uniques = []
    for k, v in self.tile_matching_transforms.items():
      already_exists = False
      for u in uniques:
        already_exists = np.allclose(v, u[1], atol = 1e-6, rtol = 1e-6)
        if already_exists:
          break
      if not already_exists:
        uniques.append((k, v))
    self.tile_matching_transforms = dict(uniques)

  def _find_tile_equivalence_classes(self):
    base_tiles = [t for t in self.tiles[:self.n_tiles]]
    patch_tiles = [t for t in self.tiles[:self.n_patch]]
    equivalents = defaultdict(set)
    # every tile is in its own set of equivalents, so add them
    for group in self.tile_groups:
      group_symmetries = []
      in_group_equivalents = {}
      # make up source and target tiles within this group
      source_tiles = [t for t in base_tiles if t.ID in group]
      target_tiles = [t for t in patch_tiles if t.ID in group]
      for tr_i, transform in self.tile_matching_transforms.items():
        matched_geoms = {}
        possible_symmetry = True
        for source in source_tiles:
          match_id = self._match_under_transform(
            source.shape, target_tiles, transform)
          if match_id == -1:
            possible_symmetry = False
          else:
            matched_geoms[source.ID] = match_id
        if possible_symmetry:
          for k, v in matched_geoms.items():
            in_group_equivalents[(tr_i, k)] = match_id % self.n_tiles
          group_symmetries.append(tr_i)
      # any transform not a symmetry of group is not a symmetry of tiling
      # so remove from consideration and also from the list of tile matches
      self.tile_matching_transforms = {
        i: self.tile_matching_transforms[i] for i in group_symmetries}
      for k, v in in_group_equivalents.items():
        equivalents[k[1]].add(v)
    self.tile_equivalence_classes = \
      list(set([tuple(sorted(v)) for v in equivalents.values()]))
    for c, eclass in enumerate(self.tile_equivalence_classes):
      for tile in self.tiles:
        if tile.ID % self.n_tiles in eclass:
          tile.equivalence_class = c
      

  def _match_under_transform(self, tile:Tile, patch:list[Tile], 
                             transform:tuple[float]) -> int:
    """Determines if there is a tile in the supplied patch onto which the
    supplied tile is mapped by the supplied symmetry.

    Args:
        tile (Tile): _description_
        patch (list[Tile]): _description_
        transform (tuple[float]): _description_

    Returns:
      int: ID of the tile in patch that matches tile under the transform if one
        exists, otherwise returns -1.
    """
    match_id = -1
    for tile2 in patch:
      match = tiling_utils.geometry_matches(
        affine.affine_transform(tile, transform), tile2.shape)
      if match:
        return tile2.ID
    return match_id

  def _label_tiles(self):
    """Labels the base tile vertices, then copies to corresponding tile in the 
    radius-1 patch.
    """
    first_letter = 0
    for group in range(self.n_tile_groups):
      tiles = [t for t in self.tiles[:self.n_tiles] if t.group == group]
      s = Symmetries(self.unique_tile_shapes[group].shape)
      vlabels = list(s.get_unique_labels(offset = first_letter)["rotations"][0])
      elabels = self._get_edge_labels_from_vertex_labels(vlabels)
      for tile in tiles:
        tile.vertex_labels = vlabels
        tile.edge_labels = elabels
      first_letter = first_letter + len(set(vlabels))
    for tile in self.tiles: #[self.n_tiles:]:
      tile.vertex_labels = self.tiles[tile.ID % self.n_tiles].vertex_labels
      tile.edge_labels = self.tiles[tile.ID % self.n_tiles].edge_labels

  def _get_edge_labels_from_vertex_labels(self, vlabels:list[str]) -> list[str]:
    r"""Given a list of vertex labels, returns a list of edge labels such that
    each distinct pair of vertex labels at the ends of an edge will be given a
    distinct lower case letter label. E.g., A B C B A will yield a b b a c from 
    AB -> a, BC -> b, CB -> b, BA -> a, AA -> c

    Args:
      vlabels (list[str]): ordered list of vertex labels.

    Returns:
      list[str]: corresponding list of edge labels.
    """
    edge_labels = [a + b for a, b in zip(vlabels, vlabels[1:] + vlabels[:1])]
    letter = ALPHABET.index(min(list(vlabels)))
    elabels = {}
    for e in edge_labels:
      if not e in elabels:
        # if e[::-1] in elabels:
        #   label = elabels[e[::-1]]
        #   elabels[e] = label + "-"
        #   elabels[e[::-1]] = label + "+"
        if e[::-1] in elabels:
          elabels[e] = elabels[e[::-1]]
        else:
          elabels[e] = alphabet[letter]
          letter = letter + 1
    return [elabels[l] for l in edge_labels]

  def _label_vertices(self):
    """Labels vertices based on cycle of tile vertex labels around them.
    """
    uniques = set()
    # first label vertices in the core tiles
    vs = set()
    for t in self.tiles[:self.n_tiles]:
      vs = vs.union(t.get_corner_IDs())
    # Note: sorting is important for repeatable labelling!
    for vi in sorted(vs):
      v = self.points[vi]
      label = ""
      for tile in v.tiles:
        label = label + \
          tile.vertex_labels[tile.corners.index(v)]
      # TODO: resolve this question: cyclic sort seems more correct,
      # but neither approach seems to work in all cases... see esp.
      # cyclic sort applied to the cheese sandwich tiling.
      v.label = "".join(self._cyclic_sort_first(list(label)))
      # v.label = min(list(label))
      uniques.add(v.label)
    for vi in sorted(vs):
      v = self.points[vi]
      v.label = ALPHABET[sorted(uniques).index(v.label)]
    # now copy to corresponding vertices in the rest of tiling
    for ti, t in enumerate(self.tiles):
      if ti >= self.n_tiles:
        t0 = self.tiles[ti % self.n_tiles]
        for v0, v in zip(t0.corners, t.corners):
          # Vertices may appear in more than one tile, only label once!
          if self.points[v.ID].label is None:
            self.points[v.ID].label = self.points[v0.ID].label

  def _label_edges(self):
    """Labels edges based on the tile edge label on each side.
    """
    uniques = set()
    labelled = set()
    # first label base tile edges from the central tile unit
    for t in self.tiles[:self.n_tiles]:
      for e in t.edges:
        rt, lt = e.right_tile, e.left_tile
        rlab, llab = rt.get_edge_label(e), lt.get_edge_label(e)
        e.label = "".join(sorted([rlab, llab]))
        uniques.add(e.label)
        labelled.add(e.ID)
    for ei in sorted(labelled):
      e = self.edges[ei]
      e.label = alphabet[sorted(uniques).index(e.label)]
    # now copy to corresponding edges in the rest of tiling
    for ti, t in enumerate(self.tiles):
      if ti >= self.n_tiles:
        t0 = self.tiles[ti % self.n_tiles]
        for e0, e in zip(t0.edges, t.edges):
          # edges may appear in more than one tile, only label once!
          if e.label is None:
            e.label = e0.label

  def _cyclic_sort_first(self, lst:Iterable) -> Iterable:
    """Returns supplied Iterable in canonical cyclic sorted form. E.g. the sequence ACABD, yields 5 possible cycles ACABD, CABDA, ABDAC, BDACA, and 
    DACAB. The lexically first of these is ABDAC, which would be returned. 

    Args:
        lst (Iterable): an Iterable of elements (usu. strings).

    Returns:
      Iterable: the lexically first of the possible cycles in the supplied 
        iterable.
    """
    cycle = lst + lst
    n = len(lst)
    cycles = [cycle[i:i + n] for i in range(n)]
    return sorted(cycles)[0]

  def add_vertex(self, pt:geom.Point) -> "Vertex":
    """Adds a Vertex at the specified location, returning it to the caller. No
    attempt is made to ensure Vertex IDs are an unbroken sequemce: a new ID is 
    generated one greater than the existing highest ID.

    Args:
        pt (geom.Point): point location of the Vertex.

    Returns:
        Vertex: the added Vertex object.
    """
    n = 0 if len(self.points) == 0 else max(self.points.keys()) + 1
    v = Vertex(pt, n)
    self.points[n] = v
    return v

  def add_edge(self, vs:list["Vertex"]) -> "Edge":
    """Adds an Edge to the edges dictionary. Edges are self indexing by the IDs
     of their end Vertices. Returns the new Edge to the caller.

    Args:
      vs (list[Vertex]): list of Vertices in the Edge to be created.

    Returns:
        Edge: the added Edge.
    """
    e = Edge(vs)
    self.edges[e.ID] = e
    return e

  def _get_tile_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries([t.shape for t in self.tiles])

  def _get_tile_centre_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries([t.centre for t in self.tiles])

  def _get_point_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries([v.point for v in self.points.values()])

  def _get_vertex_geoms(self) -> gpd.GeoSeries:
    return gpd.GeoSeries([v.point for v in self.points.values()
                          if v.is_tiling_vertex])

  def _get_edge_geoms(self, offset:float = 0.0) -> gpd.GeoSeries:
    return gpd.GeoSeries([e.get_topology().parallel_offset(offset)
                          for e in self.edges.values()])

  def generate_dual(self) -> list[geom.Polygon]:
    """Create the dual tiiing for the tiling of this Topology.
    
    TODO: make this a viable replacement for the existing dual tiling 
    generation.

    Returns:
      list[geom.Polygon]: a list of polygon objects.
    """
    for v in self.points.values():
      v.clockwise_order_incident_tiles()
    self.dual_tiles = {}
    for v in self.points.values():
      if v.is_interior() and len(v.tiles) > 2:
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
      self._get_tile_geoms().plot(
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
        edges = self._get_edge_geoms(dist / 100 if offset_edges else 0)
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
                      selector:str, type:str, **kwargs) -> "Topology":
    r"""Applies a specified transformation of Edges in the Topology whose labels
    match the label parameter, optionally applying the transform to update tile
    shapes, and optionally returning a new Topology object (or applying it to 
    this one). 
    
    Implemented in this way so that transformations can be applied one at a time
    without creating an intermediate set of new tiles, which may be invalid and
    fail. So, if you wish to apply (say) 3 transforms and generate a new 
    Topology leaving the existing one intact:
    
        new_topo = old_topo.transform_edges(True,  False, "a", ...) \
                           .transform_edges(False, False, "b", ...) \
                           .transform_edges(False, True,  "c", ...)
    
    The first transform requests a new Topology, subsequent steps do not, and it
    is only the last step which attempts to create the new tile polygons.
    
    **kwargs supply named parameters for the requested type of transformation.

    Args:
      new_topology (bool): if True returns a new Topology object, else returns 
        the current Topology modified.
      apply_to_tiles (bool): if True attempts to create new Tiles after the 
        transformation has been applied. Usually set to False, unless the last
        transformation in a pipeline, to avoid problems of topologically invalid
        tiles at intermediate steps.
      selector (str): label of edges to which to apply the transformation. Note
        that all letters in the supplied string are checked, so you can use e.g.
        "abc" to apply a transformation to edges labelled "a", "b" or "c".
      type (str): name of the type of transformation requested. Currently only
        "zigzag" is supported. Different types will expect different kwargs
        per the associated <type>_edge method.

    Returns:
      Topology: if new_topology is True a new Topology based on this one with
        after transformation, if False this Topology is returned after the
        transformation.
    """
    # TODO: check if this is the best way to make a deep copy...
    topo = pickle.loads(pickle.dumps(self)) if new_topology else self
    for e in topo.edges.values():
      if e.label in selector:
        if type == "zigzag":
          topo.zigzag_edge(e, **kwargs)
    if apply_to_tiles:
      for t in topo.tiles:
        t.set_corners_from_edges()
    topo.tileable.tiles.geometry = gpd.GeoSeries(
      [tiling_utils.get_clean_polygon(topo.tiles[i].shape) 
       for i in range(topo.n_tiles)])
    topo.tileable.setup_regularised_prototile_from_tiles()
    return topo

  def zigzag_edge(self, edge:"Edge", start:str = "A",
                  n:int = 2, h:float = 0.5, smoothness:int = 0):
    """Applies a zigzag transformation to the supplied Edge. Currently this will
    only work correctly if h is even.
    
    TODO: make it possible for odd numbers of 'peaks' to work (this may require
    allowing bidirectional Edges, i.e. storing Edges in both directions so that
    all Tile edges are drawn CW). The `start` parameter is a temporary hack for 
    this

    Args:
      edge (Edge): Edge to transform
      n (int, optional): number of zigs and zags in the edge. Defaults to 2.
      start (str, optional): label at one end of edge which is used to determine
        the sense of h, enabling C-curves with an odd number n of zigs and zags to be applied. Defaults to 'A'.
      h (float, optional): width of the zig zags relative to edge length. 
        Defaults to 0.5.
      smoothness (int, optional): spline smoothness. 0 gives a zig zag proper,
        higher values will produce a sinusoid. Defaults to 0.
    """
    v0, v1 = edge.vertices[0], edge.vertices[1]
    if n % 2 == 1 and v0.label != start:
      h = -h
    ls = tiling_utils.zigzag_between_points(v0.point, v1.point, 
                                            n, h, smoothness)
    # remove current corners
    self.points = {k: v for k, v in self.points.items()
                   if not k in edge.get_corner_IDs()[1:-1]}
    # add the new ones
    new_corners = [self.add_vertex(geom.Point(xy)) for xy in ls.coords]
    edge.corners = edge.vertices[:1] + new_corners + edge.vertices[-1:]
    if not edge.right_tile is None:
      edge.right_tile.set_corners_from_edges(False)
    if not edge.left_tile is None:
      edge.left_tile.set_corners_from_edges(False)
