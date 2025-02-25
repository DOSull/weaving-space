#!/usr/bin/env python
# coding: utf-8

import shapely.geometry as geom
# import weavingspace.tiling_utils as tiling_utils
from weavingspace import tiling_utils

class Tile(object):
  """Class to capture and manipulate essential features of polygons in a tiling.
  """
  ID: int
  """integer ID number which indexes the Tile in the containing Topology tiles 
  list."""
  base_ID: int
  """ID of corresponding Tile in the base tileable unit"""
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
  shape_group: int = None
  """the tile shape group of this tile in its containing Topology."""
  transitivity_class: int = None
  """the tile transitivity class of this tile its containing Topology"""
  
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
    # the first Vertex (which is the one where the merge is occurring). In both 
    # cases contingent on whether existing Edges are CW or CCW we may need to 
    # flip the Vertex sequence to ensure that the merge Vertex is in the middle 
    # of the new edge that will be created
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
    if not offset is None or offset == 0:
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
    return [geom.LineString([s.centroid, c]).line_interpolate_point(d) 
            for s in sides]

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
  base_ID: int = 1_000_000
  """ID of corresponding Vertex in the tileable base_unit"""
  transitivity_class: int = None
  """transitivity class of the vertex under symmetries of the tiling"""
  label: str = ""
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
    self.base_ID = self.ID
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
  base_ID: tuple[int] = (1_000_000, 1_000_000)
  """ID of corresponding edge in the base tileable"""
  transitivity_class: int = None
  """transitivity class of the edge under symmetries of the tiling"""
  label: str = ""
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

