#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass

import shapely.geometry as geom
import topojson

from weavingspace.tileable import Tileable

@dataclass
class Topology:
  """Class to represent topology of a Tileable.
  """

  unit: Tileable = None
  topology: dict = None
  vertices: dict[tuple[int], tuple[float]] = None
  edges: list[dict] = None
  faces: list[int] = None


  def __init__(self, unit:Tileable):
    self.unit = unit
    self.set_topology()
    self.set_vertices_and_edges()
    self.set_faces()
    self.set_edge_labels()

  def set_topology(self):
    patch = self.unit.get_local_patch(r = 1, include_0 = True)
    self.topology = topojson.Topology(patch, prequantize = True).to_dict()
    return

  def _get_index(self, x, iterable):
    if x in iterable:
      return iterable.index(x)
    else:
      return -1

  def set_vertices_and_edges(self):
    sx, sy = self.topology["transform"]["scale"]
    tx, ty = self.topology["transform"]["translate"]
    self.vertices = {}
    edges = []
    for arc in self.topology["arcs"]:
      pts = [tuple(arc[0])]
      for dxdy in arc[1:]:
        pts.append((pts[-1][0] + dxdy[0], pts[-1][1] + dxdy[1]))
      for p in pts:
        self.vertices[p] = geom.Point(tx + sx * p[0], ty + sy * p[1])
      edges.append(geom.LineString([self.vertices[p] for p in pts]))
    self.edges = {}
    for i, e in enumerate(edges):
      self.edges[i] = {"geometry": e}
      self.edges[-i - 1] = {"geometry": e.reverse()}

  def set_faces(self):
    self.faces = []
    poly_id = 0
    for obj in self.topology["objects"]["data"]["geometries"]:
      pts = []
      for arc in obj["arcs"][0]:
        new_pts = (self.edges[arc]["geometry"].coords)
        if arc < 0:
          self.edges[arc]["right"] = poly_id
          self.edges[-arc -1]["left"] = poly_id
        else:
          self.edges[arc]["right"] = poly_id
          self.edges[-arc - 1]["left"] = poly_id
        if len(pts) == 0:
          pts.extend(new_pts)
        else:
          pts.extend(new_pts[1:])
      self.faces.append({
        "geometry": geom.Polygon([geom.Point(p) for p in pts]),
        "element_id": obj["properties"]["element_id"]})
      poly_id = poly_id + 1

  def set_edge_labels(self):
    for e in self.edges.values():
      if "left" in e and "right" in e:
        e["label"] = self.faces[e["left"]]["element_id"] + \
                     self.faces[e["right"]]["element_id"]
      else:
        e["label"] = ""
