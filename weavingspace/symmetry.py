#!/usr/bin/env python
# coding: utf-8

from typing import Iterable
from typing import Any
from typing import Callable
from dataclasses import dataclass
import string
import copy

import numpy as np
import matplotlib.pyplot as pyplot

import shapely.geometry as geom
import shapely.affinity as affine
import shapely.ops
import geopandas as gpd

import weavingspace.tiling_utils as tiling_utils


class KMP_Matcher:
  sequence: Iterable

  def __init__(self, sequence:Iterable):
    self.sequence = sequence

  def find_matches(self, pat:Iterable[tuple[float]]) -> list[int]:
    """ Implements Knuth-Morris-Pratt string pattern matching algorithm. See:
    https://en.wikipedia.org/wiki/Knuth–Morris–Pratt_algorithm which provides
    detailed pseudo-code on which this code is directly based. See also:
    
    Knuth DE, JH Morris Jr, and VR Pratt. 1977. Fast pattern  matching in 
    strings. SIAM Journal on Computing 6(2): 323–350. doi: 10.1137/0206024.

    This implementation expects sequences of tuples of floats, although in 
    principle any objects could be contained in the Iterables.

    Args:
      pat (Iterable[tuple[float]]): The sequence to match.

    Returns:
      Iterable[int]: _description_
    """
    j, k, = 0, 0
    finds = []
    table = self._get_table(pat)
    while j < len(self.sequence):
      if self._equal_tuples(pat[k], self.sequence[j]):
        j = j + 1
        k = k + 1
        if k == len(pat):
          finds.append(j - k)
          k = table[k]
      else:
        k = table[k]
        if k < 0:
          j = j + 1
          k = k + 1
    return finds

  def _get_table(self, pattern:Iterable) -> Iterable[int]:
    """Returns the 'offsets' table used by KMP pattern matching algorithm as
    required by self.find_matches, based on the supplied pattern.

    Args:
      pattern (Iterable): The pattern to set up the table for.
      equals (function, optional): A function to use to test for equality of
        elements in patterns. Defaults to _equal_tuples().

    Returns:
        Iterable[int]: _description_
    """
    pos = 1
    cnd = 0
    T = {0: -1}
    while pos < len(pattern):
      if np.allclose(pattern[pos], pattern[cnd]):
        T[pos] = T[cnd]
      else:
        T[pos] = cnd
        while cnd >= 0 and not self._equal_tuples(pattern[pos], pattern[cnd]):
          cnd = T[cnd]
      pos = pos + 1
      cnd = cnd + 1
    T[pos] = cnd
    return tuple(T.values())

  def _equal_tuples(self, t1:Iterable[float], t2:Iterable[float]) -> bool:
    """Tests for near equality of two iterables of floats using numpy
    allclose. Wrapped so as to allow find_matches function to take 
    alternative equality tests as input.
    """
    return np.allclose(t1, t2, atol = tiling_utils.RESOLUTION * 10)

  def _round_tuple(self, t:Iterable[float], digits:int = 3) -> Iterable[float]:
    """Convenience function to round all members of an interable. Used for 
    display purposes."""
    return tuple([np.round(x, digits) for x in t])


@dataclass
class Symmetry():
  type: str = "identity"
  description: str = "none"
  angle: float = 0.0
  centre: geom.Point = geom.Point(0, 0)
  transform: tuple[float] = (1, 0, 0, 1, 0, 0)
  shift: int = 0
  
  def apply(self, geometry:Any) -> Any:
    return affine.affine_transform(geometry, self.transform)

  def get_recentred(self, centre:geom.Point) -> "Symmetry":
    new_symmetry = copy.deepcopy(self)
    new_symmetry.centre = centre
    if new_symmetry.type == "rotation":
      new_symmetry.transform = tiling_utils.get_rotation_transform(
        new_symmetry.angle, (new_symmetry.centre.x, new_symmetry.centre.y))
    else:
      new_symmetry.transform = tiling_utils.get_reflection_transform(
        new_symmetry.angle, (new_symmetry.centre.x, new_symmetry.centre.y))
    return new_symmetry

  def plot(self, ax:pyplot.axes, geometry:geom.Polygon) -> pyplot.axes:
    bb = geometry.bounds
    w, h = bb[2] - bb[0], bb[3] - bb[1]
    gpd.GeoSeries([geometry]).plot(ax = ax, fc = "#ffffff00", 
                                   ec = "k", lw = 0.5)
    if self.type == "rotation":
      if self.angle > 0:
        rotn = np.radians(self.angle) 
        arc = geom.LineString(
          [[w/6 * np.cos(a), w/6 * np.sin(a)]
            for a in np.linspace(0, rotn, 50)])
        angle = geom.LineString(
          [(w/3, 0), (0, 0), 
           (w/4 * np.cos(rotn), w/4 * np.sin(rotn))])
        ls = [affine.translate(l, self.centre.x, self.centre.y)
              for l in [arc, angle]]
        gpd.GeoSeries(ls).plot(ax = ax, ec = "b", lw = 0.35)
    else:
        ls = geom.LineString([(min(-h, -w), 0), (max(h, w), 0)])
        ls = affine.rotate(ls, self.angle)
        ls = affine.translate(ls, self.centre.x, self.centre.y)
        ls = ls.intersection(geometry)
        gpd.GeoSeries([ls]).plot(ax = ax, ec = "r", ls = "dashed", lw = 0.5)
    pyplot.axis("off")
    return ax


@dataclass
class Symmetries():

  polygon:geom.Polygon = None
  matcher:KMP_Matcher = None
  n:int = None
  p_code:list[tuple[float]] = None
  p_code_r:list[tuple[float]] = None
  rotation_shifts:list[int] = None
  reflection_shifts:list[int] = None
  symmetries:list[Symmetry] = None


  def __init__(self, polygon:geom.Polygon):
    self.polygon = polygon
    self.n = shapely.count_coordinates(self.polygon) - 1
    self.p_code = self._get_polygon_code(self.polygon)
    self.p_code_r = self._get_polygon_code(self.polygon, mirrored = True)
    self.matcher = KMP_Matcher(self.p_code + self.p_code[:-1])
    self.symmetries = self.get_symmetries()
    

  def _get_polygon_code(self, p:geom.Polygon, 
                        mirrored = False) -> list[tuple[float]]:
    r"""Returns a list of length-angle pairs to uniquely encode a polygon
    shape (up to scale). The alignment of lengths and angles is important
    and symmetry detection is sensitively dependent on it...

    The forward (i.e. unmirrored) polygon pairs are of edge lengths and the 
    angle between each edge and its successor going CW around the polygon. 
    
           A[i] ----L[i]---- A (i+1)
          /                   \
         /                     \
                             
    i.e. edge i is between angles i and i + 1, and angle i is between edge
    i-1 and i. The unmirrored encoding pairs length[i] with angle[i+1].
    
    The mirrored encoding pairs matched indexes of lengths and angles in the
    original polygon. This means the angle is still the one between and edge 
    and its successor but proceeding CCW around the polygon.
    
    See in particular for clarification.
    
    Eades P. 1988. Symmetry Finding Algorithms. In Machine Intelligence and 
    Pattern Recognition, ed. GT Toussaint, 6:41-51. Computational Morphology. 
    North-Holland. doi: 10.1016/B978-0-444-70467-2.50009-6.
    
    Args:
      p (geom.Polygon): the polygon to encode
        mirrored (bool, optional): if true encoding will be in CCC order. 
        Defaults to False.

    Returns:
      list[tuple[float]]: _description_
    """
    # get edge lengths and angles
    lengths = tiling_utils.get_side_lengths(p)
    raw_angles = tiling_utils.get_interior_angles(p)
    if mirrored:
      lengths_r = list(reversed(lengths))
      angles_r = list(reversed(raw_angles))
      return list(zip(lengths_r, angles_r))
    else:
      angles = raw_angles[1:] + raw_angles[:1]
      return list(zip(lengths, angles))


  def get_symmetries(
      self, other_polygon:geom.Polygon = None) -> tuple[list[int], list[int]]:
    """Finds rotation and reflection symmetries of the supplied polygon.
    Based on
    
      Eades P. 1988. Symmetry Finding Algorithms. In Machine Intelligence and 
      Pattern Recognition, ed. GT Toussaint, 6:41-51. Computational Morphology. 
      North-Holland. doi: 10.1016/B978-0-444-70467-2.50009-6.
    
    and also
    
      Wolter JD, TC Woo, and RA Volz. 1985. Optimal algorithms for symmetry 
      detection in two and three dimensions. The Visual Computer 1(1): 37-48. 
      doi: 10.1007/BF01901268.
    
    Details in these papers are not unambiguous. This implementation was
    developed based on them, but with a lot of trial and error to get the index
    offsets and (especially) retrieval of the reflection axes angles to work.

    Args:
      other (geom.Polygon, optional): finds transforms between this polygon and 
        another that may be supplied. Defaults to None, when the comparison 
        polygon will be that for which this Symmetries object has been  
        initialised.

    Returns:
      dict[str, list[float]]: _description_
    """
    # compose extended sequence of length-angle pairs for cyclic matching
    # S = self.p_code + self.p_code[:-1]
    # get code for the 'other' polygon, if it exists...
    if not other_polygon is None:
      p_code = self._get_polygon_code(other_polygon)
      p_code_r = self._get_polygon_code(other_polygon, mirrored = True)
    else:
      p_code = self.p_code
      p_code_r = self.p_code_r
    # ==== Rotations ====
    self.rotation_shifts = self.matcher.find_matches(p_code)
    # ==== Reflections ====
    self.reflection_shifts = self.matcher.find_matches(p_code_r)
    return self.get_rotations(self.rotation_shifts) + \
           self.get_reflections(self.reflection_shifts)


  def get_rotations(self, offsets:list[int]) -> list[Symmetry]:
    """Gets the rotations associated with this collection of symmetries.

    Returns:
      list[Symmetry]: A list of the rotation symmetries associated with this
        polygon.
    """
    c = self.polygon.centroid
    rot_angles = [np.round(i * 360 / self.n, 6) for i in offsets] 
    rotations = [
      np.round(tiling_utils.get_rotation_transform(a, (c.x, c.y)), 6) 
      for a in rot_angles]
    return [Symmetry("rotation", f"{angle}° around ({c.x},{c.y})", 
                     angle, c, rotation, offset)
            for rotation, angle, offset
            in zip(rotations, rot_angles, offsets)]
    

  def get_reflections(self, offsets:list[int]) -> list[Symmetry]:
    """Gets the reflections associated with this collection of symmetries.

    Returns:
      list[Symmetry]: A list of the reflection symmetries associated with this
        polygon.
    """
    # calculate all possible reflection axes - ordering of these is important
    # for correct picking - don't mess with this code!
    c = self.polygon.centroid
    bearings = tiling_utils.get_side_bearings(self.polygon)
    reflection_axes = []
    interior_angles = tiling_utils.get_interior_angles(self.polygon)
    for i in range(self.n):
      # angle bisectors
      reflection_axes.append(bearings[i] - interior_angles[i] / 2)
      # edge_perpendicular bisectors
      reflection_axes.append(bearings[i] - 90)  # normal to the edge
    # pick out those where matches occurred
    ref_angles = [np.round(reflection_axes[i], 6) for i in offsets]
    reflections = [
      np.round(tiling_utils.get_reflection_transform(a, (c.x, c.y)), 6)
      for a in ref_angles]
    return [Symmetry("reflection", f"in line at {angle}° thru ({c.x},{c.y})", 
                     angle, c, reflection, offset)
            for reflection, angle, offset 
            in zip(reflections, ref_angles, offsets)]


  def get_matching_transforms(self, 
                              other:geom.Polygon) -> dict[str, list[float]]:
    """Finds the transforms that will map another polygon onto this one.
    Reuses the _find_symmetries() method, but supplying the other polygon
    as that method's optional argument. Before determining the match the
    points on the perimeter of the other polygon are reordered so that the
    index positions match those of the reference polygon.

    Args:
        other (geom.Polygon): the polygon to match.

    Returns:
      dict[str,Union[str,list[float]]]: A dictionary with the following entries:
      'rotation-shifts': vertex offsets to form matches by rotation
      'reflection-shifts': vertex offsets to form matches by reflection
      'rotation-angles': list of rotation angles of matches (in degrees)
      'reflection-angles': list of angles of reflection axes (in degrees)
      'rotation-transforms': list of shapely affine transform tuples of    
        rotation matches
      'reflection-transforms': list of shapely affine transform tuples of
        reflection matches
      'pre-translation': translation to align other with this polygon as tuple
        of floats
      'pre-rotation': rotation to align other with this polygon (in degrees)
      'pre-transform': shapely affine transform tuple that brings other into
        alignment with this polygon such that the matching transforms work 
    """
    if self.n == shapely.count_coordinates(other) - 1:
      # shift them to the same centroid
      c0 = self.polygon.centroid
      trans = (c0.x - other.centroid.x, c0.y - other.centroid.y)
      other = affine.translate(other, trans[0], trans[1])
      # now see if they align at all
      # symms = self.get_symmetries(other_polygon = other)
      # rot_shifts = symms["rotation-shifts"]
      symms = self.get_symmetries(other_polygon = other)
      rot_shifts = [s.shift for s in symms if s.type == "rotation"]
      if len(rot_shifts) == 0:
        print(f"Polygons do not match at all!")
        return None
      # now align the corners so indexes of their corners match
      offset = rot_shifts[0]
      if offset != 0:
        other = tiling_utils.offset_polygon_corners(other, -offset)
        rotation = -self.get_rotations(rot_shifts)[0].angle #["angles"][0]
        other = tiling_utils.rotate_preserving_order(other, rotation, c0)
      else:
        rotation = 0
      # now determine the additional rotation needed to line them up
      p0_0 = tiling_utils.get_corners(self.polygon)[0]
      other_corners = tiling_utils.get_corners(other, repeat_first = False)
      rotation = rotation + \
        tiling_utils.get_inner_angle(other_corners[0], c0, p0_0)
      # now find store and report the matching symmetries
      result = {}
      symms = self.get_symmetries(other_polygon = other)
      result["offset"] = offset
      result["pre-translation"] = trans
      result["pre-rotation"] = rotation
      result["pre-transform"] = tiling_utils.combine_transforms(
        [tiling_utils.get_translation_transform(trans[0], trans[1]),
         tiling_utils.get_rotation_transform(rotation, (c0.x, c0.y))])
      return result
    else:
      print(f"Polygons have different numbers of sides!")
      return None


  def get_corner_labels(self) -> dict[str, list[str]]:
    """Returns all the reorderings of vertex labels corresponding to each
    symmetry.

    Returns:
      dict[str, list[str]]: A dictionary with two entries. "rotations" is
        a list of labels under the rotation symmetries, "reflections" those
        under the reflection symmetries.
    """
    labels = list(string.ascii_letters.upper())[:self.n]
    return self._get_labels_under_symmetries(labels)


  def get_unique_labels(self, offset:int = 0) -> list[str]:
    labellings = self.get_corner_labels()
    labellings = labellings["rotations"] + labellings["reflections"]
    labellings = ["".join(sorted(x)) for x in zip(*labellings)]
    n_new_labels = len(set(labellings))
    letters = list(string.ascii_letters.upper())[offset:offset + n_new_labels]
    mapping = dict()
    i = 0
    for label in labellings:
      if not label in mapping:
        mapping[label] = letters[i]
        i = i + 1
    return self._get_labels_under_symmetries([mapping[x] for x in labellings])
  
  
  def _get_labels_under_symmetries(
      self, labels:list[str]) -> dict[str, list[str]]:
    cycle = labels * 2
    under_rotation = ["".join(cycle[i:self.n + i]) 
                      for i in self.rotation_shifts]
    under_reflection = ["".join(cycle[self.n + i:i:-1]) 
                        for i in self.reflection_shifts]
    return {"rotations": under_rotation,
            "reflections": under_reflection}


  def plot(self):
    n_subplots = len(self.symmetries)
    
    nr = int(np.round(np.sqrt(n_subplots), 0))
    nc = int(np.ceil(n_subplots / nr))
    n_plots = 0

    fig = pyplot.figure(figsize = (10, 10))
    
    for s in self.symmetries:
      n_plots = n_plots + 1
      ax = fig.add_subplot(nr, nc, n_plots)
      s.plot(ax, self.polygon)

    return ax
