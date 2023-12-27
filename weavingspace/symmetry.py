#!/usr/bin/env python
# coding: utf-8

from typing import Iterable
from dataclasses import dataclass
import numpy as np

import shapely.geometry as geom
import shapely.affinity as affine
import shapely.ops
import weavingspace.tiling_utils as tiling_utils


def _equal_tuples(t1:Iterable[float], t2:Iterable[float]) -> bool:
  """Tests for near equality of two iterables of floats using numpy
  allclose. Wrapped so as to allow find_matches function to take 
  alternative equality tests as input.
  """
  return np.allclose(t1, t2, atol = tiling_utils.RESOLUTION * 10)

def _round_tuple(t:Iterable[float], digits:int = 1) -> Iterable[float]:
  """Convenience function to round all members of an interable. Used for 
  display purposes."""
  return tuple([np.round(x, digits) for x in t])


@dataclass
class Symmetries():

  polygon:geom.Polygon = None
  n:int = None
  rotations:tuple[float] = None
  reflection_axes:tuple[float] = None
  p_code:list[tuple[float]] = None
  p_code_r:list[tuple[float]] = None
  rotation_shifts:list[int] = None
  reflection_shifts:list[int] = None
  
  def __init__(self, polygon:geom.Polygon):
    self.polygon = polygon
    self.n = shapely.count_coordinates(self.polygon) - 1
    self.p_code = self._get_polygon_code(self.polygon)
    self.p_code_r = self._get_polygon_code(self.polygon, mirrored = True)
    syms = self.get_symmetries()
    self.rotation_shifts = syms["rotation-shifts"]
    self.rotations = syms["rotations"]
    self.reflection_shifts = syms["reflection-shifts"]
    self.reflection_axes = syms["reflections"]


  def find_matches(self, seq:Iterable[tuple[float]], 
                  pat:Iterable[tuple[float]],
                  equals = _equal_tuples, debug:bool = False) -> list[int]:
    """ Implements Knuth-Morris-Pratt string pattern matching algorithm. See:
    https://en.wikipedia.org/wiki/Knuth–Morris–Pratt_algorithm which provides
    detailed pseudo-code on which this code is directly based. This
    implementation expects sequences of tuples of floats, although in principle
    any objects could be contained in the Iterables. The default equality
    function tests for close matches, not exact equality.

    Args:
      seq (Iterable[tuple[float]]): The sequence of float tuples to find   
        matches in.
      pat (Iterable[tuple[float]]): The sequence to match.
      equals (function, optional): A function to use to test for equality of
        elements in patterns. Defaults to _equal_tuples().

    Returns:
        Iterable[int]: _description_
    """
    if debug:
      print(f"""
        Searching for: {[_round_tuple(t) for t in pat]} in
        {[_round_tuple(t) for t in seq]}""")
    j, k, = 0, 0
    finds = []
    table = self.get_table(pat, equals)
    while j < len(seq):
      if equals(pat[k], seq[j]):
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


  def get_table(self, pattern:Iterable,
                equals = _equal_tuples) -> Iterable[int]:
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
        while cnd >= 0 and not equals(pattern[pos], pattern[cnd]):
          cnd = T[cnd]
      pos = pos + 1
      cnd = cnd + 1
    
    T[pos] = cnd
    return tuple(T.values())


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
    Pattern Recognition, ed. GT Toussaint, 6:41–51. Computational Morphology. 
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


  def get_symmetries(self, other:geom.Polygon = None) -> dict[str, list[float]]:
    """Detects rotation and reflection symmetries of the supplied polygon.
    Based on
    
      Eades P. 1988. Symmetry Finding Algorithms. In Machine Intelligence and 
      Pattern Recognition, ed. GT Toussaint, 6:41–51. Computational Morphology. 
      North-Holland. doi: 10.1016/B978-0-444-70467-2.50009-6.
    
    and also
    
      Wolter JD, TC Woo, and RA Volz. 1985. Optimal algorithms for symmetry 
      detection in two and three dimensions. The Visual Computer 1(1): 37–48. 
      doi: 10.1007/BF01901268.
    
    Details in these papers are not unambiguous. This implementation was
    developed based on them, but with a lot of trial and error to get the index
    offsets and (especially) retrieval of the reflection axes angles to work.
    This implementation appears to work and is cleaner than either of those 
    descriptions!

    Args:
      p (geom.Polygon): The polygon for which symmetries are required.

    Returns:
      dict[str, list[float]]: Dictionary with entries:
        'rotations', a list of the angles of rotational symmetry. This will 
          always contain at least [0] the 'identity' symmetry. 
        'reflections', a list of the angles to the x-axis of lines of reflection
          symmetry. This list may be empty if no such lines exist.
    """
    # compose extended sequence of length-angle pairs for cyclic matching
    S = self.p_code + self.p_code[:-1]
    c = self.polygon.centroid
    # get code for the 'other' polygon, if it exists...
    if not other is None:
      p_code = self._get_polygon_code(other)
      p_code_r = self._get_polygon_code(other, mirrored = True)
    else:
      p_code = self.p_code
      p_code_r = self.p_code_r
    # ==== Rotations ====
    rotation_matches = self.find_matches(S, p_code)
    # ==== Reflections ====
    reflection_matches = self.find_matches(S, p_code_r)
    # calculate rotations
    rot_angles = [i * 360 / self.n for i in rotation_matches] 
    rotations = [self._get_rotation_transform(a, (c.x, c.y)) 
                 for a in rot_angles]
    # calculate all possible reflection axes - ordering of these is important
    # for correct picking - don't mess with this code!
    bearings = tiling_utils.get_side_bearings(self.polygon)
    reflection_axes = []
    interior_angles = tiling_utils.get_interior_angles(self.polygon)
    for i in range(self.n):
      # angle bisectors
      reflection_axes.append(bearings[i] - interior_angles[i] / 2)
      # edge_perpendicular bisectors
      reflection_axes.append(bearings[i] - 90)  # normal to the edge
    # pick out those where matches occurred
    ref_angles = [reflection_axes[i] for i in reflection_matches]
    reflections = [self._get_reflection_transform(a, (c.x, c.y))
                   for a in ref_angles]
    return {
      "rotation-shifts": rotation_matches,
      "rotations": rot_angles,
      "rotation-transforms": rotations,
      "reflection-shifts": reflection_matches,
      "reflections": ref_angles,
      "reflection-transforms": reflections}


  def get_matching_transforms(self, polygon:geom.Polygon):
    if self.n == shapely.count_coordinates(polygon) - 1:
      c0 = self.polygon.centroid
      c2 = polygon.centroid
      translation = (c0.x - c2.x, c0.y - c2.y)
      comparison_poly = affine.translate(polygon, translation[0], translation[1])
      p0_0 = tiling_utils.get_corners(self.polygon)[0]
      # this next move because affine transforms can sometimes reorder corners
      corners = tiling_utils.get_corners(comparison_poly, repeat_first = False)
      p2_0 = corners[0]
      rotation = (tiling_utils.get_inner_angle(p2_0, c0, p0_0))
      comparison_poly = geom.Polygon([
        affine.rotate(c, rotation, c0) for c in corners])
      syms = self.get_symmetries(other =  comparison_poly)
      syms["pre-translation"] = \
        self._get_translation_transform(translation[0], translation[1])
      syms["pre-rotation"] = \
        self._get_rotation_transform(rotation, (c0.x, c0.y))
      return syms
    else:
      print(f"Polygons have different numbers of sides!")
      return None


  def _get_translation_transform(self, dx:float, dy:float):
    return [1, 0, 0, 1, dx, dy]
  

  def _get_rotation_transform(self, angle:float, 
                              centre:tuple[float] = None):
    if centre is None or np.allclose((0, 0), centre, 
                                     atol = tiling_utils.RESOLUTION):
      a = np.radians(angle)
      return [np.cos(a), -np.sin(a), np.sin(a), np.cos(a), 0, 0]
    t1 = self._get_translation_transform(-centre[0], -centre[1])
    r = self._get_rotation_transform(angle)
    t2 = self._get_translation_transform(centre[0], centre[1])
    return self._combine_transforms([t1, r, t2])
  
  
  def _get_reflection_transform(self, angle:float,
                                centre:tuple[float] = None):
    if centre is None or np.allclose((0, 0), centre, 
                                     atol = tiling_utils.RESOLUTION):
      A = 2 * np.radians(angle)
      return (np.cos(A), np.sin(A), np.sin(A), -np.cos(A), 0, 0)
    r = self._get_reflection_transform(angle)
    t1 = self._get_translation_transform(-centre[0], -centre[1])
    t2 = self._get_translation_transform(centre[0], centre[1])
    return self._combine_transforms([t1, r, t2])


  def _combine_transforms(self, transforms:list[list[float]]) -> list[float]:
    result = np.identity(3)
    for t in transforms:
      result = self._as_numpy_matrix(t) @ result
    return self._as_shapely_transform(result)


  def _reverse_transform(self, transform:list[float]) -> list[float]:
    return self._as_shapely_transform(
      np.linalg.inv(self._as_numpy_matrix(transform)))


  def _as_shapely_transform(self, arr:np.array) -> list[float]:
    return [arr[0][0], arr[0][1], arr[1][0], arr[1][1], arr[0][2], arr[1][2]]
  
  
  def _as_numpy_matrix(self, tr:list[float]) -> np.array:
    return np.array([[tr[0], tr[1], tr[4]],
                     [tr[2], tr[3], tr[5]],
                     [    0,     0,     1]])

