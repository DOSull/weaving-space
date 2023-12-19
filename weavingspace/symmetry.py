#!/usr/bin/env python
# coding: utf-8

from typing import Iterable
from dataclasses import dataclass
import numpy as np

import shapely.geometry as geom
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
  rotations:tuple[float] = None
  reflection_axes:tuple[float] = None
  
  def __init__(self, polygon:geom.Polygon):
    self.polygon = polygon
    syms = self.get_symmetries(self.polygon)
    self.rotations = syms["rotations"]
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


  def get_symmetries(self, 
                      p:geom.Polygon) -> tuple[tuple[float], tuple[float]]:
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
    developed based on them, but with some trial and error to get the index
    offsets and (especially) the retrieval of the reflection axes angles
    right. This implementation appears to work and cleaner than either of 
    those descriptions.

    Args:
      p (geom.Polygon): The polygon for which symmetries are required.

    Returns:
      dict[str, list[float]]: Dictionary with entries:
        'rotations', a list of the angles of rotational symmetry. This will 
        always contain at least [0] the 'non-rotation' symmetry. 
        'reflections', a list of the angles to the x-axis of lines of reflection
        symmetry. This list may be empty if no such lines exist.
        'mirrors', a list of all the potential reflection symmetry axes. There 
        will 2 x n of these, bisectors of the interior angles and 
        perpendiculars to each edge. Primarily for debugging purposes (if any 
        are actual lines of reflection symmetry they will be returned in the 'reflections' entry). 
    """
    # get edge lengths and angles
    lengths = tiling_utils.get_edge_lengths(p)
    raw_angles = tiling_utils.get_interior_angles(p)
    # the forward (i.e. not-mirrored) polygon pairs edge lengths and the angle
    # between an edge and its successor going CW around the polygon, hence this
    # slice operation. See Eades in particular for clarification.
    #
    #    A[i] ----L[i]---- A (i+1)
    #   /                   \
    #  /                     \ 
    #                          
    # i.e. edge i is between angles i and i + 1, and angle i is between edge
    # i-1 and i. The unmirrored encoding pairs length[i] with angle[i+1]
    angles = raw_angles[1:] + raw_angles[:1]
    # encode polygon as length-angle pairs
    poly_code = list(zip(lengths, angles))
    # ==== Rotations ====
    # compose extended sequence of length-angle pairs for cyclic matching
    S = poly_code + poly_code[:-1]
    rotation_matches = self.find_matches(S, poly_code)
    # ==== Reflections ====
    # encode mirrored polygon: note that here the indexes of lengths and angles
    # in the original polygon are matched. This means the angle is still the one
    # between the edge and its successor now proceeding CCW around the polygon.
    lengths_r = list(reversed(lengths))
    angles_r = list(reversed(raw_angles))
    poly_r_code = list(zip(lengths_r, angles_r))
    reflection_matches = self.find_matches(S, poly_r_code)
    # calculate rotations
    n = len(lengths)
    rotations = [i * 360 / n for i in rotation_matches]
    # calculate all possible reflection axes - ordering of these is important
    # for correct picking - don't mess with this code!
    bearings = tiling_utils.get_edge_bearings(p)
    reflection_axes = []
    for i in range(n):
      reflection_axes.append(bearings[i] - raw_angles[i] / 2)  # angle bisector
      reflection_axes.append(bearings[i] - 90)  # normal to the edge
    # pick out those where matches occured
    reflections = [reflection_axes[i] for i in reflection_matches]
    return {"rotations": rotations,
            "reflections": reflections,
            "mirrors": _round_tuple(reflection_axes)}
