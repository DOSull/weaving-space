#!/usr/bin/env python
# coding: utf-8

from typing import Iterable
from dataclasses import dataclass
import string
import numpy as np
import matplotlib.pyplot as pyplot

import shapely.geometry as geom
import shapely.affinity as affine
import shapely.ops
import geopandas as gpd

import weavingspace.tiling_utils as tiling_utils


def _equal_tuples(t1:Iterable[float], t2:Iterable[float]) -> bool:
  """Tests for near equality of two iterables of floats using numpy
  allclose. Wrapped so as to allow find_matches function to take 
  alternative equality tests as input.
  """
  return np.allclose(t1, t2, atol = tiling_utils.RESOLUTION * 10)

def _round_tuple(t:Iterable[float], digits:int = 3) -> Iterable[float]:
  """Convenience function to round all members of an interable. Used for 
  display purposes."""
  return tuple([np.round(x, digits) for x in t])


@dataclass
class Symmetries():

  polygon:geom.Polygon = None
  n:int = None
  p_code:list[tuple[float]] = None
  p_code_r:list[tuple[float]] = None
  rotation_shifts:list[int] = None
  reflection_shifts:list[int] = None
  

  def __init__(self, polygon:geom.Polygon):
    self.polygon = polygon
    self.n = shapely.count_coordinates(self.polygon) - 1
    self.p_code = self._get_polygon_code(self.polygon)
    self.p_code_r = self._get_polygon_code(self.polygon, mirrored = True)
    symmetries = self._find_symmetries()
    self.rotation_shifts = symmetries["rotation-shifts"]
    self.reflection_shifts = symmetries["reflection-shifts"]
    

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


  def _find_symmetries(
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
      polygon will that for which this Symmetries object has been initialised.

    Returns:
      dict[str, list[float]]: _description_
    """
    # compose extended sequence of length-angle pairs for cyclic matching
    S = self.p_code + self.p_code[:-1]
    # get code for the 'other' polygon, if it exists...
    if not other_polygon is None:
      p_code = self._get_polygon_code(other_polygon)
      p_code_r = self._get_polygon_code(other_polygon, mirrored = True)
    else:
      p_code = self.p_code
      p_code_r = self.p_code_r
    # ==== Rotations ====
    rotation_shifts = self.find_matches(S, p_code)
    # ==== Reflections ====
    reflection_shifts = self.find_matches(S, p_code_r)
    return {"rotation-shifts": rotation_shifts,
            "reflection-shifts": reflection_shifts}


  def get_rotations(self, offsets:list[int] = None) -> dict[str, list[float]]:
    """Gets the rotations associated with this collection of symmetries.

    Returns:
      dict[str, list[float], list[float]]: Two item dictionary with 'angles' and
        'transforms' entries, each a list of floats.
    """
    if offsets is None:
      offsets = self.rotation_shifts
    c = self.polygon.centroid
    rot_angles = [np.round(i * 360 / self.n, 6) for i in offsets] 
    rotations = [np.round(self._get_rotation_transform(a, (c.x, c.y)), 6) 
                 for a in rot_angles]
    return {"angles": rot_angles, "transforms": rotations}
    

  def get_reflections(self, offsets:list[int] = None) -> dict[str, list[float]]:
    """Gets the reflections associated with this collection of symmetries.

    Returns:
      dict[str, list[float], list[float]]: Two item dictionary with 'angles' and
        'transforms' entries, each a list of floats.
    """
    if offsets is None:
      offsets = self.reflection_shifts
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
    reflections = [np.round(self._get_reflection_transform(a, (c.x, c.y)), 6)
                   for a in ref_angles]
    return {"angles": ref_angles, "transforms": reflections}


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
      symms = self._find_symmetries(other_polygon = other)
      rot_shifts = symms["rotation-shifts"]
      if len(rot_shifts) == 0:
        print(f"Polygons do not match at all!")
        return None
      # now align the corners so indexes of their corners match
      offset = rot_shifts[0]
      if offset != 0:
        other = tiling_utils.offset_polygon_corners(other, -offset)
        rotation = -self.get_rotations(rot_shifts)["angles"][0]
        other = tiling_utils.rotate_preserving_order(other, rotation, c0)
      else:
        rotation = 0
      # now determine the additional rotation needed to line them up
      p0_0 = tiling_utils.get_corners(self.polygon)[0]
      other_corners = tiling_utils.get_corners(other, repeat_first = False)
      rotation = rotation + \
        tiling_utils.get_inner_angle(other_corners[0], c0, p0_0)
      # now find store and report the matching symmetries
      symms = self._find_symmetries(other_polygon = other)
      rot_shifts = symms["rotation-shifts"]
      rot_details = self.get_rotations(rot_shifts)
      symms["offset"] = offset
      symms["rotation-angles"] = rot_details["angles"]
      symms["rotation-transforms"] = rot_details["transforms"]
      ref_shifts = symms["reflection-shifts"]
      ref_details = self.get_reflections(ref_shifts)
      symms["reflection-angles"] = ref_details["angles"]
      symms["reflection-transforms"] = ref_details["transforms"]
      symms["pre-translation"] = trans
      symms["pre-rotation"] = rotation
      symms["pre-transform"] = self._combine_transforms(
        [self._get_translation_transform(trans[0], trans[1]),
         self._get_rotation_transform(rotation, (c0.x, c0.y))])
      return symms
    else:
      print(f"Polygons have different numbers of sides!")
      return None


  def _get_translation_transform(self, dx:float, dy:float) -> list[float]:
    """Returns the shapely affine transform tuple for a translation.

    Args:
        dx (float): translation distance in x direction.
        dy (float): translation distance in y direction.

    Returns:
      list[float]: a six item list of floats, per the shapely.affinity.
      affine_transform method, see 
        https://shapely.readthedocs.io/en/stable/manual.html#affine-transformations
    """
    return [1, 0, 0, 1, dx, dy]
  

  def _get_rotation_transform(self, angle:float, 
                              centre:tuple[float] = None) -> list[float]:
    """Returns the shapely affine transform tuple for a rotation, optionally
    about a supplied centre point.

    Args:
        angle (float): the angle of rotation (in degrees).
        centre (tuple[float], optional): An option centre location. Defaults to 
        None, which will in turn be converted to (0, 0).

    Returns:
      list[float]: a six item list of floats, per the shapely.affinity.
        affine_transform method, see 
          https://shapely.readthedocs.io/en/stable/manual.html#affine-transformations
    """
    if centre is None or np.allclose((0, 0), centre, 
                                     atol = tiling_utils.RESOLUTION):
      a = np.radians(angle)
      return [np.cos(a), -np.sin(a), np.sin(a), np.cos(a), 0, 0]
    
    dx, dy = centre
    t1 = self._get_translation_transform(-dx, -dy)
    r = self._get_rotation_transform(angle)
    t2 = self._get_translation_transform(dx, dy)
    return self._combine_transforms([t1, r, t2])
  
  
  def _get_reflection_transform(
      self, angle:float, centre:tuple[float] = None) -> list[float]:
    """Returns a shapely affine transform tuple that will reflect a shape
    in a line at the specified angle, optionally through a specified centre
    point.

    Args:
      angle (float): angle to the x-axis of the line of reflection.
      centre (tuple[float], optional): point through which the line of 
        reflection passes. Defaults to None, which
        will in turn be converted to (0, 0).

    Returns:
      list[float]: a six item list of floats, per the shapely.affinity.
        affine_transform method, see 
          https://shapely.readthedocs.io/en/stable/manual.html#affine-transformations
    """
    if centre is None or np.allclose((0, 0), centre, 
                                     atol = tiling_utils.RESOLUTION):
      A = 2 * np.radians(angle)
      return (np.cos(A), np.sin(A), np.sin(A), -np.cos(A), 0, 0)
    dx, dy = centre
    t1 = self._get_translation_transform(-dx, -dy)
    r = self._get_reflection_transform(angle)
    t2 = self._get_translation_transform(dx, dy)
    return self._combine_transforms([t1, r, t2])


  def _combine_transforms(self, transforms:list[list[float]]) -> list[float]:
    """Returns a shapely affine transform list that combines the listed
    sequence of transforms applied in order.

    Args:
      transforms (list[list[float]]): sequence of transforms to combine.

    Returns:
      list[float]: a transform tuple combining the supplied transforms applied
        in order, see 
          https://shapely.readthedocs.io/en/stable/manual.html#affine-transformations
    """
    result = np.identity(3)
    for t in transforms:
      result = self._as_numpy_matrix(t) @ result
    return self._as_shapely_transform(result)


  def _reverse_transform(self, transform:list[float]) -> list[float]:
    """Returns the inverse shapely affine transform of the supplied transform.

    Args:
      transform (list[float]): the transform for which the inverse is desired.

    Returns:
      list[float]: shapely affine transform tuple that will invert the supplied
        transform.
    """
    return self._as_shapely_transform(
      np.linalg.inv(self._as_numpy_matrix(transform)))


  def _as_shapely_transform(self, arr:np.array) -> list[float]:
    """Returns the shapely affine transform list equivalent to the supplied
    numpy matrix of a conventional augmented affine transform matrix.

    Args:
      arr (np.array): augmented affine transform matrix of the desired 
        transform.

    Returns:
      list[float]: desired shapely affine transform list of floats.
    """
    return [arr[0][0], arr[0][1], arr[1][0], arr[1][1], arr[0][2], arr[1][2]]
  
  
  def _as_numpy_matrix(self, transform:list[float]) -> np.array:
    """Converts the supplied shapely affine transform list to an augmented
    affine transform matrix in numpy array form. This makes combining transforms
    much easier.

    Args:
      transform (list[float]): the transform in shapely format.

    Returns:
        np.array: the transform in numpy matrix format.
    """
    return np.array([[transform[0], transform[1], transform[4]],
                     [transform[2], transform[3], transform[5]],
                     [           0,            0,            1]])


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


  def plot(self, which:str = "both", all_in_one = False):
    
    if all_in_one:
      print("all_in_one option not yet implemented!")
    
    if which == "both":
      n_subplots = len(self.rotation_shifts) + len(self.reflection_shifts)
    elif which[:3] == "rot":
      n_subplots = len(self.rotation_shifts)
    elif which[:3] == "ref":
      n_subplots = len(self.reflection_shifts)
    else:
      print(f"Please specify symmetries to plot using the which parameter!")
      return
    
    nr = int(np.floor(np.sqrt(n_subplots)))
    nc = int(np.ceil(n_subplots / nr))
    nplots = 0
    mrr = shapely.oriented_envelope(self.polygon)
    w, h = sorted(tiling_utils.get_side_lengths(mrr)[:2])

    fig = pyplot.figure(figsize = (10, 10))
    labels = self.get_unique_labels()
    c = self.polygon.centroid

    if which == "both" or which[:3] == "rot":
      rotations = self.get_rotations()["angles"]
      for (i, rotn), label in zip(enumerate(rotations), labels["rotations"]):
        nplots = nplots + 1
        ax = fig.add_subplot(nr, nc, nplots)
        gpd.GeoSeries([self.polygon]).plot(ax = ax, fc = "lightgrey", lw = 0)
        for lbl, pt in zip(list(label), self.polygon.exterior.coords):
          ax.annotate(lbl, xy = pt, ha = "center", va = "center", fontsize = 14)
        if i > 0:
          rotn = np.radians(rotn) 
          arc = geom.LineString(
            [[w/6 * np.cos(a), w/6 * np.sin(a)]
             for a in np.linspace(0, rotn, 50)])
          angle = geom.LineString(
            [(w/3, 0), (0, 0), 
             (w/4 * np.cos(rotn), w/4 * np.sin(rotn))])
          ls = [affine.translate(l, c.x, c.y) for l in [arc, angle]]
          gpd.GeoSeries(ls).plot(ax = ax, ec = "b", lw = 0.35)
        self.get_corner_labels()
        pyplot.axis("off")
        
    if which == "both" or which[:3] == "ref":
      reflections = self.get_reflections()["angles"]
      for refn, label in zip(reflections, labels["reflections"]):
        nplots = nplots + 1
        ax = fig.add_subplot(nr, nc, nplots)
        gpd.GeoSeries([self.polygon]).plot(ax = ax, fc = "lightgrey", lw = 0)
        for lbl, pt in zip(list(label), self.polygon.exterior.coords):
          ax.annotate(lbl, xy = pt, ha = "center", va = "center", fontsize = 14)
        ls = geom.LineString([(-h/2 * 1.1, 0), (h/2 * 1.1, 0)])
        ls = affine.rotate(ls, refn)
        ls = affine.translate(ls, c.x, c.y)
        gpd.GeoSeries([ls]).plot(ax = ax, ec = "r", ls = "dashed", lw = 0.5)
        pyplot.axis("off")

    return ax


  def find_matches(self, seq:Iterable[tuple[float]], 
                   pat:Iterable[tuple[float]],
                   equals = _equal_tuples) -> list[int]:
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


