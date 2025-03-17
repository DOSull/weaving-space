#!/usr/bin/env python
# coding: utf-8

import inspect
from typing import Iterable
from typing import Any
from typing import Union
from dataclasses import dataclass
import string
import copy

import numpy as np
import matplotlib.pyplot as pyplot
import io
import PIL

import shapely.geometry as geom
import shapely.affinity as affine
import shapely.ops
import geopandas as gpd

# import weavingspace.tiling_utils as tiling_utils
from weavingspace import tiling_utils


class KMP_Matcher:
  """Class to find matching subsequences in a sequence."""

  sequence: Iterable
  """Iterable in which subsequences are to be found."""

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

class Transform:
  """Class to store details of a transform and draw it."""
  transform_type: str
  """Type of transform, 'rotation', 'reflection', 'translation' or 'identity'."""
  angle: float
  """Angle of rotation (degrees)."""
  centre: geom.Point
  """Centre of the transformation."""
  translation: tuple[float] 
  """X and Y coordinates shifts of a translation transform. A glide reflection
  may also include this."""
  transform: tuple[float]
  """Six element tuple for the transform in shapely.transform format. See
  https://shapely.readthedocs.io/en/stable/manual.html#affine-transformations and
  methods in `weavingspace.tiling_utils`."""
  offset: int

  def __init__(self, transform_type:str, angle:float, centre:geom.Point, 
               translation: tuple[float], transform:tuple[float]):
    self.transform_type = transform_type
    self.angle = angle
    self.centre = centre
    self.translation = translation
    self.transform = transform

  def __str__(self) -> str:
    if self.transform_type == "rotation":
      return f"{self.transform_type} {self.angle:.1f}° POINT ({self.centre.x:.1f} {self.centre.y:.1f}) {tuple([np.round(x, 3) for x in self.transform])}"
    elif self.transform_type == "translation":
      return f"{self.transform_type} {self.translation} {tuple([np.round(x, 3) for x in self.transform])}"
    elif self.transform_type == "reflection":
      return f"{self.transform_type} {self.angle:.1f}° {tuple([np.round(x, 3) for x in self.transform])}"
    else:
      return f"{self.transform_type} {tuple([np.round(x, 3) for x in self.transform])}"

  def __repr__(self) -> str:
    return str(self)

  def apply(self, geometry:Any) -> Any:
    """Applies this transform to supplied shapely geometry and returns result.

    Args:
      geometry (Any): a shapely geometry to transform.
    
    Returns:
      Any: the resulting transformed shapely geometry.
    """
    return affine.affine_transform(geometry, self.transform)

  def draw(self, ax:pyplot.axes, **kwargs) -> pyplot.Axes:
    """Draws this transform on the supplied axes. Arguments specific to each
    transform type are supplied as **kwargs and are documented in 
    `draw_rotation`, `draw_reflection`, and `draw_translation`.

    Args:
        ax (pyplot.axes): the axes on which to draw the transform.

    Returns:
        pyplot.Axes: the axes with the rendering of this transform added.
    """
    if self.transform_type == "rotation":
      rotn_args = list(inspect.signature(self.draw_rotation).parameters)
      rotn_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in rotn_args}
      return self.draw_rotation(ax = ax, **rotn_dict)
    elif self.transform_type == "reflection":
      refn_args = list(inspect.signature(self.draw_reflection).parameters)
      refn_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in refn_args}
      return self.draw_reflection(ax = ax, **refn_dict)
    if self.transform_type == "translation":
      trans_args = list(inspect.signature(self.draw_translation).parameters)
      trans_dict = {k: kwargs.pop(k) for k in dict(kwargs) if k in trans_args}
      return self.draw_translation(ax = ax, **trans_dict)
    elif self.transform_type == "identity":
      w = ax.get_window_extent().width
      return ax.set_title(f"{self.transform_type}", fontsize = w / 20)

  def draw_rotation(self, ax:pyplot.Axes, 
                    radius = 200, add_title = True) -> pyplot.Axes:
    x, y = self.centre.x, self.centre.y
    axis = geom.LineString([(x, y), (x + radius * 1.25, y)])
    arc = geom.LineString([
      geom.Point(x + radius * np.cos(np.radians(a)),
                  y + radius * np.sin(np.radians(a)))
      for a in np.linspace(0, self.angle, 50)])
    gpd.GeoSeries([self.centre]).plot(
      ax = ax, color = "r", markersize = 4, zorder = 5)
    gpd.GeoSeries([axis, arc]).plot(ax = ax, color = "r", lw = .5)
    if add_title:
      w = ax.get_window_extent().width
      ax.set_title(f"{self.transform_type} {np.round(self.angle, 2)}",
                   fontsize = w / 20)
    return ax

  def draw_reflection(self, ax:pyplot.Axes, w:float = 5, 
                      mirror_length = 100, add_title:bool = True) -> pyplot.Axes:
    x, y = self.centre.x, self.centre.y
    dx, dy = self.translation
    r = np.sqrt(self.translation[0] ** 2 + self.translation[1] ** 2) 
    mirror = geom.LineString([
      (x - mirror_length / 2 * np.cos(np.radians(self.angle)), 
        y - mirror_length / 2 * np.sin(np.radians(self.angle))),
      (x + mirror_length / 2 * np.cos(np.radians(self.angle)), 
        y + mirror_length / 2 * np.sin(np.radians(self.angle)))])
    gpd.GeoSeries([mirror]).plot(
      ax = ax, color = "r", lw = 1, ls = "dashdot", zorder = 5)
    no_slide = np.isclose(r, 0, rtol = 1e-6, atol = 1e-6)
    if not no_slide:
      pyplot.arrow(
        x - dx / 2, y - dy / 2, dx, dy, length_includes_head = True,
        width = w, fc = "k", ec = None, head_width = w * 6, zorder = 5)
    if add_title:
      w = ax.get_window_extent().width
      ax.set_title(f"{self.transform_type} {np.round(self.angle, 2)}",
                   fontsize = w / 20)
    return ax

  def draw_translation(self, ax:pyplot.Axes, c:geom.Point, 
                       w:float = 5, add_title:bool = True) -> pyplot.Axes:
    gpd.GeoSeries([c]).plot(ax = ax, color = "b")
    pyplot.arrow(c.x, c.y, self.translation[0], self.translation[1], lw = 0.5,
                width = w, fc = "b", ec = None, head_width = w * 6, zorder = 5)
    if add_title:
      w = ax.get_window_extent().width
      ax.set_title(f"{self.transform_type} ({self.translation[0]:.1f}, {self.translation[1]:.1f})",
                   fontsize = w / 20)
    return ax


@dataclass
class Symmetries():
  """Class to identify and store the symmetries of a supplied shapely.Polygon
  as a list of `Transform` objects.
  """
  polygon:geom.Polygon = None
  """Polygon from which these symmetries are derived."""
  matcher:KMP_Matcher = None
  """the subsequence matcher used in symmetry detection."""
  n:int = None
  """number of vertices of the polygon."""
  p_code:list[tuple[float]] = None
  """the encoding of the polygon as sequence of length, angle pairs."""
  p_code_r:list[tuple[float]] = None
  """the reversed encoding used to detect reflection symmetries."""
  rotation_shifts:list[int] = None
  """list of number of 2pi/n rotation symmetries."""
  reflection_shifts:list[int] = None
  """list of pi/n relection angle symmetries."""
  symmetries:list[Transform] = None
  """list of Transform objects with more complete information."""
  symmetry_group:str = None
  """the code denoting the symmetry group"""

  def __init__(self, polygon:geom.Polygon):
    self.polygon = polygon
    self.n = shapely.count_coordinates(self.polygon) - 1
    self.p_code = self._get_polygon_code(self.polygon)
    self.p_code_r = self._get_polygon_code(self.polygon, mirrored = True)
    self.matcher = KMP_Matcher(self.p_code + self.p_code[:-1])
    self.symmetries = self.get_symmetries()
    self.symmetry_group = self.get_symmetry_group_code()

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
    i-1 and i. The mirrored encoding pairs length[i] with angle[i+1].
    
    The mirrored encoding pairs matched indexes of lengths and angles in the
    original polygon. This means the angle is still the one between and edge 
    and its successor but proceeding CCW around the polygon.
    
    See in particular for clarification:
    
      Eades P. 1988. Symmetry Finding Algorithms. In Machine Intelligence and 
      Pattern Recognition, ed. GT Toussaint, 6:41-51. Computational Morphology. 
      North-Holland. doi: 10.1016/B978-0-444-70467-2.50009-6.
    
    Args:
      p (geom.Polygon): the polygon to encode
      mirrored (bool, optional): if true encoding will be in CCW order. 
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

  def get_symmetries(self) -> list[Transform]:
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

    Returns:
      list[Transform]: a list of Transform objects representing the polygon
        symmetries.
    """
    # ==== Rotations ====
    self.rotation_shifts = self.matcher.find_matches(self.p_code)
    # ==== Reflections ====
    self.reflection_shifts = self.matcher.find_matches(self.p_code_r)
    return self.get_rotations(self.rotation_shifts) + \
           self.get_reflections(self.reflection_shifts)

  def get_rotations(self, offsets:list[int]) -> list[Transform]:
    """Gets the rotations associated with this collection of symmetries.

    Returns:
      list[Transform]: A list of the rotation symmetries associated with this
        polygon.
    """
    c = self.polygon.centroid
    rot_angles = [np.round(i * 360 / self.n, 6) for i in offsets] 
    rotation_transforms = [
      np.round(tiling_utils.get_rotation_transform(a, (c.x, c.y)), 6) 
      for a in rot_angles]
    return [Transform("rotation", angle, c, (0, 0), transform)
            for transform, angle in zip(rotation_transforms, rot_angles)]

  def get_reflections(self, offsets:list[int]) -> list[Transform]:
    """Gets the reflections associated with this collection of symmetries.

    Returns:
      list[Transform]: A list of the reflection symmetries associated with this
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
    reflection_transforms = [
      np.round(tiling_utils.get_reflection_transform(a, (c.x, c.y)), 6)
      for a in ref_angles]
    return [Transform("reflection", angle, c, (0, 0), transform)
            for transform, angle in zip(reflection_transforms, ref_angles)]

  def get_symmetry_group_code(self):
    if len(self.reflection_shifts) == 0:
      return f"C{len(self.rotation_shifts)}"
    else:
      return f"D{len(self.rotation_shifts)}"

  def get_corner_offset(self, poly2:geom.Polygon):
    s2 = Symmetries(poly2)
    matches = self.matcher.find_matches(s2.p_code)
    if len(matches) > 0:
      return matches[0]
    matches = self.matcher.find_matches(s2.p_code_r)
    if len(matches) > 0:
      return -matches[0]
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

  def plot(self, as_image:bool = False, title:str = ""):
    fig = pyplot.figure()
    fig.suptitle(title)

    n_subplots = len(self.symmetries)
    nr = int(np.ceil(np.sqrt(n_subplots)))
    nc = int(np.ceil(n_subplots / nr))
    n_plots = 0

    for s in self.symmetries:
      n_plots = n_plots + 1
      ax = fig.add_subplot(nr, nc, n_plots)
      s.plot(ax, self.polygon)
    
    if as_image:
      buf = io.BytesIO()
      fig.savefig(buf)
      buf.seek(0)
      return PIL.Image.open(buf)
    return ax


class Shape_Matcher:
  shape: geom.Polygon
  s1: Symmetries
  other: geom.Polygon
  centre: geom.Point
  translation: tuple[float]
  matches: list[Transform]
  identity_transform: tuple[float] = (1, 0, 0, 1, 0, 0)
  
  def __init__(self, shape: geom.Polygon):
    self.shape = shape
    self.s1 = Symmetries(shape)

  def get_polygon_matches(self, shape2: geom.Polygon):
    s2 = Symmetries(shape2)
    if self.s1.symmetry_group != s2.symmetry_group:
      # print("No matches")
      return None
    else:
      match_rot, rots = self._get_rotation_matches(s2)
      match_ref, refs = self._get_reflection_matches(s2)
      if match_rot and match_ref:
        # print(f"Rotation and reflection matches found")
        return rots + refs
      elif match_rot:
        # print(f"Only rotation matches found")
        return rots
      elif match_ref:
        # print(f"Only reflection matches found")
        return refs
      else:
        # print (f"No matches found")
        return None

  def _get_rotation_matches(self, s2:Symmetries):
    matches = self.s1.matcher.find_matches(s2.p_code)
    if len(matches) == 0:
      return False, None
    else:
      transforms = []
      # get lists of polygon corners aligned correctly to measure the angle
      p1_corners = tiling_utils.get_corners(self.shape, repeat_first = False)
      p2_corners = tiling_utils.get_corners(s2.polygon, repeat_first = False)
      for m in matches:
        a, c = self.get_angle_between_polygons(p1_corners, p2_corners, m)
        if a == "translation":
          # print("One of the rotation matches is a translation.")
          transforms.append(
            Transform("translation", None, None, c, (1, 0, 0, 1, c[0], c[1])))
        elif a == "identity":
          transforms.append(
            Transform("identity", None, None, c, (1, 0, 0, 1, 0, 0)))
        elif a is None:
          continue
        else:
          transforms.append(Transform("rotation", a, c, (0, 0),
            tiling_utils.get_rotation_transform(a, (c.x, c.y))))
      return True, transforms

  def get_angle_between_polygons(self, corners1:list[geom.Point], 
      corners2:list[geom.Point], offset:int) -> tuple[Union[float,geom.Point]]:
    corners2 = corners2[offset:] + corners2[:offset]
    dists = [p1.distance(p2) for p1, p2 in zip(corners1, corners2)]
    if all([np.isclose(dists[0], d, atol = 1e-3, rtol = 1e-3) for d in dists]):
      if np.isclose(dists[0], 0, atol = 1e-3, rtol  = 1e-3):
        return "identity", (0, 0)
      else:
        return "translation", (corners1[0].x - corners2[0].x,
                               corners1[0].y - corners2[0].y)
    ordered_dists = sorted([(i, d) for i, d in enumerate(dists)], 
                           key = lambda x: x[1])
    AB = ordered_dists[-2:]
    p1A, p1B = corners1[AB[0][0]], corners1[AB[1][0]]
    p2A, p2B = corners2[AB[0][0]], corners2[AB[1][0]]
    perpAA = tiling_utils.get_straight_line(p1A, p2A, True)
    perpBB = tiling_utils.get_straight_line(p1B, p2B, True)
    centre = tiling_utils.get_intersection(perpAA, perpBB)
    if centre is None:
      angle = None
    else:
      angle = -tiling_utils.get_inner_angle(p1A, centre, p2A)
      if angle < -179:
        angle = angle + 360
      elif angle > 180:
        angle = angle - 360
    return angle, centre

  def _get_reflection_matches(self, s2:Symmetries):
    ctr1 = self.s1.polygon.centroid
    ctr2 = s2.polygon.centroid
    c = ((ctr1.x + ctr2.x) / 2, (ctr1.y + ctr2.y) / 2)
    matches = self.s1.matcher.find_matches(s2.p_code_r)
    if len(matches) == 0:
      return False, None
    reflections1 = self.s1.get_reflections(matches)
    reflections2 = s2.get_reflections(matches)
    angles = [(ref1.angle + ref2.angle) / 2
              for ref1, ref2 in zip(reflections1, reflections2)]
    trs = [tiling_utils.get_reflection_transform(angle, c)
           for angle in angles]
    ctrs2r = [affine.affine_transform(ctr2, tr) for tr in trs]
    dxdys = [(ctr1.x - ctr2r.x, ctr1.y - ctr2r.y) for ctr2r in ctrs2r]
    return True, [Transform(
      "reflection", angle, geom.Point(c), dxdy,
      tiling_utils.combine_transforms(
        [tiling_utils.get_reflection_transform(angle, c), 
         (1, 0, 0, 1, dxdy[0], dxdy[1])])) 
      for angle, dxdy in zip(angles, dxdys)]
