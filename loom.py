#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass

import numpy as np


@dataclass
class Loom:
  indices: list[tuple]
  orderings: list
  parity: int
  dimensions: tuple
  orientations: tuple
  n_axes: int
  
  def __init__(self, matrices):
    if len(matrices) == 1:
      m = matrices[0]
      self.dimensions = m.shape
      self.n_axes = len(self.dimensions)
      self.parity = None
      self.orientations = (0, -90)
      self.indices = [(i, j) for i in range(m.shape[0]) \
                            for j in range(m.shape[1])]
      self.orderings = [decode_biaxial_to_order(m[ij]) for ij in self.indices]
    else:
      nA, nB, nC = [max(m.shape) for m in matrices]
      self.dimensions = (nA, nB, nC)
      self.n_axes = len(self.dimensions)
      self.parity = (nA + nB + nC - 3) // 2
      self.orientations = (0, 120, 240)
      all_indices = [(i, j, k) for i in range(nA) for j in range(nB) for k in range(nC)]
      self.indices = [x for x in all_indices if sum(x) in (self.parity, self.parity + 1)]

      m1, m2, m3 = matrices
      mAB = np.tile(m1, (np.lcm(nA, m1.shape[0]) // m1.shape[0],
                         np.lcm(nA, m1.shape[1]) // m1.shape[1]))
      mBC = np.tile(m2, (np.lcm(nB, m2.shape[0]) // m2.shape[0],
                         np.lcm(nB, m2.shape[1]) // m2.shape[1]))
      mCA = np.tile(m3, (np.lcm(nC, m3.shape[0]) // m3.shape[0],
                         np.lcm(nC, m3.shape[1]) // m3.shape[1]))
      
      ordAB = [decode_biaxial_to_order(mAB[ij], axis = 1) for ij in [(x[1], x[0]) for x in self.indices]]
      ordBC = [decode_biaxial_to_order(mBC[ij], axis = 2) for ij in [(x[2], x[1]) for x in self.indices]]
      ordCA = [decode_biaxial_to_order(mCA[ij], axis = 3) for ij in [(x[0], x[2]) for x in self.indices]]
      self.orderings = [combine_orders(abc) for abc in zip(ordAB, ordBC, ordCA)]
    #   self.orderings = [combine_orderings(abc) for abc in zip(ordAB, ordBC, ordCA)]


bi_to_tri = (
    {1: (0, ), 2: (1, ), 3: "NA", 4: (0, 1,), 5: (1, 0,)},
    {1: (1, ), 2: (0, ), 3: "NA", 4: (1, 0,), 5: (0, 1,)},
    {1: (2, ), 2: (1, ), 3: "NA", 4: (2, 1,), 5: (1, 2,)},
    {1: (0, ), 2: (2, ), 3: "NA", 4: (0, 2,), 5: (2, 0,)}
)

def decode_biaxial_to_order(code, axis = 0):
    return bi_to_tri[axis][code]

combine = {
    "NA":   { "NA":  { "NA":  None, 
                          2:  (2, )}, 
                 1:  { "NA":  (1, )},
                 2:  {    2:  (2, )}},
    0:      { "NA":  {    0:  (0, )},
                 2:  {(2,0):  (2, 0),
                      (0,2):  (0, 2)}},
    1:      {    1:  { "NA":  (1, )},
             (1,2):  {    2:  (1, 2)},
             (2,1):  {    2:  (2, 1)}},
    (0,1):  {    1:  {    0:  (0, 1)},
             (1,2):  {(0,2):  (0, 1, 2)},
             (2,1):  {(2,0):  (2, 0, 1),
                      (0,2):  (0, 2, 1)}},
    (1,0):  {    1:  {    0:  (1, 0)},
             (1,2):  {(2,0):  (1, 2, 0),
                      (0,2):  (1, 0, 2)},
             (2,1):  {(2,0):  (2, 1, 0)}}
}

def combine_orders(orders):
    # print(f"orders: {orders}")
    try:
        result = combine[orders[0]][orders[1]][orders[2]]
    except:
        print(f"Unable to determine unique ordering on {orders}")
        return "NA"
    else:
        return result 


# def match(vals, lst, default = 100):
#     return [default if v not in lst else lst.index(v) for v in vals]
 
# def x_over_y(L, x, y):
#     return L.index(x) < L.index(y) 

# # SOME inconsistencies in ordering of tied values make this unworkable? I think?
# def combine_orderings(orderings, values = (0, 1, 2), verbose = False, default = 100):
#     ranks = [match(values, order, default) for order in orderings]
#     scores = [sum(x) for x in zip(*ranks)]
#     hiscore = len(orderings) * default
#     number_present = sum([s < hiscore for s in scores])
#     if number_present == 0:
#         result = None
#     elif len(set(scores)) != len(scores) and number_present != 1:
#         print(f"Unable to determine ordering on {orderings}")
#         result = "NA"
#     else:
#         sorted_scores = sorted(scores)
#         indexes = [sorted_scores.index(x) for x in scores][:number_present]
#         result = [values[i] for i in indexes]
#     print(f"orderings: {orderings} values: {values} result: {result}")
#     if verbose:
#         return result, ranks, scores
#     else:
#         return result

