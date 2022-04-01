#!/usr/bin/env python
# coding: utf-8

# # Tileable biaxial weave units
# Functions that can be used to generate sf data 'weave units' i.e. a tileable 
# repeating element that when tiled gives the appearance of a  biaxial woven 
# surface composed of criss-crossing rectangular elements. Implementation is 
# based on ideas discussed in variously
# 
# Glassner, A. 2002. Digital weaving. 1. IEEE Computer Graphics and Applications 
#   22 (6):108–118.
# ———. 2003a. Digital weaving. 3. IEEE Computer Graphics and Applications 
#   23 (2):80–83.
# ———. 2003b. Digital weaving. 2. IEEE Computer Graphics and Applications 
#   23 (1):77–90.
# 
# and (unpublished)
# 
# Griswold, R. 2006. Mathematical and Computational Topics in Weaving
# https://www2.cs.arizona.edu/patterns/weaving/webdocs/mo/Griswold-MO.pdf 
# (accessed 29/10/21).
# 
# where weaving is shown to be a matrix multiplication of tie-up, threading and 
# treadling matrices. An accessible introduction can be found at 
# https://www.youtube.com/watch?v=oMOSiag3dxg

import numpy as np
from loom import Loom
from render_weave_grids import make_shapes_from_coded_weave_matrix
from weaving_space_utils import get_strand_ids

# # augment matrix by adding `by` rows and columns identical to 
# # the first row and first column of 
# def augment_matrix(m, by = 1):
#     s = np.array(m.shape)
#     s_aug = s + by
#     scale = tuple(np.ceil(s_aug / s).astype(np.int32))
#     return np.tile(m, scale)[0:(tuple(s_aug)[0]), 0:(tuple(s_aug)[1])]


# def augment_matrix_with_values(m, by = 1, values = 0):
#     output = augment_matrix(m, by)
#     s1 = m.shape
#     s2 = output.shape
#     output[s1[0]:s2[0], :] = values
#     output[:, s1[1]:s2[1]] = values
#     return output


def reps_needed(x1, x2):
    n = np.lcm(x1, x2)
    return tuple(i for i in (n // x1, n // x2))


# Returns a 1/2 encoded weave matrix given tie_up, treadling and
# threading matrices. The following conditions must be satisfied to
# avoid non-conformable matrix error:
# The "twill", "random", "basket" and "plain" options should guarantee
# this, but the "this" option requires the user to make this happen
# If the warp_n and weft_n values are not factors of nrow(treadling) and
# ncol(threading) respectively, the output matrix will be repeated as
# needed to make this match
def get_pattern(tie_up, treadling, threading, warp_n, weft_n, rep = 1):
    pattern = (treadling @ tie_up @ threading > 0) + 0
    rep_warp = reps_needed(warp_n, pattern.shape[1])
    rep_weft = reps_needed(weft_n, pattern.shape[0])
    return np.tile(pattern, (rep_weft[1] * rep, rep_warp[1] * rep))


# Note that as currently written this function requires the warp and weft
# matrices to be the same size, which get_weave_pattern_matrix will ensure, b
# but which may not be the case if called from elsewhere
def encode_biaxial_weave(pattern, warp, weft):
    pattern = np.where(pattern == 1, 5, pattern)         # warp present and on top
    pattern = np.where(pattern == 0, 4, pattern)         # weft present and on top
    pattern = np.where(warp < 0, 1, pattern)             # warp absent
    pattern = np.where(weft < 0, 2, pattern)             # weft absent
    return np.where((warp < 0) & (weft < 0), 3, pattern)   # both absent


# THE WEAVES
# -- Plain --
# simple over-under weave
def make_plain_pattern(warp_n = 1, weft_n = 1):
    return make_twill_pattern(n = 1, warp_n = warp_n, weft_n = weft_n)


# -- Twills --
# twill weave with n the number of over-unders
# note this is used with n = 1 to make plain weaves
def make_twill_pattern(n = 2, warp_n = 2, weft_n = 2):
    over_under = n
    if type(over_under) == int:
        over_under = (over_under, over_under)
    tie_up = make_twill_matrix(over_under)
    threading = np.diag(np.ones(tie_up.shape[0]))
    treadling = np.diag(np.ones(tie_up.shape[1]))
    return get_pattern(tie_up, treadling, threading, warp_n, weft_n)


# returns a vector of runs of 1s and 0s per the supplied vector.
# If the length of n is odd then it is doubled to produce an even
# length over-under sequence that repeats. If we don't do this then,
# e.g, 1:3 becomes 100111 which repeated is 1001111001111, i.e. a 2-4
# over-under pattern. Doubling it makes 100111011000 which has the
# requested pattern
def make_over_under_row(n):
    over_under = n
    if type(over_under) == int:
        over_under = [n, n]
    elif len(n) % 2 != 0:
        over_under = n * 2
    x = 1
    row = []
    for y in over_under:
        row.extend([x] * y)
        x = 1 - x
    return row


# wraps a vector
# by : the number of positions to shift the row
# r  : the row
def wrap_row(by, r):
    return r[-by:] + r[:-by]


# makes a matrix like
# 1 1 0 0
# 0 1 1 0
# 0 0 1 1
# 1 0 0 1
# where the repeat runs in each row are length n
def make_twill_matrix(over_under):
    row = make_over_under_row(over_under)
    d = len(row)
    out = []
    for i in range(d):
        row = wrap_row(1, row)
        out.extend(row)
    return np.array(out).reshape(d, d)


# -- Basket weave --
def make_basket_pattern(n = 2, warp_n = 2, weft_n = 2):
    tie_up = make_basket_matrix(n)
    threading = np.diag(np.ones(tie_up.shape[0]))
    treadling = np.diag(np.ones(tie_up.shape[1]))
    return get_pattern(tie_up, treadling, threading, warp_n, weft_n)


# makes a matrix like
# 1 1 0 0
# 1 1 0 0
# 0 0 1 1
# 0 0 1 1
# where the repeat runs in each row are length n
def make_basket_matrix(n):
    return np.array((([1] * n + [0] * n) * n) + 
                    (([0] * n + [1] * n) * n)).reshape(n * 2, n * 2)


# -- Other stuff --
# This is just a pass through function. Could try to enforce
#   ncol(treadling) == nrow(tie_up) and ncol(tie_up) == nrow(threading)
# but unsure what would be an appropriate way to do this...
def make_this_pattern(tie_up, threading, treadling, warp_n = 2, weft_n = 2):
    threading = np.diag(np.ones(tie_up.shape[0]))
    treadling = np.diag(np.ones(tie_up.shape[1]))
    return get_pattern(tie_up, treadling, threading, warp_n, weft_n)


# returns a matrix giving where 1 indicates warp on top, 2
# indicates weft on top
# types:
# "plain"    : over 1 under 1 in both directions
# "twill"    : over n under n each weft thread, shifting
#              one along between rows
# n          : the over-under for twill patterns
# warp, weft : a vector of distinct values (ints or chars) where each
#              indicates a different thread colour; repeats are allowed,
#              and "-" indicates that a thread should be skipped
def get_weave_pattern_matrix(*,
        weave_type = "plain", n = 2, warp = list("ab"), weft = list("cd"),
        tie_up = make_twill_matrix((2, 2)), tr = np.diag(np.ones(2)), th = np.diag(np.ones(2))):
    
    warps = [-1 if c == "-" else i for i, c in enumerate(warp)]
    wefts = [-1 if c == "-" else i for i, c in enumerate(weft)]
    width = len(warp)
    height = len(weft)
    
    if weave_type == "plain":
        p = make_plain_pattern()
    elif weave_type == "twill":
        p = make_twill_pattern(n = n, warp_n = width, weft_n = height)
    elif weave_type == "basket":
        p = make_basket_pattern(n = n, warp_n = width, weft_n = height)
    else:
        p = make_this_pattern(tie_up = tie_up, threading = th, treadling = tr,
                               warp_n = width, weft_n = height)
    nr, nc = p.shape
    warp_threads = np.array(warps * reps_needed(nc, len(warps))[1] * nr).reshape((nr, nc))
    weft_threads = np.array(wefts * reps_needed(nr, len(wefts))[1] * nc).reshape((nc, nr)).transpose()
    
    # warp_threads = np.repeat(warps, np.prod(p.shape) / len(warps)).reshape(p.shape)
    # weft_threads = np.repeat(wefts, np.prod(p.shape) / len(wefts)).reshape(p.shape)
    # # encode to reflect missing threads
    return encode_biaxial_weave(p, warp_threads, weft_threads)


# EXTERNAL API
def get_biaxial_weave_unit(*, spacing = 10_000, aspect = 1, margin = 0,
        weave_type = "twill", n = (2, 2), strands = "ab|cd", crs = 3857,
        tie_up = make_twill_matrix((2, 2)), tr = np.diag(np.ones(2)), th = np.diag(np.ones(2))):
    strand_ids = get_strand_ids(strands)
    warp_threads = strand_ids[0]
    weft_threads = strand_ids[1]
    
    if weave_type == "basket":
        n = n[0]
    
    treadling = np.diag(np.ones(tie_up.shape[0]))
    threading = np.diag(np.ones(tie_up.shape[1]))
    
    p = get_weave_pattern_matrix(weave_type = weave_type, n = n, 
                    warp = warp_threads, weft = weft_threads, 
                    tie_up = tie_up, tr = tr, th = th)

    return make_shapes_from_coded_weave_matrix(
        Loom([p]), spacing = spacing, width = aspect, margin = margin, 
        axis1_threads = weft_threads, axis2_threads = warp_threads, crs = crs)
