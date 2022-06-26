#!/usr/bin/env python
# coding: utf-8

"""Functions to generate the matrices summarising elements of tileable 
repeating elements that when repeated across a map area give the appearance of 
a  biaxial woven surface composed of criss-crossing rectangular elements. 

Implementation is based on ideas discussed in variously
 
+ Glassner, A. 2002. Digital weaving. 1. IEEE Computer Graphics and Applications 22 (6):108–118.
+ ———. 2003a. Digital weaving. 3. IEEE Computer Graphics and Applications 23 (2):80–83.
+ ———. 2003b. Digital weaving. 2. IEEE Computer Graphics and Applications 23 (1):77–90.
 
and (unpublished)

+ Griswold, R. 2006. Mathematical and Computational Topics in Weaving
https://www2.cs.arizona.edu/patterns/weaving/webdocs/mo/Griswold-MO.pdf 
(accessed 29/10/21).
 
where weaving is shown to be a matrix multiplication of tie-up, threading and 
treadling matrices. An accessible introduction can be found at 
https://www.youtube.com/watch?v=oMOSiag3dxg
"""

from typing import Union
import numpy as np

def reps_needed(x1:int, x2:int) -> tuple[int]:
    """Returns how many repetitions of sequences of length x1 and x2 are 
    needed to make a matched length sequence 

    Args:
        x1 (int): length of first sequence.
        x2 (int): length of second sequence.

    Returns:
        tuple[int]: (number of repeats of x1, number of repeats of x2).
    """  
    n = np.lcm(x1, x2)
    return tuple(i for i in (n // x1, n // x2))


def get_pattern(tie_up:np.ndarray, treadling:np.ndarray, threading:np.ndarray,
                warp_n:int, weft_n:int, rep:int = 1) -> np.ndarray:
    """Returns a 0/1 encoded weave matrix. 
    
    Given tie_up, treadling and threading matrices. The following conditions  must be satisfied to avoid non-conformable matrix error:
    
        treadling.shape[1] == tie_up.shape[0]
        tie_up.shape[1] == threading.shape[0]
    
    The "twill", "random", "basket" and "plain" options should guarantee
    this, but the "this" option requires the user to make this happen
    If the warp_n and weft_n values are not factors of treadling.shape[0] and
    threading.shape[1] respectively, the output matrix will be repeated as
    needed to make this match

    Args:
        tie_up (np.ndarray): the tie up matrix.
        treadling (np.ndarray): the treadling matrix.
        threading (np.ndarray): the threading matrix.
        warp_n (int): the number of uniquely labelled warp threads.
        weft_n (int): the number of uniquely labelled weft threads.
        rep (int, optional): optional additional repetition of the output       
            applied in both directions. Defaults to 1.

    Returns:
        np.ndarray: _description_
    """    
    pattern = (treadling @ tie_up @ threading > 0) + 0
    rep_warp = reps_needed(warp_n, pattern.shape[1])
    rep_weft = reps_needed(weft_n, pattern.shape[0])
    return np.tile(pattern, (rep_weft[1] * rep, rep_warp[1] * rep))


# Note that as currently written this function requires the warp and weft
# matrices to be the same size, which get_weave_pattern_matrix will ensure, b
# but which may not be the case if called from elsewhere
def _encode_biaxial_weave(pattern:np.ndarray, warp:np.ndarray, 
                         weft:np.ndarray) -> np.ndarray:
    """Encodes a biaxial weave pattern as below.

        1 - warp is absent
        2 - weft is absent
        3 - both threads are absent
        4 - weft is on top
        5 - warp is on top 
    
    Args:
        pattern (np.ndarray): pattern matrix with 1 where warp is on top, 
            0 where weft is on top.
        warp (np.ndarray): warp matrix where -1 values indicate absent threads.
        weft (np.ndarray): weft matrix where -1 values indicate absent threads.

    Returns:
        np.ndarray: matrix of values encoded as above.
    """    
    # Can't recall why this particular encoding, but don't think it matters...
    pattern = np.where(pattern == 1, 5, pattern)           # warp on top
    pattern = np.where(pattern == 0, 4, pattern)           # weft on top
    pattern = np.where(warp < 0, 1, pattern)               # warp absent
    pattern = np.where(weft < 0, 2, pattern)               # weft absent
    return np.where((warp < 0) & (weft < 0), 3, pattern)   # both absent


def make_plain_pattern(warp_n:int = 1, weft_n:int = 1) -> np.ndarray:
    """Returns plain weave (checkerboard) matrix.

    Args:
        warp_n (int, optional): number of warp thread labels. Defaults to 1.
        weft_n (int, optional): number of weft thread labels. Defaults to 1.

    Returns:
        np.ndarray: 0/1 matrix where 1 = warp on top.
    """    
    return make_twill_pattern(n = 1, warp_n = warp_n, weft_n = weft_n)


def make_twill_pattern(n:Union[int, tuple[int]] = 2, 
                       warp_n:int = 2, weft_n:int = 2) -> np.ndarray:
    """Returns twill pattern matrix extended for warp and weft patterns.

    n is the number of over-unders. With n = 1 we get a plain weave.

    Args:
        n (Union[int,tuple[int]], optional): specifies over-under sequence in
            the weave. Defaults to 2.
        warp_n (int, optional): number of warp thread labels. Defaults to 2.
        weft_n (int, optional): number of weft thread labels. Defaults to 2.

    Returns:
        np.ndarray: _description_
    """    
    tie_up = make_twill_matrix(n)
    threading = np.diag(np.ones(tie_up.shape[0]))
    treadling = np.diag(np.ones(tie_up.shape[1]))
    return get_pattern(tie_up, treadling, threading, warp_n, weft_n)


def make_over_under_row(n:Union[int,tuple[int]]) -> list[int]:
    """Returns a tuple of runs of 1s and 0s.

    Returns a tuple of runs of 1s and 0s per supplied n. 
    
    If n is an integer, returned tuple will be a series of n 1s followed by
    n 0s.
    
    If n is a tuple of odd length it will be repeated to produce an even
    length over-under sequence. This is required to avoid cases such as
    e.g, (1,2,3) -> 100111 which repeated is 1001111001111, i.e. a (2,4)
    pattern. Repeating it yields 100111011000 which is the requested pattern.

    Args:
        n (Union[int,tuple[int]]): requested over-under sequence. See details.

    Returns:
        list[int]: list of runs of 1s and 0s in the requested pattern. 
    """    
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


def wrap_row(r:list, by:int = 1) -> list:
    """Wraps the list r by number of positions.
    
    Positive by will shift right.
    Negative shifts left.

    Args:
        r (list): list to wrap.
        by (int, optional): number of positions to shift by. Defaults to 1.

    Returns:
        list: wrapped list.
    """    
    return r[-by:] + r[:-by]


def make_twill_matrix(over_under:Union[int, tuple[int]]) -> np.ndarray:
    """Makes a twill 0/1 matrix.

    Makes a matrix like
    
        1 1 0 0
        0 1 1 0
        0 0 1 1
        1 0 0 1
      
    where the repeat runs in each row are lengths returned by
    make_over_under_rown(n)

    Args:
        over_under (Union[int,tuple[int]]): over-under run specification. See
            make_over_under_row().

    Returns:
        np.ndarray: a matrix of 0s and 1s.
    """    
    row = make_over_under_row(over_under)
    d = len(row)
    out = []
    for i in range(d):
        row = wrap_row(row, 1)
        out.extend(row)
    return np.array(out).reshape(d, d)


def make_basket_pattern(n:int = 2, warp_n:int = 2, 
                        weft_n:int = 2) -> np.ndarray: 
    """Returns basket pattern matrix extended for warp and weft patterns.

    Args:
        n (int, optional): over under count. Defaults to 2.
        warp_n (int, optional): number of warp thread labels. Defaults to 2.
        weft_n (int, optional): number of weft thread labels. Defaults to 2.

    Returns:
        np.ndarray: _description_
    """    
    tie_up = make_basket_matrix(n)
    threading = np.diag(np.ones(tie_up.shape[0]))
    treadling = np.diag(np.ones(tie_up.shape[1]))
    return get_pattern(tie_up, treadling, threading, warp_n, weft_n)


def make_basket_matrix(n:int) -> np.ndarray:
    """Returns a basket weave pattern matrix.

    Makes a matrix like
    
        1 1 0 0
        1 1 0 0
        0 0 1 1
        0 0 1 1
    
    where the repeat runs in each row are length n

    Args:
        n (int): dimension of the 'basket' over-under pattern.

    Returns:
        np.ndarray: 0/1 matrix in basket weave pattern.
    """    
    return np.array((([1] * n + [0] * n) * n) + 
                    (([0] * n + [1] * n) * n)).reshape(n * 2, n * 2)


def make_this_pattern(tie_up:np.ndarray, 
                      threading:np.ndarray = None,
                      treadling:np.ndarray = None,
                      warp_n:int = 1, weft_n:int = 1) -> np.ndarray:
    """Pass through returns weave pattern matrix from supplied input.

    This is just a pass through function which applies suitable size identity. Could try to enforce
    
        treadling.shape[1] == tie_up.shape[0] and 
        tie_up.shape[1] == threading.shape[0]
    
    but unsure what would be an appropriate way to do this...

    Args:
        tie_up (np.ndarray): desired weave pattern.
        threading (np.ndarray): desired threading pattern.
        treadling (np.ndarray): desired treadling pattern.
        warp_n (int, optional): number of warp thread labels. Defaults to 1.
        weft_n (int, optional): number of weft thread labels. Defaults to 1.

    Returns:
        np.ndarray: resulting 0/1 weave pattern matrix.
    """    
    threading = (np.diag(np.ones(tie_up.shape[0]))
                 if threading is None
                 else threading)
    treadling = (np.diag(np.ones(tie_up.shape[1]))
                 if treadling is None
                 else treadling)
    return get_pattern(tie_up, treadling, threading, warp_n, weft_n)


def get_weave_pattern_matrix(weave_type:str = "plain", 
                             n:Union[int, tuple[int]] = 2, 
                             warp:Union[list[str], tuple[str]] = ["a", "b"],
                             weft:Union[list[str], tuple[str]] = ["c", "d"], 
                             tie_up:np.ndarray = make_twill_matrix((2, 2)), 
                             tr:np.ndarray = None, 
                             th:np.ndarray = None) -> np.ndarray:
    """Returns encoded weave pattern matrix.
    
    See `_encode_biaxial_weave()` for the encoding.

    Allowed weave_types: "plain", "twill", "basket", and "this" (pass-thru).
    These are explained in the respective functions to which this function
    delegates construction of the base matrices before applying the encoding.

    Args:
        weave_type (str, optional): one of "plain", "twill", "basket" or
            "this". Defaults to "plain".
        n (Union[int,tuple[int]], optional): over under pattern. See 
            make_over_under_row() for details. Defaults to 2.
        warp (Union[list[str],tuple[str]], optional): list of labels for warp 
            strands. Defaults to ["a", "b"].
        weft (Union[list[str],tuple[str]], optional): list of labels for weft 
            strands. Defaults to ["c", "d"].
        tie_up (np.ndarray, optional): a weave pattern matrix to pass thru in 
            the "this" case. Defaults to make_twill_matrix((2, 2)).
        tr (np.ndarray, optional): treadling matrix for the "this" case.        
            Defaults to None.
        th (np.ndarray, optional): threading matrix for the "this" case.        
            Defaults to None.

    Returns:
        np.ndarray: encoded weave pattern matrix.
    """    
    warps = [-1 if c == "-" else i for i, c in enumerate(warp)]
    wefts = [-1 if c == "-" else i for i, c in enumerate(weft)]
    width = len(warp)
    height = len(weft)
    
    if weave_type == "plain":
        p = make_plain_pattern(warp_n = width, weft_n = height)
    elif weave_type == "twill":
        p = make_twill_pattern(n = n, warp_n = width, weft_n = height)
    elif weave_type == "basket":
        p = make_basket_pattern(n = n, warp_n = width, weft_n = height)
    else:
        p = make_this_pattern(tie_up = tie_up, threading = th, treadling = tr,
                              warp_n = width, weft_n = height)
    nr, nc = p.shape
    warp_threads = np.array(warps * 
            reps_needed(nc, len(warps))[1] * nr).reshape((nr, nc))
    weft_threads = np.array(wefts * 
            reps_needed(nr, len(wefts))[1] * nc).reshape((nc, nr)).transpose()
    # encode to reflect missing threads
    return _encode_biaxial_weave(p, warp_threads, weft_threads)

