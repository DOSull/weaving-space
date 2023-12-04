import pytest
import numpy as np

from weavingspace import weave_matrices

def test_reps_needed():
  assert weave_matrices.reps_needed(3, 4) == (4, 3)

def test_encode_biaxial_weave():
  expected = np.array([4, 5, 2,
                       5, 4, 2,
                       1, 1, 3]).reshape(3, 3)
  pattern = np.array( [0, 1, 0,
                       1, 0, 1,
                       0, 1, 0]).reshape(3, 3)
  warp = np.array([ 1,  1,  1,
                    2,  2,  2,
                   -1, -1, -1]).reshape(3, 3)
  weft = np.array([ 3,  4, -1,
                    3,  4, -1,
                    3,  4, -1]).reshape(3, 3)
  np.testing.assert_equal(
    weave_matrices._encode_biaxial_weave(pattern, warp, weft), expected)

def test_get_pattern():
  expected = np.array([0, 1, 0, 1,
                       1, 0, 1, 0,
                       0, 1, 0, 1,
                       1, 0, 1, 0]).reshape(4, 4)
  tu = np.array([0, 1,
                 1, 0]).reshape(2, 2)
  tr = np.array([1, 0,
                 0, 1]).reshape(2, 2)
  th = np.array([1, 0,
                 0, 1]).reshape(2, 2)
  np.testing.assert_equal(
    weave_matrices.get_pattern(tu, tr, th, 2, 2, 2), expected)
  np.testing.assert_equal(
    weave_matrices.get_pattern(tu, tr, th, 1, 2, 2), expected)
  np.testing.assert_equal(
    weave_matrices.get_pattern(tu, tr, th, 2, 1, 2), expected)
  np.testing.assert_equal(
    weave_matrices.get_pattern(tu, tr, th, 1, 1, 1), expected[0:2, 0:2])
  np.testing.assert_equal(
    weave_matrices.get_pattern(np.fliplr(tu), tr, th, 2, 2, 2),
    np.fliplr(expected))
  np.testing.assert_equal(
    weave_matrices.get_pattern(tu, np.fliplr(tr), th, 2, 2, 2),
    np.fliplr(expected))
  np.testing.assert_equal(
    weave_matrices.get_pattern(tu, tr, np.fliplr(th), 2, 2, 2),
    np.fliplr(expected))

def test_make_plain_pattern():
  expected = np.array([0, 1,
                       1, 0]).reshape(2, 2)
  np.testing.assert_equal(weave_matrices.make_plain_pattern(2), expected)

def test_make_twill_pattern():
  expected = np.array([0, 1, 1, 0,
                       0, 0, 1, 1,
                       1, 0, 0, 1,
                       1, 1, 0, 0]).reshape(4, 4)
  np.testing.assert_equal(
    weave_matrices.make_twill_pattern(2, 2, 2), expected)

def test_make_over_under_row():
  expected = [1, 0, 0, 1, 1, 1,
              0, 1, 1, 0, 0, 0]
  assert weave_matrices.make_over_under_row((1, 2, 3)) == expected
  assert weave_matrices.make_over_under_row(1) == [1, 0]

def test_wrap_row():
  input = [1, 2, 3, 4, 5]
  assert weave_matrices.wrap_row(input, 1) == [5, 1, 2, 3, 4]
  assert weave_matrices.wrap_row(input, -1) == [2, 3, 4, 5, 1]

def test_make_twill_matrix():
  expected = np.array([0, 1, 1, 0,
                       0, 0, 1, 1,
                       1, 0, 0, 1,
                       1, 1, 0, 0]).reshape(4, 4)
  np.testing.assert_equal(weave_matrices.make_twill_matrix(2), expected)

def test_make_basket_pattern():
  expected = np.array([1, 1, 0, 0,
                       1, 1, 0, 0,
                       0, 0, 1, 1,
                       0, 0, 1, 1]).reshape(4, 4)
  np.testing.assert_equal(weave_matrices.make_basket_pattern(2, 2, 2), expected)

def test_make_basket_matrix():
  expected = np.array([1, 1, 0, 0,
                       1, 1, 0, 0,
                       0, 0, 1, 1,
                       0, 0, 1, 1]).reshape(4, 4)
  np.testing.assert_equal(weave_matrices.make_basket_matrix(2), expected)

