import pytest
import numpy as np

from weavingspace import _weave_grid
from weavingspace import tiling_utils

grid2 = _weave_grid.WeaveGrid(
  n_axes = 2, orientations = (0, -90), spacing = 1)
grid3 = _weave_grid.WeaveGrid(
  n_axes = 3, orientations = (0, 120, 240), spacing = 1)

def test_setup_basis():
  expected2 = np.array([0, 1, 1, 0]).reshape(2, 2)
  expected3 = np.array([0, -np.sqrt(3)/2, np.sqrt(3)/2,
                        1,          -0.5,         -0.5]).reshape(2, 3) * 2 / 3
  np.testing.assert_almost_equal(grid2.basis, expected2)
  np.testing.assert_almost_equal(grid3.basis, expected3)


def test_get_coordinates():
  np.testing.assert_almost_equal(grid2.get_coordinates((1, 2)), (2, 1))
  np.testing.assert_almost_equal(grid3.get_coordinates((1, 2, 3)),
                                 (np.sqrt(3)/3, -1))


def test_get_grid_cell_at():
  assert grid2.grid_cell.equals_exact(grid2.get_grid_cell_at((0, 0)), 1e-8)
  assert len(grid2.get_grid_cell_at((1, 2)).exterior.coords) == 5
  assert grid3.grid_cell.equals_exact(grid3.get_grid_cell_at((0, 0, 0)), 1e-8)
  assert len(grid3.get_grid_cell_at((1, 2, -3)).exterior.coords) == 4


def test_get_visible_cell_strands():
  assert len([p for p in grid2.get_visible_cell_strands()
              if not p.is_empty]) == 1
  assert len([p for p in grid2.get_visible_cell_strands(width = 0.8)
              if not p.is_empty]) == 2
  assert len([p for p in grid3.get_visible_cell_strands()
              if not p.is_empty]) == 1
  assert len([p for p in grid3.get_visible_cell_strands(width = 0.8)
              if not p.is_empty]) == 3


def test_get_tile_from_cells():
  rect = tiling_utils.get_regular_polygon(1.1, 4)
  hex = tiling_utils.get_regular_polygon(1.1, 6)
  assert grid2.get_tile_from_cells(rect).bounds == (-.5, -.5, .5, .5)
  assert all(np.isclose(x, y, 3)
             for x, y in zip(grid3.get_tile_from_cells(rect).bounds,
                             (-.577, -.5, .5, .577)))
