#!/usr/bin/env python
# coding: utf-8

from itertools import chain

import numpy as np
import geopandas
from shapely.affinity import translate
from shapely.geometry import Polygon
from shapely.geometry import MultiPolygon
from shapely.ops import unary_union

from loom import Loom
from weave_grids import _WeaveGrid


def centre_offset(shape: Polygon, 
                    target:tuple[float] = (0, 0)) -> tuple[float]:
    """Returns vector required to move centroid of polygon to target. 

    Args:
        shape (Polygon): polygon to move.
        target (tuple[float], optional): target to move to. Defaults to (0, 0).

    Returns:
        tuple[float]: tuple of x, y movement required.
    """  
    shape_c = shape.centroid.coords[0]
    return (target[0] - shape_c[0], target[1] - shape_c[1])


# builds the sf associate with a given weave supplied as 'loom' which is a list
# containing the coordinates in an appropriate grid (Cartesian or triangular)
# and the orderings of the strands at each coordinate location
def make_shapes_from_coded_weave_matrix(
        loom:Loom, spacing:float = 1.0, width:float = 1.0, margin:float = 0.1,
        strand_labels:list[list[str]] = [["a"], ["b"], ["c"]], crs:int = 3857
    ) -> dict:
    """Returns weave elements GeoDataFrame and tile GeoDataFrame in a dictionary

    Args:
        loom (Loom): matrix or stack of matrices representing the weave pattern.
        spacing (float, optional): spacing of weave strands. Defaults to 1.0.
        width (float, optional): width of weave strands relative to spacing.
            Defaults to 1.0.
        margin (float, optional): inset margin applied to strands. Defaults 
            to 0.1.
        strand_labels (list[list[str]], optional): list of lists of labels
            for strands in each direction. Defaults to [["a"], ["b"], ["c"]].
        crs (int, optional): usually an EPSG code, but any object interpretable 
            by geopandas as a CRS. Defaults to 3857.

    Returns:
        dict: dictionary with contents {"weave_unit": GeoDataFrame of weave 
            elements, "tile": GeoDataFrame of the tile}
    """  
    grid = _WeaveGrid(loom.n_axes, loom.orientations, spacing)
    # expand the list of strand labels if needed in each direction
    labels = []
    for dim, thread in zip(loom.dimensions, strand_labels):
        labels.append(thread * int(np.ceil(dim // len(thread))))
    weave_polys = []
    bb_polys = []
    strand_ids = []
    for coords, strand_order in zip(loom.indices, loom.orderings):
        ids = [thread[coord] for coord, thread in zip(coords, labels)]
        cell = grid.get_grid_cell_at(coords)
        bb_polys.append(grid._gridify(cell))
        # print(f"strand_order: {strand_order} ids: {ids}")
        if strand_order is None: continue  # No strands present
        if strand_order == "NA": continue  # Inconsistency in layer order
        n_slices = [len(id) for id in ids]
        next_polys = grid.get_visible_cell_strands(width, coords, 
                                                   strand_order, n_slices)
        weave_polys.extend(next_polys)
        next_labels = [list(ids[i]) for i in strand_order]  # a list of lists
        next_labels = list(chain(*next_labels))  # flatten list of lists
        # print(f"n: {len(next_polys)} labels: {labels}")
        strand_ids.extend(next_labels)
    tile = unary_union(bb_polys)
    shift = centre_offset(tile)
    tile = translate(tile, shift[0], shift[1])
    return {
        "weave_unit": make_weave_gdf(weave_polys, strand_ids, tile, 
                                     shift, spacing, margin, crs),
        "tile": geopandas.GeoDataFrame(geometry = geopandas.GeoSeries([tile]),
                                       crs = crs)
    }


def make_weave_gdf(polys:list[Polygon|MultiPolygon], strand_ids:list[str], 
                   bb:Polygon, offset:tuple[float], spacing:float, margin:float, crs:int) -> geopandas.GeoDataFrame:
    """Makes a GeoDataFrame from weave element polygons, labels, tile, etc.

    Args:
        polys (list[Polygon | MultiPolygon]): list of weave element polygons.
        strand_ids (list[str]): list of strand labels.
        bb (Polygon): bounding shape tile.
        offset (tuple[float]): offset to centre the weave elements on the tile.
        spacing (float): spacing of weave strands.
        margin (float): inset margin to apply to strands.
        crs (int): CRS usually an EPSG code integer, but any geopandas
            acceptable CRS object will work.

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame clipped to the tile, with margin applied.
    """    
    weave = geopandas.GeoDataFrame(
        data = {"strand": strand_ids},
        geometry = geopandas.GeoSeries([translate(p, offset[0], offset[1]) for p in polys])
    )
    weave = weave[weave.strand != "-"]
    weave = weave.dissolve(by = "strand", as_index = False)
    # this buffer operation cleans up some geometry issues
    weave.geometry = weave.buffer(-0.0001 * spacing)
    weave.geometry = weave.buffer((0.0001 - margin) * spacing)
    return weave.clip(bb).set_crs(crs)


