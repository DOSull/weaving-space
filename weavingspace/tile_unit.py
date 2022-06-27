#!/usr/bin/env python
# coding: utf-8

"""The `TileUnit` subclass of `weavingspace.tileable.Tileable` implements
many 'conventional' tilings of the plane.

Examples:
    A `TileUnit` is initialised like this
    
        tile_unit = TileUnit(tiling_type = "cairo")
        
    The `tiling_type` may be one of the following
    
    + "cairo" the Cairo tiling more formally known as the Laves 
    [3<sup>2</sup>.4.3.4] tiling. The author's favourite tiling, hence it 
    has its own tiling_type.
    + "hex-dissection" a range of dissections of the regular hexagon into,
    2, 3, 4, 6, or 12 'pie slices'. The number of slices is set by 
    specifying an additional argument `n`. Slices are cut either starting 
    at the corners of  the hexagon or from the midpoints of hexagon edges, 
    by specifying an additional argument `dissection_offset` set to either 
    0 or 1 respectively.
    + "laves" a range of isohedral tilings. See [this article](https://en.wikipedia.org/wiki/List_of_Euclidean_uniform_tilings#Laves_tilings). 
    The desired tiling is specified by the additional argument `code` which 
    is a string like "3.3.4.3.4".
    + "archimedean" a range of tilings by regular polygons. See [this 
    article](https://en.wikipedia.org/wiki/Euclidean_tilings_by_convex_regular_polygons#Archimedean,_uniform_or_semiregular_tilings). Many of these are the dual tilings of 
    the Laves tilings. The desired tiling is specified by the additional 
    argument `code` which is a string like "3.3.4.3.4". Not all the 
    possible Archimedean tilings are implemented.
    + "hex-colouring" three colourings of the regular hexagon tiling, of 
    either 3, 4, or 7 colours, as specified by the argument `n`.
    + "square-colouring" one colouring of the regular square tiling, of 5 
    colours as specified by the argument `n = 5`.
    
    See [this notebook](https://github.com/DOSull/weaving-space/blob/main/weavingspace/all-the-tiles.ipynb) for exact usage, and illustrations of 
    each tiling. 
    
    Spacing and coordinate reference of the tile unit are specified by the
    `weavingspace.tileable.Tileable` superclass variables 
    `weavingspace.tileable.Tileable.spacing` and 
    `weavingspace.tileable.Tileable.crs`.

    Base tilings by squares, hexagons or triangles can also be requested 
    using
    
        tile_unit = TileUnit()  # square tiling, the default
        tile_unit = TileUnit(tile_shape = TileShape.HEXAGON)
        tile_unit = TileUnit(tile_shape = TileShape.TRIANGLE)
        
    The first two of these have only one element_id value, and so cannot be 
    used for multivariate mapping. The triangle case has two element_id 
    values so may be useful in its base form.
    
    To create custom tilings start from one of the base tiles above, and 
    explicitly set the `weavingspace.tileable.Tileable.elements` variable 
    by geometric construction of suitable shapely.geometry.Polygons. TODO: A detailed example of this usage can be found here ....
"""

import copy
from dataclasses import dataclass
import string

import geopandas as gpd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine

from weavingspace.tileable import Tileable
from weavingspace.tileable import TileShape

import weavingspace.tiling_utils as tiling_utils
import weavingspace.tiling_geometries as tiling_geometries


@dataclass
class TileUnit(Tileable):
    """Class to represent the tileable elements of a 'conventional' tiling.
    
    Args:
        tiling_type (str): tiling type as detailed above.
        dissection_offset (int): offset for "hex-dissection" tilings. See above 
            for details. Defaults to 1. 
        n (int): number of dissections or colours in "hex-dissection", 
            "hex-colouring", or "square-colouring" tiling types. Defaults to 3.
        code (str): the code for "laves" or "archimedean" tiling types. 
        Defaults to "3.3.4.3.4".

    Returns:
        _type_: _description_
    """
    tiling_type:str = None
    dissection_offset:int = 1
    n:int = 3
    code:str = "3.3.4.3.4"
        
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        if not self.tiling_type is None:
            self.tiling_type = self.tiling_type.lower()
        if self.tile_shape == TileShape.TRIANGLE:
            self._modify_tile()
            self._modify_elements()
            self.setup_vectors()
            self.setup_regularised_tile_from_elements()
        if self.regularised_tile is None:
            self.setup_regularised_tile_from_elements()


    def _setup_tile_and_elements(self) -> None:
        """Delegates setup of the unit to various functions depending
        on self.tiling_type.
        """
        if self.tiling_type == "cairo":
            tiling_geometries.setup_cairo(self)
        elif self.tiling_type == "hex-dissection":
            tiling_geometries.setup_hex_dissection(self)
        elif self.tiling_type == "laves":
            tiling_geometries.setup_laves(self)
        elif self.tiling_type == "archimedean":
            tiling_geometries.setup_archimedean(self)
        elif self.tiling_type in ("hex-colouring", "hex-coloring"):
            tiling_geometries.setup_hex_colouring(self)
        elif self.tiling_type in ("square-colouring", "square-coloring"):
            tiling_geometries.setup_square_colouring(self)
        else:
            tiling_geometries._setup_none_tile(self)
        return

    
    def _modify_elements(self) -> None:
        """It is not trivial to tile a triangle, so this function augments
        augments the elements of a triangular tile to a diamond by 180 degree 
        rotation. Operation is 'in place'.
        """
        elements = self.elements.geometry
        ids = list(self.elements.element_id)
        
        new_ids = list(string.ascii_letters[:(len(ids) * 2)])
        elements = elements.translate(0, -elements.total_bounds[1])
        twins = [affine.rotate(element, a, origin = (0, 0)) 
                 for element in elements
                 for a in range(0, 360, 180)]
        self.elements = gpd.GeoDataFrame(
            data = {"element_id": new_ids},
            geometry = gpd.GeoSeries(twins), crs = self.elements.crs)
        return None


    def _modify_tile(self) -> None:
        """It is not trivial to tile a triangular tile so this function 
        changes the tile to a diamond by manually altering the tile in place
        to be a diamond shape.
        """
        tile = self.tile.geometry[0]
        # translate to sit on x-axis
        tile = affine.translate(tile, 0, -tile.bounds[1])
        pts = [p for p in tile.exterior.coords]
        pts[-1] = (pts[1][0], -pts[1][1])
        self.tile.geometry = gpd.GeoSeries([geom.Polygon(pts)], crs = self.crs)
        self.tile_shape = TileShape.DIAMOND
        return None

    
    def _get_legend_key_shapes(self, polygon:geom.Polygon, 
                               counts:int = 25, angle:float = 0, 
                               radial:bool = False) -> list[geom.Polygon]:
        """Returns a set of shapes that can be used to make a legend key 
        symbol for the supplied polygon. In TileUnit this is a set of 'nested
        polygons.

        Args:
            polygon (geom.Polygon): the polygon to symbolise.
            n (int, optional): number of steps. Defaults to 25.
            rot (float, optional): rotation that may have to be applied.  
                Not used in the TileUnit case. Defaults to 0.

        Returns:
            list[geom.Polygon]: a list of nested polygons.
        """
        if not radial:
            n = sum(counts)
            # bandwidths = list(np.cumsum(counts))
            bandwidths = [c / n for c in counts]
            bandwidths = [bw if bw > 0.05 or bw == 0 else 0.05 
                          for bw in bandwidths]
            n = sum(bandwidths)
            bandwidths = [0] + [bw / n for bw in bandwidths]
            # # make buffer widths that will yield approx equal area 'annuli'
            # bandwidths = range(n_steps + 1)
            # sqrt exaggerates outermost annuli, which can otherwise disappear
            bandwidths = [np.sqrt(bw) for bw in bandwidths]
            distances = np.cumsum(bandwidths)
            # get the negative buffer distance that will 'collapse' the polygon
            radius = tiling_utils.get_collapse_distance(polygon)
            distances = distances * radius / distances[-1]
            nested_polys = [polygon.buffer(-d, resolution = 1, 
                                           join_style = 2) for d in distances]
            # return converted to annuli (who knows someone might set alpha < 1)
            nested_polys = [g1.difference(g2) for g1, g2 in 
                            zip(nested_polys[:-1], nested_polys[1:])]
            return [p for c, p in zip(counts, nested_polys) if c > 0]      
        else:
            n = sum(counts)
            slice_posns = list(np.cumsum(counts))
            slice_posns = [0] + [p / n for p in slice_posns]
            return [tiling_utils.get_polygon_sector(polygon, i, j) 
                    for i, j in zip(slice_posns[:-1], slice_posns[1:])]


    # Note that geopandas clip is not order preserving hence we do this
    # one polygon at a time...
    def inset_tile(self, d:float = 0) -> "TileUnit":
        """Returns a new TileUnit clipped by `self.regularised_tile` after
        a negative buffer d has been applied.

        Args:
            d (float, optional): the inset distance. Defaults to 0.

        Returns:
            TileUnit: the new TileUnit with inset applied.
        """
        inset_tile = self.regularised_tile.geometry.buffer(
                -d, resolution = 1, join_style = 2)[0]
        # the clean_geometry seems needed to stop proliferation of vertices
        new_elements = [tiling_utils.clean_polygon(inset_tile.intersection(e))
                        for e in self.elements.geometry]
        result = copy.deepcopy(self)
        result.elements.geometry = gpd.GeoSeries(new_elements)
        return result
    
    
    def scale_elements(self, sf:float = 1) -> "TileUnit":
        """Scales the elements by the specified factor, centred on (0, 0).

        Args:
            sf (float, optional): scale factor to apply. Defaults to 1.

        Returns:
            TileUnit: the scaled TileUnit.
        """
        result = copy.deepcopy(self)
        result.elements.geometry = self.elements.geometry.scale(
            sf, sf, origin = (0, 0))
        return result


