#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
import string

import geopandas as gpd
import numpy as np
import shapely.geometry as geom
import shapely.affinity as affine

from tileable import Tileable
from tileable import TileShape

import tiling_utils
import _tiling_geometries


@dataclass
class TileUnit(Tileable):
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
            self.make_regularised_tile_from_elements()
        if self.regularised_tile is None:
            self.make_regularised_tile_from_elements()


    def setup_tile_and_elements(self) -> None:
        """Delegates setup of the unit to various functions depending
        on self.tiling_type.
        """
        if self.tiling_type == "cairo":
            _tiling_geometries.setup_cairo(self)
        elif self.tiling_type == "hex-dissection":
            _tiling_geometries.setup_hex_dissection(self)
        elif self.tiling_type == "laves":
            _tiling_geometries.setup_laves(self)
        elif self.tiling_type == "archimedean":
            _tiling_geometries.setup_archimedean(self)
        elif self.tiling_type in ("hex-colouring", "hex-coloring"):
            _tiling_geometries.setup_hex_colouring(self)
        elif self.tiling_type in ("square-colouring", "square-coloring"):
            _tiling_geometries.setup_square_colouring(self)
        else:
            _tiling_geometries.setup_none_tile(self)
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
        changes the tile to a diamond by copying and joining a 180 degree
        rotated copy. Operation is 'in-place'.
        """
        tile = self.tile.geometry[0]
        # translate to sit on x-axis
        tile = affine.translate(tile, 0, -tile.bounds[1])
        # make rotated copies
        # buffering applied to ensure union 'sticks'
        twins = [affine.rotate(tile, a, origin = (0, 0)) \
                    .buffer(self.fudge_factor, resolution = 1, join_style = 2) 
                 for a in range(0, 360, 180)]
        # and here we undo the buffer
        merged_tile = gpd.GeoSeries(twins).unary_union.buffer(
                -self.fudge_factor, resolution = 1, join_style = 2)
        self.tile_shape = TileShape.DIAMOND
        self.tile.geometry = gpd.GeoSeries([merged_tile])
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
    def inset_tile(self, d:float = 0) -> None:
        inset_tile = self.regularised_tile.geometry.buffer(
                -d, resolution = 1, join_style = 2)[0]
        # the clean_geometry seems needed to stop proliferation of vertices
        new_elements = [tiling_utils.clean_polygon(inset_tile.intersection(e))
                        for e in self.elements.geometry]
        self.elements.geometry = gpd.GeoSeries(new_elements)
        return
    
    
    def scale_elements(self, sf:float = 1) -> None:
        self.elements.geometry = self.elements.geometry.scale(
            sf, sf, origin = (0, 0))
        return


