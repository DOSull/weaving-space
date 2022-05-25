#!/usr/bin/env python
# coding: utf-8

from typing import Union
from dataclasses import dataclass
import itertools
import copy

import numpy as np
import geopandas as gpd
import pandas as pd

import shapely.affinity as affine
import shapely.geometry as geom
import shapely.wkt as wkt

from weave_units import WeaveUnit
from tile_units import Tileable
from tile_units import TileShape

import tiling_utils


@dataclass
class TileGrid:
    tile:gpd.GeoSeries = None
    to_tile:gpd.GeoSeries = None
    grid_type:TileShape = None
    extent:gpd.GeoSeries = None
    centre:tuple[float] = None
    points:gpd.GeoSeries = None
    
    def __init__(self, tile:gpd.GeoSeries, to_tile:gpd.GeoSeries, 
                 grid_type:TileShape = TileShape.RECTANGLE, 
                 to_hex:bool = True) -> None:
        self.grid_type = grid_type
        self.tile = tile
        if self.grid_type == TileShape.TRIANGLE:
            self.tile, self.grid_type = self._modify_triangle_tile(to_hex)
        self.to_tile = gpd.GeoSeries([to_tile.unary_union])
        self.extent, self.centre = self._get_extent()
        self.points = self._get_points()
        
    
    def _get_extent(self) -> gpd.GeoSeries:
        mrr = self.to_tile.geometry[0].minimum_rotated_rectangle
        mrr_centre = geom.Point(mrr.centroid.coords[0])
        mrr_corner = geom.Point(mrr.exterior.coords[0])
        radius = mrr_centre.distance(mrr_corner)
        return gpd.GeoSeries([mrr_centre.buffer(radius)]), mrr_centre
    
        
    def _get_points(self) -> gpd.GeoSeries:
        if self.grid_type in (TileShape.RECTANGLE, ):
            pts = self._get_rect_centres()
        elif self.grid_type in (TileShape.HEXAGON, TileShape.TRIHEX):
            pts = self._get_hex_centres()
        elif self.grid_type in (TileShape.TRIDIAMOND, ):
            pts = self._get_diamond_centres()
        tr = affine.translate  # for efficiency here
        tiles = [tr(self.tile.geometry[0], p[0], p[1]) 
                 for p in list(pts)]
        tiles = [t for t in tiles if self.extent[0].intersects(t)]
        return gpd.GeoSeries([t.centroid for t in tiles])
    
        

    def _get_width_height_left_bottom(self, 
                                      gs:gpd.GeoSeries
                                    ) -> tuple[float]:
        """Returns width, height, left and bottom limits of a GeoSeries

        Args:
            gs (geopandas.GeoSeries): GeoSeries for which limits are required.

        Returns:
            tuple: four float values of width, height, left and bottom of gs.
        """    
        extent = gs.total_bounds
        return extent[2] - extent[0], extent[3] - extent[1], extent[0], extent[1]


    def _get_grid(self, ll: tuple[float], nums: tuple[int], 
                  tdim: tuple[float]) -> np.ndarray:
        """Returns rectilinear grid of x,y coordinate pairs.

        Args:
            ll (tuple[float]): lower left corner coordinates of the grid as 
                (x, y). 
            nums (tuple[int]): grid extent as (number of columns, number 
                of rows).
            tdim (tuple[float]): grid resolution as (column width, column 
                height)

        Returns:
            np.ndarray: a matrix of nums[0] * nums[1] rows and 2 columns, 
                each row
            containing an x, y coordinate pair.
        """    
        return np.array(np.meshgrid(np.arange(nums[0]) * tdim[0] + ll[0],
                                    np.arange(nums[1]) * tdim[1] + ll[1])
                        ).reshape(2, nums[0] * nums[1]).transpose()
        

    def _get_rect_centres(self) -> np.ndarray:
        """Returns a rectangular grid of translation vectors that will 'fill' to_tile_gs polygon with the tile_gs polygon (which should be rectangular).

        Returns:
            np.ndarray: A 2 column array each row being an x, y translation vector.
        """    
        tt_w, tt_h, tt_x0, tt_y0 = \
            self._get_width_height_left_bottom(self.extent)
        tile_w, tile_h, tile_x0, tile_y0 = \
            self._get_width_height_left_bottom(self.tile)
        # number of tiles in each direction
        nx = int(np.ceil(tt_w / tile_w))
        ny = int(np.ceil(tt_h / tile_h))
        # origin is inset from the lower lwft corner based
        x0 = (tt_w - (nx * tile_w)) / 2 + tt_x0
        y0 = (tt_h - (ny * tile_h)) / 2 + tt_y0
        return self._get_grid((x0, y0), (nx + 1, ny + 1), (tile_w, tile_h))


    def _get_hex_centres(self) -> np.ndarray:
        """Returns a hexagonal grid of translation vectors that will 'fill' 
        to_tile_gs with the tile_gs polygon (which should be hexagonal).
        
        Returns:
            np.ndarray: A 2 column array each row being an x, y translation vector.
        """    
        tt_w, tt_h, tt_x0, tt_y0  = \
            self._get_width_height_left_bottom(self.extent)
        tile_w, tile_h, tile_x0, tile_y0 = \
            self._get_width_height_left_bottom(self.tile)
        nx = int(np.ceil(tt_w / (tile_w * 3 / 2)))
        ny = int(np.ceil(tt_h / tile_h))
        # the effective width of two columns of hexagonal tiles is 3w/2
        x0 = (tt_w - (nx * tile_w * 3 / 2)) / 2 + tt_x0
        y0 = (tt_h - (ny * tile_h)) / 2 + tt_y0
        # get two offset rectangular grids and combine them
        g1 = self._get_grid((x0, y0 + tile_h / 4), 
                            (nx + 1, ny + 1), 
                            (tile_w * 3 / 2, tile_h))
        g2 = self._get_grid((x0 + tile_w * 3 / 4, y0 - tile_h / 4), 
                            (nx, ny), 
                            (tile_w * 3 / 2, tile_h))
        return np.append(g1, g2).reshape((g1.shape[0] + g2.shape[0], 2))

    
    # Actually returns rhombus centres
    # 
    #     /\
    #    /  \
    #    \  /
    #     \/
    # 
    def _get_diamond_centres(self) -> np.ndarray:
        tt_w, tt_h, tt_x0, tt_y0  = \
            self._get_width_height_left_bottom(self.extent)
        tile_w, tile_h, tile_x0, tile_y0 = \
            self._get_width_height_left_bottom(self.tile)
        nx = int(np.ceil(tt_w / tile_w))
        ny = int(np.ceil(tt_h / tile_h))
        x0 = (tt_w - (nx * tile_w)) / 2 + tt_x0
        y0 = (tt_h - (ny * tile_h)) / 2 + tt_y0
        g1 = self._get_grid((x0, y0), 
                            (nx + 1, ny + 1), 
                            (tile_w, tile_h))
        g2 = self._get_grid((x0 - tile_w / 2, y0 - tile_h / 2),
                            (nx + 1, ny + 1), 
                            (tile_w, tile_h))
        return np.append(g1, g2).reshape(g1.shape[0] + g2.shape[0], 2)


@dataclass
class Tiling:
    tile_unit:Union[WeaveUnit, Tileable] = None
    tile_shape:str = ""
    region:gpd.GeoDataFrame = None
    region_id_var:str = None
    grid:TileGrid = None
    tiles:gpd.GeoDataFrame = None

    def __init__(self, unit:Union[WeaveUnit, Tileable], 
                 region:gpd.GeoDataFrame, id_var:str) -> None:
        self.tile_shape = unit.tile_shape
        self.tile_unit = unit
        self.region = region
        self.region_id_var = ("ID" if id_var is None else id_var)
        self.grid = TileGrid(self.tile_unit.tile.geometry,
                             self.region.geometry, self.tile_shape)
        self.tiles = self.make_tiling()


    def get_tiled_map(self, id_var:str = None, rotation:float = 0.,             
                      prioritise_tiles:bool = False) -> gpd.GeoDataFrame:
        id_var = (self.region_id_var
                  if id_var is None
                  else id_var)
        tiled_map = self.rotated(rotation)
        region_vars = list(self.region.columns)
        region_vars.remove("geometry")
        region_vars.remove(id_var)
        
        if prioritise_tiles: # respect tile sides over zone boundaries
            # make column with unique ID for every element in the tiling
            tiled_map["tileUID"] = list(range(tiled_map.shape[0]))
            # overlay with the zones from the region to be tiled
            tiled_map = tiled_map.overlay(self.region, keep_geom_type = False)
            # determine areas of overlaid tile elements and drop the data
            tiled_map["area"] = tiled_map.geometry.area
            tiled_map = tiled_map.drop(columns = region_vars)
            # make a lookup by largest area element to the zone ID
            lookup = tiled_map.iloc[tiled_map.groupby("tileUID")["area"].agg(
                pd.Series.idxmax)][["tileUID", id_var]]
            tiled_map = tiled_map.drop(columns = [id_var])
            # now join the lookup and from there the region data
            tiled_map = tiled_map \
                .merge(lookup, on = "tileUID")   \
                .merge(self.region.drop(columns = ["geometry"]), on = id_var) 
        else:
            tiled_map = self.region.overlay(tiled_map)
        
        # make a dissolve variable from element_id and id_var
        tiled_map["diss_var"] = (tiled_map.element_id + 
                                 tiled_map[id_var].astype(str))
        return tiled_map \
            .dissolve(by = "diss_var", as_index = False) \
            .drop(["diss_var"], axis = 1)
    
    
    def _rotate_gdf_to_geoseries(
            self, gdf:gpd.GeoDataFrame, 
            angle:float, centre:tuple = (0, 0)
        ) -> tuple[gpd.GeoSeries, tuple[float]]:
        """Rotates the geometries in a GeoDataFrame as a single collection.
        
        Rotation is about the supplied centre (if supplied) or about the centroid of the GeoDataFrame (if not). This allows for reversal of 
        a rotation. [Note that this might not be a required precaution!]

        Args:
            gdf (geopandas.GeoDataFrame): GeoDataFrame to rotate
                angle (float): angle of rotation (degrees).
            centre (tuple, optional): desired centre of rotation. Defaults 
                to (0, 0).

        Returns:
            tuple: a geopandas.GeoSeries and a tuple (point) of the centre of 
                the rotation.
        """    
        centre = (
            gdf.geometry.unary_union.centroid.coords[0] 
            if centre is None 
            else centre)
        return gdf.geometry.rotate(angle, origin = centre), centre


    def make_tiling(self) -> gpd.GeoDataFrame:
        """Tiles the region with a weave unit tile, returning a GeoDataFrame

        Returns:
            geopandas.GeoDataFrame: a GeoDataFrame of the region tiled with the
                weave unit.
        """
        # we assume the geometry column is called geometry so make it so...
        if self.region.geometry.name != "geometry":
            self.region.rename_geometry("geometry", inplace = True)

        # chain list of lists of GeoSeries geometries to list of geometries 
        # previously used self._translate_geoms, but I don't think we need it
        tiles = itertools.chain(*[
            self.tile_unit.elements.geometry.translate(p.x, p.y)
            for p in self.grid.points])
        # replicate the element ids
        ids = list(self.tile_unit.elements.element_id) * len(self.grid.points)
        tiles_gs = gpd.GeoSeries(tiles)
        tiles_gdf = gpd.GeoDataFrame(data = {"element_id": ids},
                                     geometry = tiles_gs, 
                                     crs = self.tile_unit.crs)
        # assemble and return as a GeoDataFrame
        tiles_gdf.geometry = tiling_utils.gridify(tiles_gdf.geometry)
        return tiles_gdf
        
    
    def rotated(self, rotation:float = None):
        if self.tiles is None:
            self.tiles = self.make_tiling()
        if rotation is None or rotation == 0:
            return self.tiles
        return gpd.GeoDataFrame(
            data = {"element_id": self.tiles.element_id}, crs = self.tiles.crs,
            geometry = self.tiles.geometry.rotate(rotation, 
                                                  origin = self.grid.centre))
        
    
    def get_local_patch(self, n:int = 3) -> gpd.GeoDataFrame:
        patch_extent = gpd.GeoDataFrame(
            data = {"id": [1]}, crs = self.tile_unit.crs,
            geometry = self.tile_unit.tile.geometry.scale(n, n, origin = (0, 0))
        )
        patch = Tiling(self.tile_unit, patch_extent, id_var = "id")
        tiles_to_keep = []
        ids = []
        for p, id in zip(patch.tiles.geometry, patch.tiles.element_id):
            if p.within(patch_extent.geometry[0]):
                tiles_to_keep.append(p)
                ids.append(id)
        return gpd.GeoDataFrame(
            data = {"element_id": ids}, crs = self.tile_unit.crs,
            geometry = gpd.GeoSeries(tiles_to_keep))
    
    
