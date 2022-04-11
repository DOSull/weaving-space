#!/usr/bin/env python
# coding: utf-8

from typing import Union
from dataclasses import dataclass
import itertools
import logging

import numpy as np
import geopandas as gpd
import shapely.affinity as affine
import shapely.geometry as geom

from triaxial_weave_units import get_triaxial_weave_unit
from biaxial_weave_units import get_biaxial_weave_unit


@dataclass
class WeaveUnit:
    """ Small data class containing elements of a weave unit.
    
    Attributes:
        elements: a GeoDataFrame of strand geometries.
        tile: a GeoDataFrame of the weave_unit tileable polygon (either a
            rectangle or a hexagon).
    """  
    elements:gpd.GeoDataFrame = None
    tile:gpd.GeoDataFrame = None
    tile_shape:str = "rectangle"
    weave_type:str = "plain"
    spacing:float = 10_000.
    aspect:float = 1.
    margin:float = 0.
    n:Union[int, tuple[int]] = (2, 2)
    strands:str = "a|b|c"
    tie_up:np.ndarray = None
    tr:np.ndarray = None
    th:np.ndarray = None
    crs:int = 3857
    
    def __init__(self, **kwargs):
        """Constructor for WeaveUnit. Parameters are passed through to get_weave_unit function and also stored as instance attributes.
        
        Args:
            weave_type (str, optional): the type of weave pattern, one of 
                "plain",  "twill", "basket", "this", "cube" or "hex". Defaults
                to "plain".
            spacing (float, optional): spacing of threads in the weave in units 
                of the CRS. Defaults to 10000.
            aspect (float, optional): width of strands relative to the spacing. 
                Defaults to 1.
            margin (float, optional): margin applied to 'shrink' strand 
                elements, relative to the spacing. Defaults to 0.
            n (tuple of ints): number of over-under strands in biaxial weaves. 
                Only one item is required in a plain weave. Twill and basket patterns expect an even number of elements in the tuple. Defaults to (2, 2).
            strands (str, optional): specification of the strand labels 
                along each axis. Defaults to "a|b|c".
            tie_up (numpy.ndarray, optional): used when type is "this" to
                specify a desired weave pattern. See: Glassner A, 2002, “Digital weaving. 1” IEEE Computer Graphics and Applications 22(6) 108–118 DOI: 10.1109/MCG.2002.1046635. Defaults to None.
            tr (numpy.ndarray, optional): used when type is "this" to specify 
                the treadling matrix. See: Glassner 2002. Defaults to None.
            th (numpy.ndarray, optional): used when type is "this" to specify
                the threading matrix. See: Glassner 2002. Defaults to None.
            crs (int, optional): coordinate reference system. Usually an integer
                EPSG code, but any CRS object interpretable by geopandas will
                work. Defaults to 3857 (for Web Mercator).
        """
        unit = self._get_weave_unit(**kwargs)
        self.elements = unit["weave_unit"]
        self.tile = unit["tile"]
        for k, v in kwargs.items():
            self.__dict__[k] = v
        self.tile_shape = ("hexagon" 
                           if self.weave_type in ("hex", "cube")
                           else "rectangle") 


    def _get_weave_unit(self, weave_type:str = "plain", spacing:float = 10000, 
            aspect:float = 1, margin:float = 0, n:tuple[int] = (2, 2), 
            strands:str = "a|b|c", tie_up:np.ndarray = None, 
            tr:np.ndarray = None, th:np.ndarray = None, crs:int = 3857
        ) -> dict:
        """Returns dictionary with weave unit and tile GeoDataFrames
        
        Args:
            weave_type (str, optional): the type of weave pattern, one of 
                "plain",  "twill", "basket", "this", "cube" or "hex". Defaults
                to "plain".
            spacing (float, optional): spacing of threads in the weave in units 
                of the CRS. Defaults to 10000.
            aspect (float, optional): width of strands relative to the spacing. 
                Defaults to 1.
            margin (float, optional): margin applied to 'shrink' strand 
                elements, relative to the spacing. Defaults to 0.
            n (tuple of ints): number of over-under strands in biaxial weaves. 
                Only one item is required in a plain weave. Twill and basket patterns expect an even number of elements in the tuple. Defaults to (2, 2).
            strands (str, optional): specification of the strand labels 
                along each axis. Defaults to "a|b|c".
            tie_up (numpy.ndarray, optional): used when type is "this" to
                specify a desired weave pattern. See: Glassner A, 2002, “Digital weaving. 1” IEEE Computer Graphics and Applications 22(6) 108–118 DOI: 10.1109/MCG.2002.1046635. Defaults to None.
            tr (numpy.ndarray, optional): used when type is "this" to specify 
                the treadling matrix. See: Glassner 2002. Defaults to None.
            th (numpy.ndarray, optional): used when type is "this" to specify
                the threading matrix. See: Glassner 2002. Defaults to None.
            crs (int, optional): coordinate reference system. Usually an integer
                EPSG code, but any CRS object interpretable by geopandas will
                work. Defaults to 3857 (for Web Mercator).

        Returns:
            dict: dictionary with contents {"weave_unit": GeoDataFrame of weave 
                elements, "tile": GeoDataFrame of the tile}.
        """
        self._parameter_info(margin, aspect)

        if weave_type in ("hex", "cube"):
            unit = get_triaxial_weave_unit(weave_type = weave_type,
                                        spacing = spacing, aspect = aspect,
                                        margin = margin, strands = strands, 
                                        crs = crs)
        else:
            unit = get_biaxial_weave_unit(weave_type = weave_type, n = n,
                                        spacing = spacing, aspect = aspect,
                                        margin = margin, strands = strands, 
                                        crs = crs, tie_up = tie_up, )
        return unit


    def _parameter_info(self, margin: float, aspect: float) -> None:
        """Outputs logging message concerning the supplied margin and aspect settings.

        Args:
            margin (float): weave unit margin.
            aspect (float): weave unit aspect.
        """    
        
        if aspect == 0:
            logging.info("""Setting aspect to 0 is probably not a great plan.""")

        if aspect < 0 or aspect > 1:
            logging.warning("""Values of aspect outside the range 0 to 1 won't 
                            produce tiles that will look like weaves, but they might be pretty anyway! Values less than -1 seem particularly promising, especially with opacity set less than 1.""")

        # maximum margin that will produce a weave-able tile
        max_margin = (1 - aspect) / 2
        if margin > max_margin:
            logging.warning(f"""With aspect set to {aspect:.3f} the largest margin 
                            that will work is {max_margin:.3f}. Lower values are required to produce proper tileable weaves. Specifically, with too wide a margin, strands in adjacent tiles will not 'join up' when tiled. Higher values will make nice tilings with broken strands, which aren't 'proper' weaves. The best alternative is to make the weave unit with margin = 0, then apply a negative buffer after you have tiled your map.""")   
        return None
  

@dataclass
class TileGrid:
    tile:gpd.GeoSeries = None
    to_tile:gpd.GeoSeries = None
    is_hex_grid:bool = None
    extent:gpd.GeoSeries = None
    centre:tuple[float] = None
    points:gpd.GeoSeries = None
    
    def __init__(self, tile, to_tile, hexes:bool = False):
        self.tile = tile
        self.to_tile = gpd.GeoSeries([to_tile.unary_union])
        self.is_hex_grid = hexes
        self.extent, self.centre = self._get_extent()
        self.points = self._get_points()
        
    
    def _get_extent(self) -> gpd.GeoSeries:
        mrr = self.to_tile.geometry[0].minimum_rotated_rectangle
        mrr_centre = geom.Point(mrr.centroid.coords[0])
        mrr_corner = geom.Point(mrr.exterior.coords[0])
        radius = mrr_centre.distance(mrr_corner)
        return gpd.GeoSeries([mrr_centre.buffer(radius)]), mrr_centre
    
        
    def _get_points(self) -> gpd.GeoSeries:
        pts = (self._get_hex_centres()
               if self.is_hex_grid
               else self._get_rect_centres())
        tr = affine.translate
        tiles = [ tr(self.tile.geometry[0], p[0], p[1]) 
                  for p in list(pts) ]
        tiles = [ t for t in tiles if self.extent[0].intersects(t) ]
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
    
    # EXPERIMENTAL
    # Based on 4 offset rectangular grids
    # 
    #     /\
    #    /__\
    #    \  /
    #     \/
    # 
    def _get_tri_centres(self) -> np.ndarray:
        tt_w, tt_h, tt_x0, tt_y0  = \
            self._get_width_height_left_bottom(self.extent)
        tile_w, tile_h, tile_x0, tile_y0 = \
            self._get_width_height_left_bottom(self.tile)
        nx = int(np.ceil(tt_w / tile_w))
        ny = int(np.ceil(tt_h / tile_h / 2))
        x0 = (tt_w - (nx * tile_w)) / 2 + tt_x0
        y0 = (tt_h - (ny * tile_h * 2)) / 2 + tt_y0
        # this time there are four of them
        g1 = self._get_grid((x0, y0), 
                            (nx + 1, ny + 1), 
                            (tile_w, tile_h * 2))
        g2 = self._get_grid((x0, y0 + tile_h * 2 / 3),
                            (nx + 1, ny + 1), 
                            (tile_w, tile_h * 2))
        g3 = self._get_grid((x0 + tile_w / 2, y0 - tile_h / 3),
                            (nx + 1, ny + 1), 
                            (tile_w, tile_h * 2))
        g4 = self._get_grid((x0 + tile_w / 2, y0 + tile_h),
                            (nx + 1, ny + 1), 
                            (tile_w, tile_h * 2))
        return np.append(np.append(np.append(g1, g2), g3), g4).reshape(
            g1.shape[0] + g2.shape[0] + g3.shape[0] + g4.shape[0], 2)


class Tiling:
    tile:WeaveUnit = None
    region:gpd.GeoDataFrame = None
    grid:TileGrid = None
    tiles:gpd.GeoDataFrame = None

    def __init__(self, unit:WeaveUnit, 
                 region:gpd.GeoDataFrame) -> None:
        self.tile = unit
        self.region = region
        self.grid = TileGrid(self.tile.tile.geometry,
                             self.region.geometry, 
                             self.tile.tile_shape == "hexagon")
        self.tiles = self.make_tiling()


    def _translate_geoms(self, gs:gpd.GeoSeries, 
                dx:float = 0., dy:float = 0.
            ) -> list[Union[geom.Polygon, geom.MultiPolygon]]:
        """Translates geometries in supplied GeoSeries by (dx, dy).
        
        This is needed in place of GeoSeries.translate because we have 
        to unpack the geometries to a list so that we can chain them into a list and bundle back up into a GeoSeries in the get_tiling function.

        Args:
            gs (geopandas.GeoSeries): GeoSeries to translate.
            dx (float, optional): x transation. Defaults to 0.
            dy (float, optional): x transation. Defaults to 0.

        Returns:
            list[Polygon|MultiPolygon]: _description_
        """ 
        tr = affine.translate   
        return [ tr(s, dx, dy) for s in gs ]


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
        tiles = itertools.chain(*[self._translate_geoms(
                                    self.tile.elements.geometry, p.x, p.y) 
                        for p in self.grid.points])
        # replicate the strand ids
        ids = list(self.tile.elements.strand) * len(self.grid.points)
        tiles_gs = gpd.GeoSeries(tiles)
        
        # assemble and return as a GeoDataFrame
        return gpd.GeoDataFrame(data = {"strand": ids},
                                      geometry = tiles_gs, 
                                      crs = self.tile.crs)
        
    
    def rotated(self, rotation:float = None):
        if self.tiles is None:
            self.tiles = self.make_tiling()
        if rotation is None or rotation == 0:
            return self.tiles
        return gpd.GeoDataFrame(
            data = {"strand": self.tiles.strand}, crs = self.tiles.crs,
            geometry = self.tiles.geometry.rotate(rotation, 
                                                  origin = self.grid.centre))