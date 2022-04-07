#!/usr/bin/env python
# coding: utf-8

from itertools import chain
from dataclasses import dataclass
import logging

import numpy as np
import geopandas 
from shapely.affinity import translate
from shapely.geometry import MultiPolygon
from shapely.geometry import Polygon

from triaxial_weave_units import get_triaxial_weave_unit
from biaxial_weave_units import get_biaxial_weave_unit


@dataclass
class WeaveUnit:
    """ Small data class containing elements of a weave unit.
    
    Attributes:
        elements: a GeoDataFrame of strand geometries.
        tile: a GeoDataFrame of the weave_unit tileable polygon (either a
            rectangle or a hexagon).
        transform: a numpy.ndarray containing any affine transform required to
            tile the unit as a square.
        type: a string indicating the type of weave, one of "plain", "twill",
            "basket", "this", "cube" or "hex".
    """  
    elements:geopandas.GeoDataFrame
    tile:geopandas.GeoDataFrame
    transform:np.ndarray = np.diag(np.ones(2))
    type:str = "plain"
    
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
        self.transform = np.diag(np.ones(2))
        for k, v in kwargs.items():
            self.__dict__[k] = v


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
            unit = get_biaxial_weave_unit(weave_type = weave_type, 
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
  

def _get_width_height_left_bottom(gs:geopandas.GeoSeries) -> tuple[float]:
    """Returns width, height, left and bottom limits of a GeoSeries

    Args:
        gs (geopandas.GeoSeries): GeoSeries for which limits are required.

    Returns:
        tuple: four float values of width, height, left and bottom of gs.
    """    
    extent = gs.total_bounds
    return extent[2] - extent[0], extent[3] - extent[1], extent[0], extent[1]


def _get_grid(ll: tuple[float], nums: tuple[int], 
              tdim: tuple[float]) -> np.ndarray:
    """Returns rectilinear grid of x,y coordinate pairs.

    Args:
        ll (tuple[float]): lower left corner coordinates of the grid as (x, y).
        nums (tuple[int]): grid extent as (number of columns, number of rows).
        tdim (tuple[float]): grid resolution as (column width, column height)

    Returns:
        np.ndarray: a matrix of nums[0] * nums[1] rows and 2 columns, each row
        containing an x, y coordinate pair.
    """    
    return np.array(np.meshgrid(np.arange(nums[0]) * tdim[0] + ll[0],
                                np.arange(nums[1]) * tdim[1] + ll[1])
                    ).reshape(2, nums[0] * nums[1]).transpose()
      

def _get_rect_centres(tile_gs:geopandas.GeoSeries, 
                      to_tile_gs:geopandas.GeoSeries) -> np.ndarray:
    """Returns a rectangular grid of translation vectors that will 'fill' to_tile_gs polygon with the tile_gs polygon (which should be rectangular).

    Args:
        tile_gs (geopandas.GeoSeries): the tile (which should be a rectangle).
        to_tile_gs (geopandas.GeoSeries): the region to be tiled.

    Returns:
        np.ndarray: A 2 column array each row being an x, y translation vector.
    """    
    tt_w, tt_h, tt_x0, tt_y0 = _get_width_height_left_bottom(to_tile_gs)
    tile_w, tile_h, tile_x0, tile_y0 = _get_width_height_left_bottom(tile_gs)
    nx = int(np.ceil(tt_w / tile_w))
    ny = int(np.ceil(tt_h / tile_h))
    x0 = ((nx * tile_w) - tt_w) / 2 + tile_x0 + tt_x0
    y0 = ((ny * tile_h) - tt_h) / 2 + tile_y0 + tt_y0
    return _get_grid((x0, y0), (nx, ny), (tile_w, tile_h))


def _get_hex_centres(tile_gs:geopandas.GeoSeries, 
                     to_tile_gs:geopandas.GeoSeries) -> np.ndarray:
    """Returns a hexagonal grid of translation vectors that will 'fill' 
    to_tile_gs with the tile_gs polygon (which should be hexagonal).

    Args:
        tile_gs (geopandas.GeoSeries): the tile (which should be hexagonal).
        to_tile_gs (geopandas.GeoSeries): the region to be tiled.

    Returns:
        np.ndarray: A 2 column array each row being an x, y translation vector.
    """    
    tt_w, tt_h, tt_x0, tt_y0  = _get_width_height_left_bottom(to_tile_gs)
    tile_w, tile_h, tile_x0, tile_y0 = _get_width_height_left_bottom(tile_gs)
    nx = int(np.ceil(tt_w / (tile_w * 3 / 2))) + 1
    ny = int(np.ceil(tt_h / tile_h)) + 1
    # the effective width of two columns of hexagonal tiles is 3w/2
    x0 = ((nx * tile_w * 3 / 2) - tt_w) / 2 + tile_x0 + tt_x0
    y0 = ((ny * tile_h) - tt_h) / 2 + tile_y0 + tt_y0
    # get two offset rectangular grids and combine them
    g1 = _get_grid((x0, y0 + tile_h / 4), 
                   (nx, ny), 
                   (tile_w * 3 / 2, tile_h))
    g2 = _get_grid((x0 + tile_w * 3 / 4, y0 - tile_h / 4), 
                   (nx, ny), 
                   (tile_w * 3 / 2, tile_h))
    return np.append(g1, g2).reshape((g1.shape[0] + g2.shape[0], 2))


def _translate_geoms(gs:geopandas.GeoSeries, 
                     dxy:tuple = (0, 0)) -> list[Polygon|MultiPolygon]:
    """Translates geometries in supplied GeoSeries by the vector dxy.
    
    This is needed in place of GeoSeries.translate because we have to unpack 
    the geometries to a list so that we can chain them into a list and bundle
    back up into a GeoSeries in the get_tiling function.

    Args:
        gs (geopandas.GeoSeries): GeoSeries to translate.
        dxdy (tuple, optional): _description_. Defaults to (0, 0).

    Returns:
        list[Polygon|MultiPolygon]: _description_
    """    
    return [translate(s, dxy[0], dxy[1]) for s in gs]


def _rotate_gdf_to_geoseries(
                    gdf:geopandas.GeoDataFrame, angle:float, 
                    centre:tuple = None) -> tuple[geopandas.GeoSeries, 
                                                 tuple[float]]:
    """Rotates the geometries in a GeoDataFrame as a single collection.
    
    Rotation is about the supplied centre (if supplied) or about the centroid of the GeoDataFrame (if not). This allows for reversal of a rotation. 
    [Note that this might not be a required precaution!]

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to rotate
        angle (float): angle of rotation (degrees).
        about (tuple, optional): desired centre of rotation. Defaults to None.

    Returns:
        tuple: a geopandas.GeoSeries and a tuple (point) of the centre of the
            rotation.
    """    
    centre = (
        gdf.geometry.unary_union.centroid.coords[0] 
        if centre is None 
        else centre)
    return gdf.geometry.rotate(angle, origin = centre), centre


# takes WeaveUnit and region GeoDataFrame and tiles the units over the region
def get_tiling(unit: WeaveUnit, 
               region:geopandas.GeoDataFrame, 
               rotation = 0) -> geopandas.GeoDataFrame:
    """Tiles the region with a weave unit tile, returning a GeoDataFrame

    Args:
        unit (WeaveUnit): WeaveUnit containing the tile.
        region (geopandas.GeoDataFrame): the region to tile.
        rotation (int, optional): a desired rotation. Defaults to 0.

    Returns:
        geopandas.GeoDataFrame: a GeoDataFrame of the region tiled with the
            weave unit.
    """
    # we assume that the geometry column is called geometry, so make it so...
    if region.geometry.name != "geometry":
        region.rename_geometry("geometry", inplace = True)

    tile_gs = unit.tile.geometry
    if rotation != 0:
        to_tile_gs, centre = _rotate_gdf_to_geoseries(region, rotation)
    else:
        to_tile_gs, centre = region.geometry, None
        
    if unit.weave_type in ("hex", "cube"):
        shifts = _get_hex_centres(tile_gs, to_tile_gs)
    else:
        shifts = _get_rect_centres(tile_gs, to_tile_gs)

    # TODO: filter shifts by some buffered version of the region, so as to save 
    # time (possibly) on generating and then clipping the tiled weave units

    # chain the list of lists of GeoSeries geometries to a list of geometries 
    tiles = chain(*[_translate_geoms(unit.elements.geometry, dxdy) 
                    for dxdy in shifts])
    # similarly replicate the list of strand ids
    ids = list(unit.elements.strand) * shifts.shape[0]
    
    # assemble into a GeoSeries and reverse rotation if required
    tiles_gs = geopandas.GeoSeries(tiles)
    if centre is not None:
        tiles_gs = tiles_gs.rotate(-rotation, origin = centre)
    
    # assemble and return as a GeoDataFrame
    return geopandas.GeoDataFrame(data = {"strand": ids},
                                  geometry = tiles_gs, 
                                  crs = unit.crs)