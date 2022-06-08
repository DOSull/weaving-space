#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
import itertools
from typing import Any
from typing import Union

import numpy as np
import geopandas as gpd
import pandas as pd

import matplotlib
import matplotlib.colors
import matplotlib.pyplot as pyplot

import shapely.affinity as affine
import shapely.geometry as geom

from tile_units import Tileable
from tile_units import TileUnit
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
        elif self.grid_type in (TileShape.HEXAGON, ):
            pts = self._get_hex_centres()
        elif self.grid_type in (TileShape.DIAMOND, ):
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
        return (extent[2] - extent[0], extent[3] - extent[1], 
                extent[0], extent[1])


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
        """Returns a rectangular grid of translation vectors that will 'fill'
        to_tile_gs polygon with the tile_gs polygon (which should be
        rectangular).

        Returns:
            np.ndarray: A 2 column array each row being an x, y translation
            vector.
        """    
        tt_w, tt_h, tt_x0, tt_y0 = \
            tiling_utils.get_width_height_left_bottom(self.extent)
        tile_w, tile_h, tile_x0, tile_y0 = \
            tiling_utils.get_width_height_left_bottom(self.tile)
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
            np.ndarray: A 2 column array each row being an x, y translation
            vector.
        """    
        tt_w, tt_h, tt_x0, tt_y0  = \
            tiling_utils.get_width_height_left_bottom(self.extent)
        tile_w, tile_h, tile_x0, tile_y0 = \
            tiling_utils.get_width_height_left_bottom(self.tile)
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
        """Reurns an diamond grid of translation vectors that will 'fill'
        to_tile_gs polygon with the tile_gs polygon (which should be
        rectangular).

        Returns:
            np.ndarray: A 2 column array each row being an x, y translation
            vector.
        """
        tt_w, tt_h, tt_x0, tt_y0 = \
            tiling_utils.get_width_height_left_bottom(self.extent)
        tile_w, tile_h, tile_x0, tile_y0 = \
            tiling_utils.get_width_height_left_bottom(self.tile)
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
    tile_unit:Tileable = None
    tile_shape:str = ""
    region:gpd.GeoDataFrame = None
    region_id_var:str = None
    grid:TileGrid = None
    tiles:gpd.GeoDataFrame = None
    rotation:float = 0.

    def __init__(self, unit:Tileable, region:gpd.GeoDataFrame, 
                 id_var:str) -> None:
        """A class to persist a tiling by filling an area relative to 
        a region that is sufficient to apply the tiling at any rotation.

        Args:
            unit (Tileable): the tile_unit to use
            region (gpd.GeoDataFrame): the region to be tiled
            id_var (str): a unique identifier variable in the region
        """
        self.tile_shape = unit.tile_shape
        self.tile_unit = unit
        self.region = region
        self.region_id_var = ("ID" if id_var is None else id_var)
        self.grid = TileGrid(self.tile_unit.tile.geometry,
                             self.region.geometry, self.tile_shape)
        self.tiles = self.make_tiling()


    def get_tiled_map(self, id_var:str = None, rotation:float = 0.,             
                      prioritise_tiles:bool = False) -> gpd.GeoDataFrame:
        """Returns a geodataframe filling a region at the requested rotation.

        Args:
            id_var (str, optional): the variable the distinguishes areas in the
                region to be tiled. None will be overwritten by the variable    
                name set on initialisation of the Tiling. Defaults to None.
            rotation (float, optional): An optional rotation to apply. Defaults 
                to 0. orientatijnto 
            prioritise_tiles (bool, optional): When True tiles will not be 
                broken at boundaries in the region dataset. Defaults to False.

        Returns:
            gpd.GeoDataFrame: a GeoDataFrame that contains the source region
                data along with the tile unit element_id variable.
        """
        # if no id_var is supplied overwrite it with the class id_var
        id_var = (self.region_id_var if id_var is None else id_var)
        tiled_map = self.rotated(rotation)
        # compile a list of the variable names we are NOT going to change
        # i.e. everything except the geometry and the id_var
        region_vars = list(self.region.columns)
        region_vars.remove("geometry")
        region_vars.remove(id_var)
        
        if prioritise_tiles: # maintain tile continuity across zone boundaries
            # make column with unique ID for every element in the tiling
            tiled_map["tileUID"] = list(range(tiled_map.shape[0]))
            # overlay with the zones from the region to be tiled
            tiled_map = tiled_map.overlay(self.region)  
            # determine areas of overlaid tile elements and drop the data
            # we join the data back later, so dropping makes that easier
            tiled_map["area"] = tiled_map.geometry.area
            tiled_map = tiled_map.drop(columns = region_vars)
            # make a lookup by largest area element to the region id variable
            lookup = tiled_map.iloc[tiled_map.groupby("tileUID")["area"].agg(
                pd.Series.idxmax)][["tileUID", id_var]]
            # remove the id_var before we replace it with a new one
            tiled_map = tiled_map.drop(columns = [id_var])
            # now join the lookup and from there the region data
            tiled_map = tiled_map \
                .merge(lookup, on = "tileUID") \
                .merge(self.region.drop(columns = ["geometry"]), on = id_var) 
        else:
            tiled_map = self.region.overlay(tiled_map)
        
        # make a dissolve variable from element_id and id_var
        tiled_map["diss_var"] = (tiled_map.element_id + 
                                 tiled_map[id_var].astype(str))
        tiled_map = tiled_map \
            .dissolve(by = "diss_var", as_index = False) \
            .drop(["diss_var"], axis = 1)
        
        tm = TiledMap()
        tm.tiling = self
        tm.tiled_map = tiled_map
        return tm
    
    
    def _rotate_gdf_to_geoseries(
            self, gdf:gpd.GeoDataFrame, 
            angle:float, centre:tuple = (0, 0)
        ) -> tuple[gpd.GeoSeries, tuple[float]]:
        """Rotates the geometries in a GeoDataFrame as a single collection.
        
        Rotation is about the supplied centre (if supplied) or about the
        centroid of the GeoDataFrame (if not). This allows for reversal of 
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
        # assemble and return as a GeoDataFrame
        tiles_gdf = gpd.GeoDataFrame(data = {"element_id": ids},
                                     geometry = tiles_gs, 
                                     crs = self.tile_unit.crs)
        # unclear if we need the below or not...
        # tiles_gdf.geometry = tiling_utils.gridify(tiles_gdf.geometry)
        return tiles_gdf
        
    
    def rotated(self, rotation:float = None) -> gpd.GeoDataFrame:
        """Returns the stored tiling rotated.

        Args:
            rotation (float, optional): Rotation angle in degrees. 
                Defaults to None.

        Returns:
            gpd.GeoDataFrame: Rotated tiling.
        """
        if self.tiles is None:
            self.tiles = self.make_tiling()
        self.rotation = rotation
        if self.rotation == 0:
            return self.tiles
        return gpd.GeoDataFrame(
            data = {"element_id": self.tiles.element_id}, crs = self.tiles.crs,
            geometry = self.tiles.geometry.rotate(rotation, 
                                                  origin = self.grid.centre))


@dataclass
class TiledMap:
    tiling: Tiling = None
    tiled_map: gpd.GeoDataFrame = None
    legend_elements: gpd.GeoDataFrame = None
    legend_key_gdf: gpd.GeoDataFrame = None
    variables: dict[str] = None
    colourmaps: dict[str:Union[str,dict]] = None
    legend: bool = True
    legend_zoom: float = 1.0
    scheme: str = "equalinterval"
    k: int = 7
    bins: list[Any] = None
    figsize: tuple[float] = (20, 15)
    dpi: float = 72
        

    def render(self, **kwargs) -> tuple[pyplot.Axes]:
        to_remove = set()  # keep track of kwargs we use to setup TiledMap
        for k, v in kwargs.items():
            if k in self.__dict__:
                self.__dict__[k] = v
                to_remove.add(k)
        for k in to_remove:
            del kwargs[k]

        if self.legend:
            # this sizing stuff is rough and ready for now
            # and possibly forever... 
            reg_w, reg_h, *_ = \
                tiling_utils.get_width_height_left_bottom(
                    self.tiled_map.geometry)
            tile_w, tile_h, *_ = \
                tiling_utils.get_width_height_left_bottom(
                    self.tiling.tile_unit.elements.geometry)
            approx_count = reg_w / tile_w * reg_h / tile_h
            sf = np.sqrt(approx_count) / 2
            gskw = {"height_ratios": [sf * tile_h, reg_h],
                    "width_ratios": [reg_w + sf * tile_w, sf * tile_w]}

            fig, axes = pyplot.subplot_mosaic(
                [["map", "legend"],
                 ["map", "."]], 
                gridspec_kw = gskw, figsize = self.figsize, 
                layout = "tight", **kwargs)
        else:
            fig, axes = pyplot.subplots(
                1, 1, figsize = self.figsize, 
                layout = "tight", **kwargs)

        if self.variables is None:
            # get any floating point columns available
            default_columns = \
                self.tiled_map.select_dtypes(
                    include = ("float64", "int64")).columns
            self.variables = dict(zip(
                self.tiled_map.element_id.unique(), 
                list(default_columns)))
            print(f"No variables specified, picked the first {len(self.variables)} numeric ones available.")        
        elif isinstance(self.variables, (list, tuple)):
            self.variables = dict(zip(
                self.tiling.tile_unit.elements.element_id.unique(),
                self.variables))
            print(f"Only a list of variables specified, assigning to available element_ids.")
                    
        if self.colourmaps is None:
            self.colourmaps = {}
            for var in self.variables.values():
                if self.tiled_map[var].dtype == pd.CategoricalDtype:
                    self.colourmaps[var] = "tab20"
                    print(f"For categorical data, you should specify colour mapping explicitly.")
                else:
                    self.colourmaps[var] = "Reds"
        
        self.plot_map(axes, **kwargs)
        return fig
    
    
    def plot_map(self, axes, **kwargs):
        if self.legend:
            axes["map"].set_axis_off()
            self.plot_subsetted_gdf(axes["map"], self.tiled_map, **kwargs)
            self.get_legend(ax = axes["legend"], **kwargs)
        else:
            axes.set_axis_off()
            self.plot_subsetted_gdf(axes, self.tiled_map, **kwargs)
        return None
    
    
    def plot_subsetted_gdf(self, ax, gdf, **kwargs):
        groups = gdf.groupby("element_id")
        for id, var in self.variables.items():
            subset = groups.get_group(id)
            # Handle custom color assignments via 'cmaps' parameter.
            # Result is setting 'cmap' variable used in plot command afterwards.
            if (isinstance(self.colourmaps, 
                           (str, matplotlib.colors.Colormap,
                            matplotlib.colors.LinearSegmentedColormap,
                            matplotlib.colors.ListedColormap))):
                cmap = self.colourmaps   # user wants one palette for all ids
            elif (len(self.colourmaps) == 0):
                cmap = 'Reds'  # set a default... here, to Brewer's 'Reds'
            elif (var not in self.colourmaps):
                cmap = 'Reds'  # var has no color specified in dict, use default
            elif (isinstance(self.colourmaps[var],
                             (str, matplotlib.colors.Colormap,
                              matplotlib.colors.LinearSegmentedColormap,
                              matplotlib.colors.ListedColormap))):
                cmap = self.colourmaps[var]  # user speced colors for this var
            elif (isinstance(self.colourmaps[var], dict)):
                colormap_dict = self.colourmaps[var]
                data_unique_sorted = subset[var].unique()
                data_unique_sorted.sort()
                cmap = matplotlib.colors.ListedColormap(
                    [colormap_dict[x] for x in data_unique_sorted])
            else:
                raise Exception(f"Color map for '{var}' is not a known type, but is {str(type(self.colourmaps[var]))}")

            subset.plot(ax = ax, column = var, cmap = cmap, **kwargs)
    
    
    def to_file(self, fname:str = None) -> str:
        return fname
    
    
    def get_legend(self, ax: pyplot.Axes = None, **kwargs) -> None:
        # turn off axes (which seems also to make it impossible
        # to set a background colour)
        ax.set_axis_off()

        legend_elements = self.tiling.tile_unit._get_legend_elements()
        # this is a bit hacky, but we will apply the rotation at the end
        # so for TileUnits which don't need it, reverse that now
        if isinstance(self.tiling.tile_unit, TileUnit):
            legend_elements.rotation = -self.tiling.rotation
        
        legend_key = self.get_legend_key_gdf(legend_elements)
        # set a zoom 
        bb = [x / self.legend_zoom for x in legend_key.geometry.total_bounds]
        ax.set_xlim(bb[0], bb[2])
        ax.set_ylim(bb[1], bb[3])
        # not using this at the moment, but if we want to colour the 
        # background here is how when axes are set off
        # ax.axhspan(bb[1], bb[3], fc = "w", lw = 0)

        # now plot background; we include the central tiles, since in
        # the weave case these may not match the legend elements
        self.tiling.tile_unit.get_local_patch(
            r = 2, include_0 = True).geometry.rotate(
                self.tiling.rotation, origin = (0, 0)).plot(
                ax = ax, fc = "#9F9F9F3F", ec = "#5F5F5F", lw = 0.5)

        # plot the legend key elements (which include the data)
        self.plot_subsetted_gdf(ax, legend_key, lw = 0, **kwargs)
        
        # now add the annotations - for this we go back to the legend elements
        legend_elements.geometry = legend_elements.geometry.rotate(
            self.tiling.rotation, origin = (0, 0))
        
        for id, tile, rotn in zip(legend_elements.element_id,
                                  legend_elements.geometry,
                                  legend_elements.rotation):
            c = tile.centroid
            ax.annotate(self.variables[id], xy = (c.x, c.y), 
                    ha = "center", va = "center", rotation_mode = "anchor", 
                    # adjust rotation to favour text reading left to right
                    rotation = (rotn + self.tiling.rotation + 90) % 180 - 90, 
                    bbox = {"lw": 0, "fc": "#ffffff40"})


    def get_legend_key_gdf(self, elements:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """Returns a GeoDataFrame of tile elements dissected and with
        data assigned to the slice so a map of them can stand for a legend.
        The 'dissection' is handled differently by WeaveUnits and TileUnits
        and delegated to get_legend_key_shapes().

        Args:
            elements (gpd.GeoDataFrame): the legend elements

        Returns:
            gpd.GeoDataFrame:  with element_id, variables and rotation
                attributes, and geometries of Tileable elements sliced into a 
                colour ramp or set of nested tiles.
        """
        key_tiles = []   # set of tiles to form a colour key (e.g. a ramp)
        ids = []         # element_ids applied to the keys
        unique_ids = []  # list of each element_id used in order 
        vals = []        # values form the data assigned to the key tiles
        rots = []        # rotation of each key tile
        subsets = self.tiled_map.groupby("element_id")
        for id, geom, rot in zip(elements.element_id,
                                 elements.geometry,
                                 elements.rotation):
            subset = subsets.get_group(id)
            d = subset[self.variables[id]]
            # if the data are categorical then it's complicated...
            if d.dtype == pd.CategoricalDtype:
                # desired order of categorical variable is the 
                # color maps dictionary keys
                num_cats = len(self.colourmaps[self.variables[id]])
                val_order = dict(zip(
                    self.colourmaps[self.variables[id]].keys(),
                    range(num_cats)))
                # compile counts of each category
                coded_data_counts = [0] * num_cats
                for v in list(d):
                    coded_data_counts[val_order[v]] += 1
                # make list of the categories containing appropriate 
                # counts of each in the order needed using a reverse lookup
                order_val = list(val_order.keys())
                data_vals = []
                for i, n in enumerate(coded_data_counts):
                    data_vals.extend([order_val[i]] * n)
            else: # any other data is easy!
                data_vals = sorted(d)
            vals.extend(data_vals)
            n = len(data_vals)
            key = self.tiling.tile_unit._get_legend_key_shapes(geom, n, rot)
            key_tiles.extend(key)
            ids.extend([id] * n)
            unique_ids.append(id)
            rots.extend([rot] * n)
        # finally make up a data table with all the data in all the
        # columns (each set of data only gets used in the subset it
        # applies to). This allows us to reuse the tiling_utils.
        # plot_subsetted_gdf() function
        key_data = {}
        for id in unique_ids:
            key_data[self.variables[id]] = vals
        
        key_gdf = gpd.GeoDataFrame(
            data = key_data | {"element_id": ids, "rotation": rots}, 
            crs = self.tiled_map.crs,
            geometry = gpd.GeoSeries(key_tiles))
        key_gdf.geometry = key_gdf.rotate(self.tiling.rotation, origin = (0, 0))
        return key_gdf
        
        
    def explore(self) -> None:
        return None
