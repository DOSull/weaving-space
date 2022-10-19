#!/usr/bin/env python
# coding: utf-8

"""Classes for tiling maps. `weavingspace.tile_map.Tiling` and 
`weavingspace.tile_map.TiledMap` are exposed in the  public API and 
respectively enable creation of a tiling and plotting of the tiling as a 
multivariate map. 
"""

from dataclasses import dataclass
from typing import Union
import itertools
import copy

import numpy as np
import geopandas as gpd
import pandas as pd

from matplotlib.figure import Figure
import matplotlib.colors
import matplotlib.pyplot as pyplot

import shapely.geometry as geom
import shapely.affinity as affine

from weavingspace.tileable import Tileable
from weavingspace.tileable import TileShape
from weavingspace.tile_unit import TileUnit

import weavingspace.tiling_utils as tiling_utils

from time import perf_counter

@dataclass
class _TileGrid():
    """A class to represent the translation centres of a tiling.
    
    We store the grid as a GeoSeries of Point objects to make it
    simple to plot it in map views if required.
    
    Implementation relies on transforming the translation vectors into a square
    space where the tile spacing is unit squares, then transforming this back 
    into the original map space. Some member variables of the class are in 
    the transformed grid generation space, some in the map space.

    Args:
        tile (TileUnit): the base tile in map space.
        to_tile (gpd.GeoSeries): geometry of region to be tiled in map space.
        transform (tuple[float]): the forward transformation from map space to
            the grid generation space. Stored in shapely's (a, b, d, e, dx, dy)
            format.
        inverse_transform (tuple[float]): the inverse transformation from
            grid generation space to map space. Stored in shapely's (a, b, d, 
            e, dx, dy) format.
        centre (tuple[float]): centre point of the extent in map space.
        points (gpd.GeoSeries): shapely.geometry.Points recording translation 
            vectors of the tiling in map space.
        _extent (gpd.GeoSeries): geometry of the (circular) extent of the   
            tiling transformed to grid generation space.
        _at_centroids (bool): if True the grid will consist of the centroids of
            the spatial units in the to_tile region, allowing a simple way to
            use a tile unit as a point symbol.
    """
    tile:TileUnit = None
    to_tile:gpd.GeoSeries = None
    transform:tuple[float] = None
    inverse_transform:tuple[float] = None
    centre:tuple[float] = None
    points:gpd.GeoSeries = None
    _extent:gpd.GeoSeries = None
    _at_centroids:bool = False
    
    def __init__(self, tile:TileUnit, to_tile:gpd.GeoSeries, 
                 at_centroids:bool = False):
        self.tile = tile
        # self.to_tile = tiling_utils.clean_polygon(
        #     gpd.GeoSeries([to_tile.unary_union]))
        self.to_tile = self._get_area_to_tile(to_tile)
        self.inverse_transform, self.transform = self.get_transforms()
        self.extent, self.centre = self._get_extent()
        self._at_centroids = at_centroids
        if self._at_centroids:
            self.points = to_tile.centroid
        else:
            self.points = self.get_grid()
        self.points.crs = self.tile.crs
    
    
    def _get_area_to_tile(self, to_tile) -> geom.Polygon:
        bb = to_tile.total_bounds
        poly = geom.Polygon(((bb[0], bb[1]), 
                             (bb[2], bb[1]), 
                             (bb[2], bb[3]), 
                             (bb[0], bb[3])))
        return gpd.GeoSeries([poly])


    def _get_extent(self) -> tuple[gpd.GeoSeries, geom.Point]:
        """Returns the extent and centre of the grid.
        
        Extent is in the grid-generation space.

        Returns:
            tuple[gpd.GeoSeries, geom.Point]: the extent of the grid and its 
                centre.
        """
        mrr = self.to_tile[0].minimum_rotated_rectangle
        mrr_centre = geom.Point(mrr.centroid.coords[0])
        mrr_corner = geom.Point(mrr.exterior.coords[0])
        radius = mrr_centre.distance(mrr_corner)
        ext = affine.affine_transform(mrr_centre.buffer(radius), 
                                      self.transform)
        return gpd.GeoSeries([ext]), (mrr_centre.x, mrr_centre.y)


    def get_transforms(self) -> tuple[float]:
        """Returns the forward and inverse transforms from map space to
        grid generation space.
        
        In grid generation space the translation vectors are (1, 0) and (0, 1)
        so we can simply form a matrix from two translation vectors and invert
        it to get the forward transform. The inverse transform is the matrix
        formed from the vectors.

        Returns:
            tuple[float]: shapely affine transform tuples (a, b, d, e, dx, dy).
                See https://shapely.readthedocs.io/en/stable/manual.html#affine-transformations for details.
        """
        v = self.tile.get_vectors()
        vector_array = np.array([[v[0][0], v[1][0]], 
                                 [v[0][1], v[1][1]]])
        inv_tfm = np.linalg.inv(vector_array)
        return (self._np_to_shapely_transform(vector_array), 
                self._np_to_shapely_transform(inv_tfm))


    def get_grid(self) -> gpd.GeoSeries:
        """Generates the grid transformed into map space
        
        Obtain dimensions of the transformed region, then set down a uniform 
        grid.
        
        Returns:
            gpd.GeoSeries: the grid as a collection of geom.Points.
        """
        _w, _h, _l, _b = tiling_utils.get_width_height_left_bottom(self.extent)
        w = int(np.ceil(_w))
        h = int(np.ceil(_h))
        l = _l - (w - _w) / 2
        b = _b - (h - _h) / 2
        mesh = np.array(np.meshgrid(np.arange(w) + l,
                                    np.arange(h) + b)).reshape((2, w * h)).T
        pts = [geom.Point(p[0], p[1]) for p in mesh]
        return gpd.GeoSeries([p for p in pts if p.within(self.extent[0])]) \
            .affine_transform(self.inverse_transform)


    def _np_to_shapely_transform(self, mat:np.ndarray) -> tuple[float]:
        """Converts a numpy affine transform matrix to shapely format.

        Args:
            mat (np.ndarray): numpy array to convert.

        Returns:
            tuple[float]: shapely affine transform tuple
        """
        return (mat[0][0], mat[0][1], mat[1][0], mat[1][1], 0, 0)


@dataclass
class Tiling:
    """Class that applies a `Tileable` object to a region to be mapped.
    
    The result of the tiling procedure is stored in the `tiles` variable and
    covers a region sufficient that the tiling can be rotated to any desired 
    angle. 
    """
    tile_unit:Tileable = None
    tile_shape:TileShape = None
    region:gpd.GeoDataFrame = None
    grid:_TileGrid = None
    tiles:gpd.GeoDataFrame = None
    rotation:float = 0.0

    def __init__(self, unit:Tileable, region:gpd.GeoDataFrame, id_var = None,
                 tile_margin:float = 0, elements_sf:float = 1, 
                 elements_margin:float = 0, as_icons:bool = False) -> None:
        """Class to persist a tiling by filling an area relative to 
        a region sufficient to apply the tiling at any rotation.
        
        The Tiling constructor allows a number of adjustments to the supplied
        `weavingspace.tileable.Tileable` object:
        
        + `tile_margin` values greater than 0 will introduce spacing of  
        the specified distance between elements on the boundary of each tile
        by applying the `TileUnit.inset_tile()` method. Note that this 
        operation does not make sense for `WeaveUnit` objects,
        and may not preserve the equality of tile element areas.
        + `elements_sf` values less than one scale down tile elements by 
        applying the `TileUnit.scale_elements()` method. Does not make sense 
        for `WeaveUnit` objects.
        + `elements_margin` values greater than one apply a negative buffer of
        the specified distance to every element in the tile by applying the
        `Tileable.inset_elements()` method. This option is applicable to both 
        `WeaveUnit` and `TileUnit` objects.

        Args:
            unit (Tileable): the tile_unit to use.
            region (gpd.GeoDataFrame): the region to be tiled.
            tile_margin (float, optional): values greater than one apply an 
                inset margin to the tile unit. Defaults to 0.
            elements_sf (float, optional): scales the tile unit elements. 
                Defaults to 1.
            elements_margin (float, optional): applies a negative buffer to 
                the tile unit elements. Defaults to 0.
            as_icons (bool, optional): if True tiles will only be placed at
                the region's zone centroids, one per zone. Defaults to 
                False. 
        """
        self.tile_unit = unit
        self.rotation = self.tile_unit.rotation
        if elements_margin > 0:
            self.tile_unit = self.tile_unit.inset_elements(elements_margin) 
        if elements_sf != 1:
            if isinstance(self.tile_unit, TileUnit):
                self.tile_unit = self.tile_unit.scale_elements(elements_sf)
            else:
                print(f"""Applying scaling to elements of a WeaveUnit does not make sense. Ignoring elements_sf setting of {elements_sf}.""")
        if tile_margin > 0:
            if isinstance(self.tile_unit, TileUnit):
                self.tile_unit = self.tile_unit.inset_tile(tile_margin)
            else:
                print(f"""Applying a tile margin to elements of a WeaveUnit does not make sense. Ignoring tile_margin setting of {tile_margin}.""")        
        self.region = region
        self.region.sindex
        if id_var != None:
            print("""id_var is no longer required and will be deprecated soon. A temporary unique index attribute is added and removed when generating the tiled map.""")
        if as_icons:
            self.grid = _TileGrid(self.tile_unit, self.region.geometry, True)
        else:
            self.grid = _TileGrid(self.tile_unit, self.region.geometry)
        self.tiles = self.make_tiling()
        self.tiles.sindex


    def get_tiled_map(self, rotation:float = 0.,
                      prioritise_tiles:bool = True,
                      tiles_or_elements:str = "elements",
                      ragged_edges:bool = True, 
                      use_centroid_lookup_approximation = False,
                      debug = False) -> "TiledMap":
        """Returns a `TiledMap` filling a region at the requested rotation.
        
        HERE BE DRAGONS! This function took a lot of trial and error to get 
        right, so modify with CAUTION!
        
        The `proritise_tiles = True` option means that the tiling will not 
        break up the elements in `TileUnit`s at the boundaries between areas 
        in the mapped region, but will instead ensure that tile elements remain
        complete, picking up their data from the region zone which they overlap
        the most. The `tiles_or_elements` flag (experimental!) potentially 
        allows prioritisation at the level of tile unit, not element. Possibly
        just a 'hook' for now to be checked later.
        
        The exact order in which operations are performed affects performance. 
        For example, the final clipping to self.region when ragged_edges = False
        is _much_ slower if it is carried out before the dissolving of tile 
        elements into the region zones. So... again... modify CAREFULLY!
        
        Args:
            rotation (float, optional): An optional rotation to apply. Defaults 
                to 0.
            prioritise_tiles (bool, optional): if True tiles will not be 
                broken at boundaries in the region dataset. Defaults to True.
            tiles_or_elements (str, optional): if "tiles" then the tile 
                prioritisation applies to tiles; if "elements" then it applies
                to tile unit elements. Defaults to "elements".
            ragged_edges (bool, optional): if True tiles at the edge of the 
                region will not be cut by the region extent - ignored if 
                prioritise_tiles is False when edges will always be clipped to 
                the region extent. Defaults to True.
            use_centroid_lookup_approximation (bool, optional): if True use
                element centroids for lookup of region data - ignored if 
                prioritise_tiles is False when it is irrelevant. Defaults to
                False.
            debug (bool, optional): if True prints timing messages. Defaults
                to False.

        Returns:
            TiledMap: a TiledMap of the source region.
        """
        if debug:
            t1 = perf_counter()
        
        id_var = self._setup_region_DZID()
        tiled_map = self.rotated(rotation)
        # compile a list of the variable names we are NOT going to change
        # i.e. everything except the geometry and the id_var
        region_vars = list(self.region.columns)
        region_vars.remove("geometry")
        region_vars.remove(id_var)

        if debug:
            t2 = perf_counter()
            print(f"STEP 1: prep data (rotation if requested): {t2 - t1:.3f}")
            
        if prioritise_tiles:  # maintain tile continuity across zone boundaries
            # select only tiles inside a spacing buffer of the region
            # make column with unique ID for every element in the tiling
            if tiles_or_elements == "tiles":
                # the join ID is based on the tile
                tiled_map["joinUID"] = self.tiles.tile_id
            else:
                # the join ID is unique per element
                tiled_map["joinUID"] = list(range(tiled_map.shape[0]))

            if use_centroid_lookup_approximation:
                t5 = perf_counter()
                tile_pts = copy.deepcopy(tiled_map)
                tile_pts.geometry = tile_pts.centroid
                lookup = tile_pts.sjoin(
                    self.region, how = "inner")[["joinUID", id_var]]
            else:
                # determine areas of overlapping tile elements and drop the data
                # we join the data back later, so dropping makes that easier
                # overlaying in region.overlay(tiles) seems to be faster??
                # TODO: also... this part is performance-critical, think about # fixes -- possibly including the above centroid-based approx
                overlaps = self.region.overlay(tiled_map, make_valid = False)
                if debug:
                    t3 = perf_counter()
                    print(f"STEP A2: overlay zones with tiling: {t3 - t2:.3f}")
                overlaps["area"] = overlaps.geometry.area
                if debug:
                    t4 = perf_counter()
                    print(f"STEP A3: calculate areas: {t4 - t3:.3f}")
                overlaps.drop(columns = region_vars, inplace = True)
                if debug:
                    t5 = perf_counter()
                    print(f"STEP A4: drop columns prior to join: {t5 - t4:.3f}")
                # make a lookup by largest area element to region id
                lookup = overlaps \
                    .iloc[overlaps.groupby("joinUID")["area"] \
                    .agg(pd.Series.idxmax)][["joinUID", id_var]]
            # now join the lookup and from there the region data
            if debug:
                t6 = perf_counter()
                print(f"STEP A5: build lookup for join: {t6 - t5:.3f}")
            tiled_map = tiled_map \
                .merge(lookup, on = "joinUID") \
                .merge(self.region.drop(columns = ["geometry"]), on = id_var) 
            if debug:
                t7 = perf_counter()
                print(f"STEP A6: perform lookup join: {t7 - t6:.3f}")

        else:  # here we overlay
            tiled_map = self.region.overlay(tiled_map)
            t7 = perf_counter()
            if debug:
                print(f"STEP B2: overlay tiling with zones: {t7 - t2:.3f}")
        
        tiled_map.drop(columns = [id_var, "joinUID"], inplace = True)
        self.region.drop(columns = [id_var], inplace = True)
        
        # if we've retained tiles and want 'clean' edges, then clip
        # note that this step is slow: geopandas unary_unions the clip layer
        if prioritise_tiles and not ragged_edges:
            tiled_map.sindex
            tiled_map = tiled_map.clip(self.region)
            if debug:
                print(f"""STEP A7/B3: clip map to region: 
                      {perf_counter() - t7:.3f}""")

        tm = TiledMap()
        tm.tiling = self
        tm.map = tiled_map
        return tm


    def _setup_region_DZID(self) -> str:
        """Creates a new guaranteed-unique attribute in the self.region 
        dataframe, and returns its name.
        
        Avoids a name clash with any existing attribute in the dataframe.

        Returns:
            str: name of the added attribute.
        """
        dzid = "DZID"
        i = 0
        while dzid in self.region.columns: 
            dzid = "DZID" + str(i)
            i = i + 1
        self.region[dzid] = list(range(self.region.shape[0]))
        return dzid
    
    
    def _rotate_gdf_to_geoseries(
            self, gdf:gpd.GeoDataFrame, 
            angle:float, centre:tuple = (0, 0)
        ) -> tuple[gpd.GeoSeries, tuple[float]]:
        """Rotates the geometries in a GeoDataFrame as a single collection.
        
        Rotation is about the supplied centre or about the centroid of the 
        GeoDataFrame (if not). This allows for reversal of  a rotation. [Note 
        that this might not be a required precaution!]

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
        tiles = itertools.chain(*[
            self.tile_unit.elements.geometry.translate(p.x, p.y)
            for p in self.grid.points])
        # replicate the element ids
        ids = list(self.tile_unit.elements.element_id) * len(self.grid.points)
        tile_ids = sorted(list(range(len(self.grid.points))) * 
                          self.tile_unit.elements.shape[0])
        tiles_gs = gpd.GeoSeries(tiles)
        # assemble and return as a GeoDataFrame
        tiles_gdf = gpd.GeoDataFrame(data = {"element_id": ids,
                                             "tile_id": tile_ids},
                                     geometry = tiles_gs, 
                                     crs = self.tile_unit.crs)
        # unclear if we need the below or not...
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
            data = {"element_id": self.tiles.element_id,
                    "tile_id": self.tiles.tile_id}, 
            crs = self.tiles.crs,
            geometry = self.tiles.geometry.rotate(rotation, 
                                                  origin = self.grid.centre))


@dataclass
class TiledMap:
    """Class representing a tiled map. Should not be accessed directly, but 
    will be created by calling `Tiling.get_tiled_map()`. After creation the 
    variables and colourmaps attributes can be set, and then 
    `TiledMap.render()` called to make a map. Settable attributes are explained 
    in documentation of the `TiledMap.render()` method. 
    
    Examples:
        Recommended usage is as follows. First, make a `TiledMap` from a `Tiling` object.
        
            tm = tiling.get_tiled_map(...)
            
        Some options in the `Tiling` constructor affect the map appearance. See 
        `Tiling` for details.
        
        Once a `TiledMap` object exists, set options on it, either when calling 
        `TiledMap.render()` or explicitly, i.e.
        
            tm.render(opt1 = val1, opt2 = val2, ...)
            
        or
        
            tm.opt1 = val1
            tm.opt2 = val2
            tm.render()
        
        Option settings are persistent, i.e. unless a new `TiledMap` object is 
        created the option settings have to be explicitly reset to default 
        values on subsequent calls to `TiledMap.render()`.
        
        The most important options are the `variables` and `colourmaps` 
        settings.
        
        `variables` is a dictionary mapping `weavingspace.tileable.Tileable` 
        element_ids (usually "a", "b", etc.) to variable names in the data. For 
        example, 
        
            tm.variables = dict(zip(["a", "b"], ["population", "income"]))
            
        `colourmaps` is a dictionary mapping dataset variable names to the 
        matplotlib colourmap to be used for each. For example,
        
            tm.colourmaps = dict(zip(tm.variables.values(), ["Reds", "Blues"]))
        
        See [this notebook](https://github.com/DOSull/weaving-space/blob/main/weavingspace/example-tiles-cairo.ipynb) for simple usage. TODO: This 
        more complicated example shows how categorical maps can be created.
    """
    # these will be set at instantion by Tiling.get_tiled_map()
    tiling:Tiling = None  #: the Tiling with the required tiles
    map:gpd.GeoDataFrame = None  #: the GeoDataFrame on which this map is based 
    
    # variables and colourmaps should be set before calling self.render()
    variables:dict[str,str] = None  #: lookup from element_id to variable names
    colourmaps:dict[str,Union[str,dict]] = None  #: lookup from variables to cmaps
    
    # the below parameters can be set either before calling self.render() 
    # or passed in as parameters to self.render()
    # these are solely TiledMap.render() options
    legend:bool = True  #: whether or not to show a legend
    legend_zoom:float = 1.0  #: <1 zooms out from legend to show more context
    legend_dx:float = 0.  #: horizontal shift of legend relative to the map
    legend_dy:float = 0.  #: vertical shift of legend relative to the map
    use_ellipse:bool = False  #: if True clips legend with an ellipse
    ellipse_magnification:float = 1.0  #: magnification to apply to clip ellipse
    radial_key:bool = False  #: if True use radial key even for ordinal/ratio data (normally these will be shown by concentric elements)
    draft_mode:bool = False  #: if True plot only the map coloured by element_id
    
    # the parameters below are geopandas.plot options which we intercept to
    # ensure they are applied appropriately when we plot a GDF
    scheme:str = "equalinterval"  #: geopandas scheme to apply
    k:int = 100  #: geopandas number of classes to apply
    figsize:tuple[float] = (20, 15)  #: maptlotlib figsize 
    dpi:float = 72  #: dots per inch for bitmap formats
        

    def render(self, **kwargs) -> Figure:
        """Renders the current state to a map.
        
        Note that TiledMap objects will usually be created by calling 
        `Tiling.get_tiled_map()`.
        
        Args:
            variables (dict[str,str]): Mapping from element_id values to
                variable names. Defaults to None.
            colourmaps (dict[str,Union[str,dict]]): Mapping from variable
                names to colour map, either a colour palette as used by
                geopandas/matplotlib, a fixed colour, or a dictionary mapping
                categorical data values to colours. Defaults to None.
            legend (bool): If True a legend will be drawn. Defaults to True.
            legend_zoom (float): Zoom factor to apply to the legend. Values <1 
                will show more of the tile context. Defaults to 1.0.
            legend_dx (float): x shift to apply to the legend position.
                Defaults to 0.0.
            legend_dy (float): x and y shift to apply to the legend position. 
                Defaults to 0.0.
            use_ellipse (bool): If True applies an elliptical clip to the   
                legend. Defaults to False.
            ellipse_magnification (float): Magnification to apply to ellipse
                clipped legend. Defaults to 1.0.
            radial_key (bool): If True legend key for TileUnit maps will be
                based on radially dissecting the tiles. Defaults to False.
            draft_mode (bool): If True a map of the tiled map coloured by
                tile element_ids (and with no legend) is returned. Defaults
                to False.
            scheme (str): passed to geopandas.plot for numeric data. Defaults to
                "equalinterval".
            k (int): passed to geopandas.plot for numeric data. Defaults to 100.
            figsize (tuple[float,floar]): plot dimensions passed to geopandas.
                plot. Defaults to (20,15).
            dpi (float): passed to pyplot.plot. Defaults to 72.
            **kwargs: other settings to pass to pyplot/geopandas.plot. 
            
        Returns:
            matplotlib.figure.Figure: figure on which map is plotted.
        """
        pyplot.rcParams['pdf.fonttype'] = 42
        pyplot.rcParams['pdf.use14corefonts'] = True
        matplotlib.rcParams['pdf.fonttype'] = 42

        to_remove = set()  # keep track of kwargs we use to setup TiledMap
        for k, v in kwargs.items():
            if k in self.__dict__:
                self.__dict__[k] = v
                to_remove.add(k)
        # remove them so we don't pass them on to pyplot and get errors
        for k in to_remove:
            del kwargs[k]
            
        if self.draft_mode:
            fig = pyplot.figure(figsize = self.figsize) 
            ax = fig.add_subplot(111)
            self.map.plot(ax = ax, column = "element_id", cmap = "tab20", 
                          **kwargs)
            return fig

        if self.legend:
            # this sizing stuff is rough and ready for now, possibly forever... 
            reg_w, reg_h, *_ = \
                tiling_utils.get_width_height_left_bottom(self.map.geometry)
            tile_w, tile_h, *_ = \
                tiling_utils.get_width_height_left_bottom(
                    self.tiling.tile_unit._get_legend_elements().rotate(
                        self.tiling.rotation, origin = (0, 0)))
            sf_w, sf_h = reg_w / tile_w / 3, reg_h / tile_h / 3
            gskw = {"height_ratios": [sf_h * tile_h, reg_h - sf_h * tile_h],
                    "width_ratios": [reg_w, sf_w * tile_w]}

            fig, axes = pyplot.subplot_mosaic(
                [["map", "legend"],
                 ["map", "."]], 
                gridspec_kw = gskw, figsize = self.figsize, 
                layout = "constrained", **kwargs)
        else:
            fig, axes = pyplot.subplots(
                1, 1, figsize = self.figsize, 
                layout = "constrained", **kwargs)
            
        if self.variables is None:
            # get any floating point columns available
            default_columns = \
                self.map.select_dtypes(
                    include = ("float64", "int64")).columns
            self.variables = dict(zip(
                self.map.element_id.unique(), 
                list(default_columns)))
            print(f"""No variables specified, picked the first
                  {len(self.variables)} numeric ones available.""")        
        elif isinstance(self.variables, (list, tuple)):
            self.variables = dict(zip(
                self.tiling.tile_unit.elements.element_id.unique(),
                self.variables))
            print(f"""Only a list of variables specified, assigning to 
                  available element_ids.""")
                    
        if self.colourmaps is None:
            self.colourmaps = {}
            for var in self.variables.values():
                if self.map[var].dtype == pd.CategoricalDtype:
                    self.colourmaps[var] = "tab20"
                    print(f"""For categorical data, you should specify colour 
                          mapping explicitly.""")
                else:
                    self.colourmaps[var] = "Reds"
        
        self._plot_map(axes, **kwargs)
        return fig
    
    
    def _plot_map(self, axes:pyplot.Axes, **kwargs) -> None:
        """Plots map to the supplied axes.

        Args:
            axes (pyplot.Axes): axes on which maps will be drawn.
        """
        bb = self.map.geometry.total_bounds
        if self.legend:
            axes["map"].set_axis_off()
            axes["map"].set_xlim(bb[0], bb[2])
            axes["map"].set_ylim(bb[1], bb[3])
            self._plot_subsetted_gdf(axes["map"], self.map, **kwargs)
            self.plot_legend(ax = axes["legend"], **kwargs)
            if (self.legend_dx != 0 or self.legend_dx != 0):
                box = axes["legend"].get_position()
                box.x0 = box.x0 + self.legend_dx
                box.x1 = box.x1 + self.legend_dx
                box.y0 = box.y0 + self.legend_dy
                box.y1 = box.y1 + self.legend_dy
                axes["legend"].set_position(box)
        else:
            axes.set_axis_off()
            axes.set_xlim(bb[0], bb[2])
            axes.set_ylim(bb[1], bb[3])
            self._plot_subsetted_gdf(axes, self.map, **kwargs)
        return None
    
    
    def _plot_subsetted_gdf(self, ax:pyplot.Axes, 
                           gdf:gpd.GeoDataFrame, **kwargs) -> None:
        """Plots a gpd.GeoDataFrame multiple times based on a subsetting
        attribute (assumed to be "element_id").
        
        NOTE: used to plot both the main map _and_ the legend.

        Args:
            ax (pyplot.Axes): axes to plot to.
            gdf (gpd.GeoDataFrame): the GeoDataFrame to plot.

        Raises:
            Exception: if self.colourmaps cannot be parsed exception is raised.
        """
        groups = gdf.groupby("element_id")
        for id, var in self.variables.items():
            subset = groups.get_group(id)
            # Handle custom color assignments via 'cmaps' parameter.
            # Result is setting 'cmap' variable used in plot command afterwards.
            if (isinstance(self.colourmaps[var], dict)):
                colormap_dict = self.colourmaps[var]
                data_unique_sorted = sorted(subset[var].unique())
                cmap = matplotlib.colors.ListedColormap(
                    [colormap_dict[x] for x in data_unique_sorted])
                subset.plot(ax = ax, column = var, cmap = cmap, **kwargs)
            else:
                if (isinstance(self.colourmaps, 
                                (str, matplotlib.colors.Colormap,
                                matplotlib.colors.LinearSegmentedColormap,
                                matplotlib.colors.ListedColormap))):
                    cmap = self.colourmaps   # one palette for all ids
                elif (len(self.colourmaps) == 0):
                    cmap = 'Reds'  # set a default... here, to Brewer's 'Reds'
                elif (var not in self.colourmaps):
                    cmap = 'Reds'  # no color specified in dict, use default
                elif (isinstance(self.colourmaps[var],
                                (str, matplotlib.colors.Colormap,
                                matplotlib.colors.LinearSegmentedColormap,
                                matplotlib.colors.ListedColormap))):
                    cmap = self.colourmaps[var]  # specified colors for this var
                else:
                    raise Exception(f"""Color map for '{var}' is not a known 
                                    type, but is {str(type(self.colourmaps[var])
                                    )}""")

                subset.plot(ax = ax, column = var, cmap = cmap, 
                            scheme = self.scheme, k = self.k, **kwargs)
    
    
    def to_file(self, fname:str = None) -> None:
        """Outputs the tiled map to a layered GPKG file. 
        
        Currently delegates to `weavingspace.tiling_utils.write_map_to_layers()`.

        Args:
            fname (str, optional): Filename to write. Defaults to None.
        """
        tiling_utils.write_map_to_layers(self.map, fname)
        return None

    
    def plot_legend(self, ax: pyplot.Axes = None, **kwargs) -> None:
        """Plots a legend for this tiled map.

        Args:
            ax (pyplot.Axes, optional): axes to draw legend. Defaults to None.
        """
        # turn off axes (which seems also to make it impossible
        # to set a background colour)
        ax.set_axis_off()

        legend_elements = self.tiling.tile_unit._get_legend_elements()
        # this is a bit hacky, but we will apply the rotation to text 
        # annotation so for TileUnits which don't need it, reverse that now
        if isinstance(self.tiling.tile_unit, TileUnit):
            legend_elements.rotation = -self.tiling.rotation
        
        legend_key = self._get_legend_key_gdf(legend_elements)
                
        legend_elements.geometry = legend_elements.geometry.rotate(
            self.tiling.rotation, origin = (0, 0))
        
        if self.use_ellipse:
            ellipse = tiling_utils.get_bounding_ellipse(
                legend_elements.geometry, 
                mag = self.ellipse_magnification)
            bb = ellipse.total_bounds
            c = ellipse.unary_union.centroid
        else:
            bb = legend_elements.geometry.total_bounds
            c = legend_elements.geometry.unary_union.centroid

        # apply legend zoom - NOTE that this must be applied even
        # if self.legend_zoom is not == 1...
        ax.set_xlim(c.x + (bb[0] - c.x) / self.legend_zoom, 
                    c.x + (bb[2] - c.x) / self.legend_zoom)
        ax.set_ylim(c.y + (bb[1] - c.y) / self.legend_zoom, 
                    c.y + (bb[3] - c.y) / self.legend_zoom)

        # plot the legend key elements (which include the data)
        self._plot_subsetted_gdf(ax, legend_key, lw = 0, **kwargs)
        
        for id, tile, rotn in zip(self.variables.keys(),
                                  legend_elements.geometry,
                                  legend_elements.rotation):
            c = tile.centroid
            ax.annotate(self.variables[id], xy = (c.x, c.y), 
                    ha = "center", va = "center", rotation_mode = "anchor", 
                    # adjust rotation to favour text reading left to right
                    rotation = (rotn + self.tiling.rotation + 90) % 180 - 90, 
                    bbox = {"lw": 0, "fc": "#ffffff40"})

        # now plot background; we include the central tiles, since in
        # the weave case these may not match the legend elements
        context_tiles = self.tiling.tile_unit.get_local_patch(r = 2, 
            # include_0 = isinstance(self.tiling.tile_unit, WeaveUnit)) \
            include_0 = True) \
                .geometry.rotate(self.tiling.rotation, origin = (0, 0))
        # for reasons escaping all reason... invalid polygons sometimes show up 
        # here I think because of the rotation /shrug... in any case, this 
        # sledgehammer should fix it
        context_tiles = gpd.GeoSeries([g.simplify(1e-6)
                                       for g in context_tiles.geometry],
                                      crs = self.tiling.tile_unit.crs)
        if self.use_ellipse:
            context_tiles.clip(ellipse, keep_geom_type = False).plot(
                    ax = ax, fc = "#9F9F9F3F", lw = 0.0)
            tiling_utils.get_boundaries(context_tiles.geometry).clip(
                ellipse, keep_geom_type = True).plot(
                    ax = ax, ec = "#5F5F5F", lw = 1)
        else:
            context_tiles.plot(ax = ax, fc = "#9F9F9F3F", 
                               ec = "#5F5F5F", lw = 0.0)
            tiling_utils.get_boundaries(context_tiles.geometry).plot(
                ax = ax, ec = "#5F5F5F", lw = 1)


    def _get_legend_key_gdf(
            self, elements:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """Returns a GeoDataFrame of tile elements dissected and with
        data assigned to the slice so a map of them can stand as a legend.
        
        'Dissection' is handled differently by `WeaveUnit` and `TileUnit` 
        objects and delegated to either `WeaveUnit._get_legend_key_shapes()` or `TileUnit._get_legend_key_shapes()`.

        Args:
            elements (gpd.GeoDataFrame): the legend elements.

        Returns:
            gpd.GeoDataFrame:  with element_id, variables and rotation
                attributes, and geometries of Tileable elements sliced into a 
                colour ramp or set of nested tiles.
        """
        key_tiles = []   # set of tiles to form a colour key (e.g. a ramp)
        ids = []         # element_ids applied to the keys
        unique_ids = []  # list of each element_id used in order 
        vals = []        # the data assigned to the key tiles
        rots = []        # rotation of each key tile
        subsets = self.map.groupby("element_id")
        for (id, var), geom, rot in zip(self.variables.items(),
                                 elements.geometry,
                                 elements.rotation):
            subset = subsets.get_group(id)
            d = subset[var]
            radial = False
            # if the data are categorical then it's complicated...
            if d.dtype == pd.CategoricalDtype:
                radial = True and self.radial_key
                # desired order of categorical variable is the 
                # color maps dictionary keys
                cmap = self.colourmaps[var]
                num_cats = len(cmap)
                val_order = dict(zip(cmap.keys(), range(num_cats)))
                # compile counts of each category
                freqs = [0] * num_cats
                for v in list(d):
                    freqs[val_order[v]] += 1
                # make list of the categories containing appropriate 
                # counts of each in the order needed using a reverse lookup
                data_vals = list(val_order.keys())
                data_vals = [data_vals[i] for i, f in enumerate(freqs) if f > 0]
            else: # any other data is easy!
                data_vals = sorted(d)
                freqs = [1] * len(data_vals)
            key = self.tiling.tile_unit._get_legend_key_shapes(
                geom, freqs, rot, radial)
            key_tiles.extend(key)
            vals.extend(data_vals)
            n = len(data_vals)
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
            crs = self.map.crs,
            geometry = gpd.GeoSeries(key_tiles))
        key_gdf.geometry = key_gdf.rotate(self.tiling.rotation, origin = (0, 0))
        return key_gdf
        
        
    def explore(self) -> None:
        """TODO: add wrapper to make tiled web map via geopandas.explore.
        """
        return None
