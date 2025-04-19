import marimo

__generated_with = "0.12.10"
app = marimo.App(
    width="full",
    app_title="MapWeaver",
    layout_file="layouts/tiling-explorer-2.grid.json",
    css_file="",
    html_head_file="splash.html",
)


@app.cell(hide_code=True)
def _(centred, mo):
    mo.vstack([
        mo.image(src="mw.png").style(centred),
        mo.md(f"<span title='Weaving maps of complex data'>2025.04.20-09:00</span>").style({'background-color':'rgba(255,255,255,0.5'}).center(),
        mo.md(f"<span title='Requires weavingspace 0.0.6.42'>0.0.6.42</span>").style({'background-color':'rgba(255,255,255,0.5','font-style':'italic'}).center(),
    ])
    return


@app.cell(hide_code=True)
def module_imports():
    import io                     # for in-memory loading of inputs and outputs
    from pathlib import Path
    from typing import Union
    import matplotlib as mpl
    import math
    import pandas as pd
    import geopandas as gpd
    from shapely import is_valid
    import weavingspace as wsp
    from weavingspace import CMAPS_SEQUENTIAL
    from weavingspace import CMAPS_DIVERGING
    return (
        CMAPS_DIVERGING,
        CMAPS_SEQUENTIAL,
        Path,
        Union,
        gpd,
        io,
        is_valid,
        math,
        mpl,
        pd,
        wsp,
    )


@app.cell(hide_code=True)
def globals(CMAPS_DIVERGING, CMAPS_SEQUENTIAL, get_colour_ramp, gpd):
    ok_message = "STATUS All good!"
    dummy_data_file = "https://raw.githubusercontent.com/DOSull/weaving-space/refs/heads/main/examples/data/dummy-data.json"
    builtin_gdf = gpd.read_file(dummy_data_file, engine="fiona")
    available_palettes = [pal for pal in CMAPS_SEQUENTIAL + CMAPS_DIVERGING if pal[-2:] != "_r" and 
                         pal not in ["berlin", "managua", "vanimo"]] # these don't seem to be in MPL3.8.4
    # make a bunch of colour ramps and save them in a dictionary so we only
    # have to make them at notebook initialisation
    color_ramps = {k: get_colour_ramp(k) 
                   for k in available_palettes + [p + "_r" for p in available_palettes]}
    return (
        available_palettes,
        builtin_gdf,
        color_ramps,
        dummy_data_file,
        ok_message,
    )


@app.cell(hide_code=True)
def marimo_states(
    available_palettes,
    builtin_gdf,
    dummy_data_file,
    mo,
    ok_message,
):
    # mo.state() variables so that some of app attributes persist in use
    # the input data source
    get_input_data, set_input_data = mo.state(dummy_data_file)
    # the GeoDataframe we are tiling - using state means we can restore
    # current GDF if a new one has issues 
    get_gdf, set_gdf = mo.state(builtin_gdf)
    # variables currently selected by the user
    get_variables, set_variables = mo.state([])
    # palettes currently selected by the user
    get_palettes, set_palettes = mo.state(available_palettes)
    # palette reverse settings
    get_reversed, set_reversed = mo.state([False] * 18)
    # app status message
    get_status_message, set_status_message = mo.state(ok_message)
    return (
        get_gdf,
        get_input_data,
        get_palettes,
        get_reversed,
        get_status_message,
        get_variables,
        set_gdf,
        set_input_data,
        set_palettes,
        set_reversed,
        set_status_message,
        set_variables,
    )


@app.cell(hide_code=True)
def upload_data(mo, set_input_data, tool_tip):
    _file_browser = mo.ui.file(filetypes=[".geojson", ".json"], 
                    on_change=set_input_data, label=f"Upload")
    mo.md(f"{tool_tip(_file_browser, "Your data should be polygons. Currently only GeoJSON formatted data is readable.")}").left()
    return


@app.cell(hide_code=True)
def read_gdf(
    builtin_gdf,
    get_gdf,
    get_input_data,
    get_numeric_variables,
    gpd,
    io,
    is_valid,
    ok_message,
    set_gdf,
    set_status_message,
    set_variables,
    wsp,
):
    # flags to control recovery if load fails for some reason
    _done = False
    _did_repairs = False
    # if input data is a string then it's the initial setting
    # otherwise it's the upload button result unless that is length 0
    if isinstance(get_input_data(), str) or len(get_input_data()) == 0:
        set_gdf(builtin_gdf)
    else:
        # store current GDF in case we need to restore it
        _old_gdf = get_gdf()
        try:
            # after much experimentation this seems to be the most reliable 
            # way to load JSON from bytes contents of the upload button
            _new_gdf = gpd.read_file(io.BytesIO(get_input_data()[0].contents).read().decode())
        except Exception as e:
            # would be better to figure out what exceptions can happen and trap
            # the exact type, but we'll keep an eye on that and update
            # meanwhile... reset to the current GDF and flag error
            set_gdf(_old_gdf)
            set_status_message("ERROR! Exception in uploading data")
            print(e.args)
            raise e
        else:
            # if not enough variables for tiling then flag and reset to current GDF
            _n = len(get_numeric_variables(_new_gdf))
            if _n < 2:
                set_status_message("WARNING! One or fewer variables, data not loaded")
                set_gdf(_old_gdf)
                _done = True # we're done
            if not _done:
                # check no topology issues and attempt repair
                if not all(is_valid(_new_gdf.geometry)):
                    _new_gdf.geometry = wsp.tiling_utils.repair_polygon(_new_gdf.geometry)
                    # recheck and if repair failed report error
                    if not all(is_valid(_new_gdf.geometry)):
                        set_status_message("ERROR! Geometries not valid, try another dataset")
                        set_gdf(_old_gdf)
                        _done = True
                    else:
                        _did_repairs = True
            if not _done:
                # check data are projected and reproject if not
                if not _new_gdf.crs.is_projected:
                    _new_gdf = _new_gdf.to_crs(3857)
                set_variables(get_numeric_variables(_new_gdf))
                if _did_repairs:
                    set_status_message("WARNING! Repaired geometry errors in your data")
                else:
                    set_status_message(ok_message)
                set_gdf(_new_gdf)
    return


@app.cell(hide_code=True)
def _(mo, tool_tip):
    download_type = mo.ui.dropdown(options=["GeoJSON", "GeoPackage", "SVG", "PNG"], value="GeoJSON")
    mo.md(f"{tool_tip(download_type, 'Set the file format for downloaded map data')}")
    return (download_type,)


@app.cell(hide_code=True)
def _(
    Path,
    download_type,
    get_input_data,
    io,
    mo,
    result,
    tiled_map,
    tiling_map,
    tool_tip,
):
    mo.stop(tiling_map)
    if isinstance(get_input_data(), str) or len(get_input_data()) == 0:
        _fname = f"{Path(get_input_data()).with_suffix('').stem}-map-weaver-map"
    else:
        _fname = f"{Path(get_input_data()[0].name).stem}-map-weaver-map"

    if download_type.value == "GeoJSON":
        _download_button = mo.download(data=tiled_map.map.to_json().encode('utf-8'), 
                                       filename=f'{_fname}.geojson', 
                                       mimetype='text/plain', 
                                       label='Download')
    elif download_type.value == "GeoPackage":
        # here we have to use a BytesIO stream to write to
        with io.BytesIO() as _f:
            tiled_map.map.to_file(_f, driver="GPKG", engine="fiona")
            _download_button = mo.download(data=_f, 
                                           filename=f'{_fname}.gpkg', 
                                           mimetype="application/geopackage+sqlite3",
                                           label='Download')
    elif download_type.value == "SVG":
        with io.BytesIO() as _f:
            result.savefig(_f, format="svg")
            _download_button = mo.download(data=_f,
                                           filename=f'{_fname}.svg',
                                           mimetype="image/svg+xml",
                                           label="Download")
    else:
        with io.BytesIO() as _f:
            result.savefig(_f, format="png", dpi=300)
            _download_button = mo.download(data=_f,
                                           filename=f'{_fname}.png',
                                           mimetype="image/png",
                                           label="Download")

    mo.md(f'{tool_tip(_download_button, "Download the tiled map")}')
    return


@app.cell
def _(get_gdf, get_modded_tile_unit, wsp):
    tiling_map = True # flag to block some other cells (potentially...)
    tiled_map = wsp.Tiling(get_modded_tile_unit(), get_gdf()).get_tiled_map(join_on_prototiles=False)
    tiling_map = False # stop blocking other cells
    return tiled_map, tiling_map


@app.cell(hide_code=True)
def _(get_selected_colour_palettes, get_tile_ids, get_variables, tiled_map):
    tiled_map.ids_to_map = [id for id, v in zip(get_tile_ids(), get_variables()) if v != "---"]
    tiled_map.vars_to_map = [v for v in get_variables() if v != "---"]
    tiled_map.colors_to_use = [c for c, v in zip(get_selected_colour_palettes(), get_variables()) if v != "---"]
    result = tiled_map.render(legend=False)
    result
    return (result,)


@app.cell(hide_code=True)
def set_number_of_variables(mo, tool_tip):
    num_tiles = mo.ui.slider(steps=range(2, 13), 
                             # steps=[x for x in range(2, 16)] + [18], # if we figure out why >12 is a problem for the app...
                             value=4, 
                             debounce=True, 
                             show_value=True)
    mo.md(f"""
    ### General settings
    #### Set number of tiling elements {tool_tip(num_tiles, 'Choose the number of distinct tiles you want to use to symbolise data.')}
    """)
    return (num_tiles,)


@app.cell(hide_code=True)
def variable_palette_map_header(mo):
    mo.md(f"### Variables &lrarr; palettes")
    return


@app.cell(hide_code=True)
def status_panel(
    get_gdf,
    get_numeric_variables,
    get_status_message,
    mo,
    num_tiles,
    set_status_message,
):
    set_status_message("STATUS All good!")
    if len(get_numeric_variables(get_gdf())) < num_tiles.value:
        set_status_message(f"WARNING! More tiles ({num_tiles.value}) than variables ({len(get_numeric_variables(get_gdf()))})")

    if "ERROR" in get_status_message():
        _bkgd, _col = "red", "white"
    elif "WARNING" in get_status_message():
        _bkgd, _col = "pink", "black"
    else: # OK
        _bkgd, _col = "lightgreen", "black"

    mo.md(f"<span style='background-color:{_bkgd};color:{_col};font-face:sans-serif;padding:2px;'>{get_status_message()}</span>").center()
    return


@app.cell(hide_code=True)
def build_variable_and_palette_dropdowns(
    available_palettes,
    get_gdf,
    get_numeric_variables,
    get_palettes,
    get_reversed,
    get_variables,
    mo,
    num_tiles,
    set_palettes,
    set_reversed,
    set_variables,
    tiling_map,
):
    mo.stop(tiling_map)
    if num_tiles.value < len(get_variables()):
        set_variables(get_variables()[:num_tiles.value])
    elif num_tiles.value > len(get_variables()):
        # add as many additional variables as we have to work with
        n_to_add = min(num_tiles.value, len(get_numeric_variables(get_gdf()))) - len(get_variables())
        to_add = [v for v in get_numeric_variables(get_gdf()) if v not in get_variables()][:n_to_add]
        set_variables(get_variables() + to_add)

    _chosen_palettes = get_palettes()[:num_tiles.value]
    _chosen_reversed = get_reversed()[:num_tiles.value]
    variables = mo.ui.array(
        [mo.ui.dropdown(options=get_numeric_variables(get_gdf()) + ["---"], value=v) for v in get_variables()], 
        on_change=set_variables) 
    pals = mo.ui.array(
        [mo.ui.dropdown(options=available_palettes, value=p) for p in _chosen_palettes], 
        on_change=set_palettes)
    rev_pals = mo.ui.array(
        [mo.ui.switch(r) for r in _chosen_reversed], on_change=set_reversed)
    return n_to_add, pals, rev_pals, to_add, variables


@app.cell(hide_code=True)
def build_var_palette_mapping(
    color_ramps,
    get_tile_ids,
    mo,
    pals,
    rev_pals,
    tool_tip,
    variables,
):
    _cols = [pal.value + ("_r" if rev.value else "") for pal, rev in zip(pals, rev_pals)]
    mo.md("\n".join([
        "&nbsp;&nbsp;".join([f"#### Tiles `{tile_id}`",
                             f"{tool_tip(var, f"Variable for tiles with id {tile_id}")} &rarr;",
                             f"{tool_tip(pal, f"Palette for variable {var.value}")}",
                             f"<span style='position:relative;top:5px;'>{tool_tip(rev, 'Reverse ramp')}</span>", 
                             f"<span style='display:inline-block;object-fit:cover;height:24px;position:relative;bottom:22px;'>{color_ramps[col]}</span>",
    ]) for tile_id, var, pal, rev, col in zip(get_tile_ids(), variables, pals, rev_pals, _cols)]))
    return


@app.cell(hide_code=True)
def set_spacing_limits(get_spacings, mo):
    _spacing_steps, _spacing_value = get_spacings()
    spacing = mo.ui.slider(steps=_spacing_steps, value=_spacing_value,
                           show_value=True, debounce=True)
    return (spacing,)


@app.cell
def _(mo, tile_or_weave):
    mo.md(f"### {tile_or_weave.value.capitalize()} modifiers")
    return


@app.cell(hide_code=True)
def tiling_modifier_ui_elements(mo):
    tile_rotate = mo.ui.slider(steps=range(-90, 91, 1), 
                               value=0, 
                               show_value=True, 
                               debounce=True)
    tile_scale_x = mo.ui.slider(start=1, 
                                stop=4, 
                                step=0.1, 
                                value=1, 
                                show_value=True, 
                                debounce=True)
    tile_scale_y = mo.ui.slider(start=1, 
                                stop=4, 
                                step=0.1, 
                                value=1, 
                                show_value=True, 
                                debounce=True)
    tile_skew_x = mo.ui.slider(steps=range(-40, 41, 1),
                               value=0, 
                               show_value=True, 
                               debounce=True)
    tile_skew_y = mo.ui.slider(steps=range(-40, 41, 1),
                               value=0, 
                               show_value=True, 
                               debounce=True)
    p_inset = mo.ui.slider(start=0, 
                           stop=10, 
                           step = 0.1, 
                           value=0, 
                           show_value=True, 
                           debounce=True)
    t_inset = mo.ui.slider(start=0, 
                           stop=5, 
                           step = 0.1, 
                           value=0, 
                           show_value=True, 
                           debounce=True)
    return (
        p_inset,
        t_inset,
        tile_rotate,
        tile_scale_x,
        tile_scale_y,
        tile_skew_x,
        tile_skew_y,
    )


@app.cell(hide_code=True)
def setup_tiling_modifiers(
    mo,
    p_inset,
    spacing,
    t_inset,
    tile_or_weave,
    tile_rotate,
    tile_scale_x,
    tile_scale_y,
    tile_skew_x,
    tile_skew_y,
    tiling_map,
    tool_tip,
):
    mo.stop(tiling_map)
    if tile_or_weave.value == "tiling":
        _str = f"""
            #### Spacing {tool_tip(spacing, 'In units of the map CRS, the approximate dimension of the repeating group.')}
            #### Rotate by {tool_tip(tile_rotate, 'Rotate tiling (degrees). Note that the tile group is rotated before any scaling or skew transforms are applied.')}
            #### Scale left-right {tool_tip(tile_scale_x, 'Scale in the x direction.')}
            #### Scale up-down {tool_tip(tile_scale_y, 'Scale in the y direction.')}
            #### Skew left-right {tool_tip(tile_skew_x, 'Skew in the x direction (degrees).')}
            #### Skew up-down {tool_tip(tile_skew_y, 'Skew in the y direction (degrees).')}
            #### Tile group inset {tool_tip(p_inset, 'Inset the tile group (% spacing).')}
            #### Tiles inset {tool_tip(t_inset, 'Inset individual tiles (% spacing).')}
            """
    else:
        _str = f"""
            #### Spacing {tool_tip(spacing, 'In units of the map CRS, the distance between strand centre lines.')}
            #### Rotate by {tool_tip(tile_rotate, "'Rotate weave (degrees). Note that the weave is rotated before any scaling or skew transforms are applied.'")}
            #### Scale left-right {tool_tip(tile_scale_x, "'Scale in the x direction.'")}
            #### Scale up-down {tool_tip(tile_scale_y, "'Scale in the y direction.'")}
            #### Skew left-right {tool_tip(tile_skew_x, "'Skew in the x direction (degrees).'")}
            #### Skew up-down {tool_tip(tile_skew_y, "'Skew in the y direction (degrees).'")}
            #### Strands inset {tool_tip(t_inset, "'Inset strands (% width).'")}
            """
    mo.md(_str)
    return


@app.cell(hide_code=True)
def tiling_or_weave_chooser(mo, num_tiles, tilings_by_n, tool_tip):
    _options = list(set([v["type"] for v in tilings_by_n[num_tiles.value].values()]))
    tile_or_weave = mo.ui.dropdown(options=_options, value="tiling", label="#### Pick tiling or weave")
    mo.md(f"{tool_tip(tile_or_weave, 'Choose tiling or a weave tiling')}")
    return (tile_or_weave,)


@app.cell(hide_code=True)
def tiling_type_chooser(mo, num_tiles, tile_or_weave, tilings_by_n, tool_tip):
    _options = [k for k, v in tilings_by_n[num_tiles.value].items() if v["type"] == tile_or_weave.value]

    family = mo.ui.dropdown(options=_options, 
                            value=_options[0], 
                            label=f"#### {tile_or_weave.value.capitalize()} type")
    mo.md(f"{tool_tip(family, "Choose tiling family")}")
    return (family,)


@app.cell(hide_code=True)
def tiling_design_plot(get_modded_tile_unit, plot_tiles):
    plot_tiles(get_modded_tile_unit())
    return


@app.cell(hide_code=True)
def setup_chosen_tiling_options(
    family,
    mo,
    num_tiles,
    tile_or_weave,
    tilings_by_n,
):
    if "slice" in family.value:
        _offset = mo.ui.slider(steps=[x / 100 for x in range(101)], 
                               value=0,
                               label="#### Offset",
                               show_value=True, 
                               debounce=True) 
    elif "dissect" in family.value:
        _offset = mo.ui.number(start=0, 
                               stop=1,
                               value=0,
                               label="#### Offset", 
                               debounce=True)
        if "hex" in family.value:
            _offset_angle = mo.ui.slider(steps=[x for x in range(-50, 86)],
                                         value=0,
                                         label="#### Inner angle",
                                         show_value=True,
                                         debounce=True)
        else:
            _offset_angle = mo.ui.slider(steps=[x for x in range(-30, 71)],
                                         value=0,
                                         label="#### Inner angle",
                                         show_value=True,
                                         debounce=True)
    elif "star1" in family.value:
        _point_angle = mo.ui.slider(steps=[x for x in range(10, 121)],
                                    value=30,
                                    label="#### Point angle",
                                    show_value=True,
                                    debounce=True)
    elif "weave" in family.value:
        _aspect = mo.ui.slider(steps=[x / 12 for x in range(1,13)], 
                               value=3/4, 
                               label="#### Strand width",
                               show_value=True, 
                               debounce=True)
        if not "cube" in family.value:
            _over_under = mo.ui.text(value=tilings_by_n[num_tiles.value][family.value]["n"],
                                     label="#### Over-under pattern")

    if tile_or_weave.value == "tiling":
        if "slice" in family.value:
            tile_spec = mo.ui.dictionary({"offset": _offset,})
            tooltips = ["0 starts at the base tile corners, 1 at the mid-point along segments equally dividing the base tile perimeter into the requested number of tiles."]
        elif "dissect" in family.value:
            tile_spec = mo.ui.dictionary({"offset": _offset, "offset_angle": _offset_angle,})
            tooltips = ["0 starts at the base tile corners, 1 at the mid-point along segments equally dividing the base tile perimeter into the requested number of tiles.", "Angle by which inner polygon is rotated relative to the outer. Reverse the angle for a different look."]
        elif "star1" in family.value:
            tile_spec = mo.ui.dictionary({"point_angle": _point_angle})
            tooltips = ["Set the interior angle of the star points."]
        else:
            tile_spec = None
    else:
        if "plain" in family.value or "cube" in family.value:
            tile_spec = mo.ui.dictionary({"aspect": _aspect})
            tooltips = ["The width of the weave strands relative to strand spacing. A value of 1 will fill the map with no gaps. Progressively smaller values will leave 'holes' in the woven pattern and will eventually give the appearance of a cross hatch."]
        else:
            tile_spec = mo.ui.dictionary({"aspect": _aspect, "over_under": _over_under,})
            tooltips = [
                "The width of the weave strands relative to strand spacing.", 
                "Comma-separated sequence of how many strands in the other direction each ribbon should go over, then under. Use values that are factors of the numbers of strands in each direction to avoid very large and complex repeating units."
            ]
    return tile_spec, tooltips


@app.cell(hide_code=True)
def additional_tiling_options(
    mo,
    tile_or_weave,
    tile_spec,
    tiling_map,
    tool_tip,
    tooltips,
):
    mo.stop(tiling_map)
    if tile_spec is None:
        _show_options = mo.md(f"#### No {tile_or_weave.value} options to set")
    else:
        _show_options = mo.md("\n".join(
            [f"### Set {tile_or_weave.value} options"] + 
            [f"#### {tool_tip(v, tt)}" for (k, v), tt in zip(tile_spec.items(), tooltips)]
        ))
    _show_options
    return


@app.cell(hide_code=True)
def _(mo):
    _radius = mo.ui.slider(steps=range(5), 
                           value=0, 
                           show_value=True)
    _show_prototile = mo.ui.switch(value=False)
    _show_reg_prototile = mo.ui.switch(value=False)
    _show_vectors = mo.ui.switch(value=False)
    _show_ids = mo.ui.switch(value=True)
    _show_scale = mo.ui.switch(value=False)

    view_settings = mo.ui.dictionary({
        "radius": _radius,
        "show_ids": _show_ids,
        "show_prototile": _show_prototile,
        "show_reg_prototile": _show_reg_prototile,
        "show_vectors": _show_vectors,
        "show_scale": _show_scale,
    })
    return (view_settings,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(f"### Design view options")
    return


@app.cell(hide_code=True)
def design_view_ui_elements(mo, tool_tip, view_settings):
    mo.md(f"""
    #### {tool_tip(view_settings['show_ids'], 'Show the tiling element labels used to match tiles to variables in the map data.')} Show tile IDs
    #### {tool_tip(view_settings['show_prototile'], 'Show in fine black outline the simple tile (usually a square or hexagon) which forms the basis of the pattern.')} Show base tile
    #### {tool_tip(view_settings['show_vectors'], 'Show the translations that map repeating tiles in the pattern onto one another.')} Show vectors
    #### {tool_tip(view_settings['show_reg_prototile'], 'Show in a red outline the repeating set tile group that pieces together jigsaw-like to form the pattern.')} Show &lsquo;jigsaw piece&rsquo;
    #### {tool_tip(view_settings['show_scale'], 'Give an indication of scale in map units.')} Show scale
    #### Tile group 'shells' to show {tool_tip(view_settings['radius'], 'The number of &lsquo;shells&rsquo; of the tiling to show around the base tile group.')}
    """)
    return


@app.cell(hide_code=True)
def _(
    family,
    get_gdf,
    get_over_under,
    get_tile_ids,
    num_tiles,
    spacing,
    tile_or_weave,
    tile_spec,
    tilings_by_n,
    wsp,
):
    def get_base_tile_unit() -> wsp.Tileable:
        """
        Returns a Tileable object based on the current configuration of the app UI
        """
        # call TileUnit or WeaveUnit with appropriate settings
        spec = tilings_by_n[num_tiles.value][family.value]
        if tile_or_weave.value == "tiling":
            result = wsp.TileUnit(
                tiling_type=spec["tiling_type"],
                spacing=spacing.value,
                n=spec["n"] if "n" in spec else None,
                code=spec["code"] if "code" in spec else None,
                offset=tile_spec["offset"].value if "offset" in spec else None,
                offset_angle=tile_spec["offset_angle"].value if "offset_angle" in spec else None,
                point_angle=tile_spec["point_angle"].value if "point_angle" in spec else None,
                crs=get_gdf().crs)
        else:
            result = wsp.WeaveUnit(
                weave_type=spec["weave_type"],
                spacing=spacing.value,
                strands=spec["strands"],
                n=get_over_under(tile_spec["over_under"].value) \
                    if spec["weave_type"] in ["twill", "basket"] else 1,
                aspect=tile_spec["aspect"].value,
                crs=get_gdf().crs)
        # remove tiles for which no variables are available
        slice = [id in get_tile_ids() for id in result.tiles.tile_id]
        result.tiles = result.tiles.loc[slice]
        return result
    return (get_base_tile_unit,)


@app.cell(hide_code=True)
def _(
    get_base_tile_unit,
    p_inset,
    spacing,
    t_inset,
    tile_or_weave,
    tile_rotate,
    tile_scale_x,
    tile_scale_y,
    tile_skew_x,
    tile_skew_y,
    tile_spec,
):
    def get_modded_tile_unit():
        """
        Returns a Tileable object based on the current configuration of the app UI
        including any modifications (rotations, skews, scaling etc.)
        """
        # make any requested modifications to the tile or weave unit
        if tile_or_weave.value == "tiling":
            return get_base_tile_unit() \
               .transform_rotate(tile_rotate.value) \
               .transform_scale(tile_scale_x.value, tile_scale_y.value) \
               .transform_skew(tile_skew_x.value, tile_skew_y.value) \
               .inset_tiles(t_inset.value * spacing.value / 100) \
               .inset_prototile(p_inset.value * spacing.value / 100)
        else:
            return get_base_tile_unit() \
               .transform_rotate(tile_rotate.value) \
               .transform_scale(tile_scale_x.value, tile_scale_y.value) \
               .transform_skew(tile_skew_x.value, tile_skew_y.value) \
               .inset_tiles(t_inset.value * tile_spec["aspect"].value * spacing.value / 100)
    return (get_modded_tile_unit,)


@app.cell(hide_code=True)
def _(Union):
    def get_over_under(pattern:str) -> Union[int,tuple[int]]:
        """
        Returns either an integer or tuple of integers of even length 
        based on the supplied comma-separated string of integer values
        """
        # convert string over under pattern to tuple[int]
        # if any invalid characters in string return a useful default
        if any([not c in "0123456789," for c in pattern]): return (2, 2)
        numbers = [int(s) for s in pattern.split(",")]
        # has to be an even number of elements so trim if needed
        length = 2 * len(numbers) // 2
        if length == 0:
            return (2, 2)
        else:
            return tuple(numbers[:length])
    return (get_over_under,)


@app.cell(hide_code=True)
def _(
    get_selected_colour_palettes,
    get_tile_ids,
    get_variables,
    mpl,
    num_tiles,
    view_settings,
    wsp,
):
    def plot_tiles(tiles:wsp.Tileable) -> mpl.axis:
        """
        Returns a matplotlib axis with the Tileable supplied plotted on it
        """
        cols = [mpl.colormaps.get(p)(2/3) for p in get_selected_colour_palettes()]
        while len(cols) < num_tiles.value:
            cols = cols + ["#ffffff00"]
        cm = mpl.colors.ListedColormap(cols)
        plot = tiles.plot(r=view_settings["radius"].value, 
                          show_vectors=view_settings["show_vectors"].value, 
                          show_prototile=view_settings["show_prototile"].value,
                          show_reg_prototile=view_settings["show_reg_prototile"].value,
                          show_ids=view_settings["show_ids"].value,
                          cmap=cm, figsize=(4.5, 4), r_alpha=0.5)
        if view_settings["show_scale"].value:
            plot.xaxis.set_visible(True)
            plot.yaxis.set_visible(False)
            plot.set_frame_on(False)
        else:
            plot.axis("off")
        _blanks = [id for id, v in zip(get_tile_ids(), get_variables()) if v == "---"]
        if len(_blanks) > 0:
            tiles.tiles.loc[tiles.tiles["tile_id"].isin(_blanks)].plot(ax = plot.axes, fc="#00000000", hatch="////", ec="w")
        return plot
    return (plot_tiles,)


@app.cell(hide_code=True)
def _(get_gdf, get_numeric_variables):
    def get_tile_ids() -> list[str]:
        """
        Returns a list of lower case letters starting ['a', 'b', ...] length limited
        by the number of available numeric variables in the data
        """
        return list("abcdefghijkl")[:len(get_numeric_variables(get_gdf()))]
    return (get_tile_ids,)


@app.cell(hide_code=True)
def _(num_tiles, pals, rev_pals):
    def get_selected_colour_palettes() -> list[str]:
        """
        Returns a list of the currently selected colourmap names, including
        _r suffix where the reverse switch has been set
        """
        return [(p if not r else p + "_r") 
                for p, r in zip(pals.value[:num_tiles.value], 
                                rev_pals.value[:num_tiles.value])]
    return (get_selected_colour_palettes,)


@app.cell(hide_code=True)
def _(io, mo, mpl):
    def get_colour_ramp(pal_name:str="Reds", rev:bool=False) -> mo.image:
        """
        Returns a marimo image to represent the colour ramp specified by the supplied
        base name (`pal`) and whether or not it is reversed (the `rev`) setting when
        _r will be appended to the base name
        """
        # returns an image representing the colour ramp for the supplied palette name
        fig, ax = mpl.pyplot.subplots()
        # this code based on code in matplotlib docs at
        # https://matplotlib.org/stable/users/explain/colors/colormaps.html
        xy = [[x / 256 for x in range(257)] for i in range(2)]
        ax.imshow(xy, aspect=30, cmap=mpl.colormaps.get(pal_name + ("_r" if rev else "")))
        ax.set_axis_off()
        # write to memory and embed in a marimo image object
        with io.BytesIO() as buf:
            mpl.pyplot.savefig(buf, dpi=24, pad_inches=0, bbox_inches="tight")
            mpl.pyplot.close(fig)
            img = mo.image(buf)
        return img
    return (get_colour_ramp,)


@app.cell(hide_code=True)
def _(gpd, pd):
    def get_numeric_variables(_gdf:gpd.GeoDataFrame) -> list[str]:
        """
        Returns a list of the column names in the supplied GeoDataFrame which
        are numeric
        """
        return [col for col in _gdf.columns if not "geom" in col 
                and pd.api.types.is_numeric_dtype(_gdf[col].dtype)]
    return (get_numeric_variables,)


@app.cell(hide_code=True)
def _(get_gdf, math):
    def get_spacings() -> list[int]:
        """
        Returns a list of int values to apply as steps for the spacing UI element
        """
        # returns reasonable rounded set of spacing options given the bounds of the data
        _bb = get_gdf().total_bounds
        _width, _height = _bb[2] - _bb[0], _bb[3] - _bb[1]
        _max = 10 ** math.floor(math.log10(max(_width, _height))) // 5
        _mid = _max // 2
        _min = _max // 20
        _stepsize1 = _min // 10
        _stepsize2 = _min // 2
        _steps = [x for x in range(_min, _mid, _stepsize1)] + \
                 [x for x in range(_mid, _max + 1, _stepsize2)]
        # if tile_or_weave.value == "weave":
        #     _steps = [x // 2 for x in _steps]
        #     _mid = _mid // 2
        return _steps, 3 * _mid // 4
    return (get_spacings,)


@app.cell(hide_code=True)
def tool_tip():
    def tool_tip(ele:str, tip:str) -> str:
        """
        Returns a HTML <span> string putting a tooltip around the supplied element 
        """
        # convenience function to add a tooltip to supplied string
        return f'<span title="{tip}">{ele}</span>'
    return (tool_tip,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        f"""
        ### Available tilings
        Organised by number of variables to symbolise
        """
    )
    return


@app.cell(hide_code=True)
def setup_tilings_dictionary():
    tilings_by_n = {
      2: {
        "plain weave a|b": dict(type="weave", weave_type="plain", strands="a|b", n="1"),
        "twill weave a|b": dict(type="weave", weave_type="twill", strands="a|b", n="2"),
        "twill weave a|b-": dict(type="weave", weave_type="twill", strands="a|b-", n="2"),
        "twill weave a|b 3": dict(type="weave", weave_type="twill", strands="a|b", n="3"),
        "twill weave a|b- 3": dict(type="weave", weave_type="twill", strands="a|b-", n="3"),
        "twill weave a|b 4": dict(type="weave", weave_type="twill", strands="a|b", n="4"),
        "twill weave a|b 1,2": dict(type="weave", weave_type="twill", strands="a|b", n="1,2"),
        "twill weave a|b- 1,2": dict(type="weave", weave_type="twill", strands="a|b-", n="1,2"),
        "twill weave a|b 1,2,2,1": dict(type="weave", weave_type="twill", strands="a|b", n="1,2,2,1"),
        "twill weave a|b 2,3": dict(type="weave", weave_type="twill", strands="a|b", n="2,3"),
        "twill weave a|b 2,3,3,2": dict(type="weave", weave_type="twill", strands="a|b", n="2,3,3,2"),
        "archimedean 4.8.8": dict(type="tiling", tiling_type="archi", code="4.8.8"),
        "square-slice 2": dict(type="tiling", tiling_type="square-slice", n=2, offset=0),
        "crosses 2": dict(type="tiling", tiling_type="cross", n=2),
        "star1 44": dict(type="tiling", tiling_type="star1", code="44", point_angle=30),
        "hex-colouring 2": dict(type="tiling", tiling_type="hex-col", n=2),
        "square-colouring 2": dict(type="tiling", tiling_type="square-col", n=2),
        "hex-slice 2": dict(type="tiling", tiling_type="hex-slice", n=2, offset=0),
        "laves 3.3.3.4.4": dict(type="tiling", tiling_type="lave", code="3.3.3.4.4"),
      },
      3: {
        "cube weave a--|b--|c--": dict(type= "weave", weave_type="cube", strands="a--|b--|c--"),
        "plain weave ab|c": dict(type="weave", weave_type="plain", strands="ab|c", n="1"),
        "plain weave ab-|c": dict(type="weave", weave_type="plain", strands="ab-|c", n="1"),
        "twill weave ab|c-": dict(type="weave", weave_type="twill", strands="ab|c-", n="2"),
        "basket weave ab|c-": dict(type="weave", weave_type="basket", strands="ab|c-", n="2"),
        "hex-slice 3": dict(type="tiling", tiling_type="hex-slice", n=3, offset=0),
        "hex-colouring 3": dict(type="tiling", tiling_type="hex-col", n=3),
        "crosses 3": dict(type="tiling", tiling_type="cross", n=3),
        "star1 33": dict(type="tiling", tiling_type="star1", code="33", point_angle=30),
        "square-colouring 3": dict(type="tiling", tiling_type="square-col", n=3),
        "archimedean 3.6.3.6": dict(type="tiling", tiling_type="archi", code="3.6.3.6"),
        "archimedean 3.12.12": dict(type="tiling", tiling_type="archi", code="3.12.12"),
        "square-slice 3": dict(type="tiling", tiling_type="square-slice", n=3, offset=0),
        "archimedean 3.3.3.4.4": dict(type="tiling", tiling_type="archi", code="3.3.3.4.4"),
        "hex-dissection 3": dict(type="tiling", tiling_type="hex-dissect", n=3, offset=0, offset_angle=0),
        "square-dissection 3": dict(type="tiling", tiling_type="square-dissect", n=3, offset=0, offset_angle=0),
        "star2 66": dict(type="tiling", tiling_type="star2", code="66"),
        "star1 36": dict(type="tiling", tiling_type="star1", code="36", point_angle=30),
        "star1 63": dict(type="tiling", tiling_type="star1", code="63", point_angle=30),
      },
      4: {
        "laves 3.3.4.3.4": dict(type="tiling", tiling_type="laves", code="3.3.4.3.4"),
        "plain weave ab|cd": dict(type="weave", weave_type="plain", strands="ab|cd", n="1"),
        "plain weave ab-|cd": dict(type="weave", weave_type="plain", strands="ab-|cd", n="1"),
        "plain weave ab-|cd-": dict(type="weave", weave_type="plain", strands="ab-|cd-", n="1"),
        "basket weave ab|cd": dict(type="weave", weave_type="basket", strands="ab|cd", n="2"),
        "basket weave ab|cd 3": dict(type="weave", weave_type="basket", strands="ab|cd", n="3"),
        "twill weave ab|cd": dict(type="weave", weave_type="twill", strands="ab|cd", n="2"),
        "twill weave ab|cd 3": dict(type="weave", weave_type="twill", strands="ab|cd", n="3"),
        "twill weave ab|cd 1,2": dict(type="weave", weave_type="twill", strands="ab|cd", n="1,2"),
        "twill weave ab|cd 1,2,2,1": dict(type="weave", weave_type="twill", strands="ab|cd", n="1,2,2,1"),
        "basket weave ab-|cd": dict(type="weave", weave_type="basket", strands="ab-|cd", n="3"),
        "twill weave ab-|cd": dict(type="weave", weave_type="twill", strands="ab-|cd", n="3"),
        "basket weave ab-|cd-": dict(type="weave", weave_type="basket", strands="ab-|cd-", n="3"),
        "twill weave ab-|cd-": dict(type="weave", weave_type="twill", strands="ab-|cd-", n="3"),
        "crosses 4": dict(type="tiling", tiling_type="cross", n=4),
        "square-slice 4": dict(type="tiling", tiling_type="square-slice", n=4, offset=0),
        "square-colouring 4": dict(type="tiling", tiling_type="square-col", n=4),
        "hex-colouring 4": dict(type="tiling", tiling_type="hex-col", n=4),
        "star2 64": dict(type="tiling", tiling_type="star2", code="64"),
        "star2 464": dict(type="tiling", tiling_type="star2", code="464"),
        "star2 844": dict(type="tiling", tiling_type="star2", code="844"),
        "hex-slice 4": dict(type="tiling", tiling_type="hex-slice", n=4, offset=0),
        "hex-dissection 4": dict(type="tiling", tiling_type="hex-dissect", n=4, offset=0, offset_angle=0),
      },
      5: {
        "square-colouring 5": dict(type="tiling", tiling_type="square-col", n=5),
        "crosses 5": dict(type="tiling", tiling_type="crosses", n=5),
        "plain weave abc|de": dict(type="weave", weave_type="plain", strands="abc|de", n="1"),
        "plain weave abc-|de": dict(type="weave", weave_type="plain", strands="abc-|de", n="1"),
        "plain weave abc-|de-": dict(type="weave", weave_type="plain", strands="abc-|de-", n="1"),
        "twill weave abc|de": dict(type="weave", weave_type="twill", strands="abc|de", n="3"),
        "twill weave abc|de-": dict(type="weave", weave_type="twill", strands="abc|de-", n="3"),
        "twill weave abc-|de-": dict(type="weave", weave_type="twill", strands="abc-|de-", n="3"),
        "basket weave abc|de": dict(type="weave", weave_type="basket", strands="abc|de", n="3"),
        "basket weave abc|de-": dict(type="weave", weave_type="basket", strands="abc|de-", n="3"),
        "basket weave abc-|de-": dict(type="weave", weave_type="basket", strands="abc-|de-", n="3"),
        "hex-colouring 5": dict(type="tiling", tiling_type="hex-col", n=5),
        "hex-slice 5": dict(type="tiling", tiling_type="hex-slice", n=5, offset=0),
        "square-slice 5": dict(type="tiling", tiling_type="square-slice", n=5, offset=0),
        "square-dissection 5": dict(type="tiling", tiling_type="square-dissect", n=5, offset=0, offset_angle=0),
      },
      6: {
        "hex-slice 6": dict(type="tiling", tiling_type="hex-slice", n=6, offset=0),
        "square-slice 6": dict(type="tiling", tiling_type="square-slice", n=6, offset=0),
        "square-colouring 6": dict(type="tiling", tiling_type="square-col", n=6),
        "laves 3.3.3.3.6": dict(type="tiling", tiling_type="lave", code="3.3.3.3.6"),
        "laves 3.12.12": dict(type="tiling", tiling_type="lave", code="3.12.12"),
        "cube weave a-b|c-d|e-f": dict(type= "weave", weave_type="cube", strands="a-b|c-d|e-f"),
        "plain weave abc|def": dict(type="weave", weave_type="plain", strands="abc|def", n="1"),
        "plain weave abc-|def": dict(type="weave", weave_type="plain", strands="abc-|def", n="1"),
        "plain weave abc-|def-": dict(type="weave", weave_type="plain", strands="abc-|def-", n="1"),
        "basket weave abc|def": dict(type="weave", weave_type="basket", strands="abc|def", n="3"),
        "twill weave abc|def": dict(type="weave", weave_type="twill", strands="abc|def", n="3"),
        "basket weave abc-|def": dict(type="weave", weave_type="basket", strands="abc-|def", n="3"),
        "twill weave abc-|def": dict(type="weave", weave_type="twill", strands="abc-|def", n="3"),
        "basket weave abc-|def-": dict(type="weave", weave_type="basket", strands="abc-|def-", n="4"),
        "twill weave abc-|def-": dict(type="weave", weave_type="twill", strands="abc-|def-", n="4"),
        "archimedean 3.3.4.3.4": dict(type="tiling", tiling_type="archi", code="3.3.4.3.4"),
        "archimedean 3.4.6.4": dict(type="tiling", tiling_type="archi", code="3.4.6.4"),
        "archimedean 4.6.12": dict(type="tiling", tiling_type="archi", code="4.6.12"),
        "crosses 6": dict(type="tiling", tiling_type="cross", n=6),
        "star2 45": dict(type="tiling", tiling_type="star2", code="45"),
        "star2 663": dict(type="tiling", tiling_type="star2", code="663"),
        "hex-colouring 6": dict(type="tiling", tiling_type="hex-col", n=6),
      },
      7: {
        "hex-colouring 7": dict(type="tiling", tiling_type="hex-col", n=7),
        "crosses 7": dict(type="tiling", tiling_type="cross", n=7),
        "plain weave abcd|efg": dict(type="weave", weave_type="plain", strands="abcd|efg", n="1"),
        "plain weave abcd|efg-": dict(type="weave", weave_type="plain", strands="abcd|efg-", n="1"),
        "plain weave abcd-|efg-": dict(type="weave", weave_type="plain", strands="abcd-|efg-", n="1"),
        "basket weave abcd|efg": dict(type="weave", weave_type="basket", strands="abcd|efg", n="3"),
        "basket weave abcd|efg-": dict(type="weave", weave_type="basket", strands="abcd|efg-", n="4"),
        "basket weave abcd|efg- 2": dict(type="weave", weave_type="basket", strands="abcd|efg-", n="2"),
        "twill weave abcd|efg": dict(type="weave", weave_type="twill", strands="abcd|efg", n="3"),
        "twill weave abcd|efg-": dict(type="weave", weave_type="twill", strands="abcd|efg-", n="4"),
        "twill weave abcd|efg- 2": dict(type="weave", weave_type="twill", strands="abcd|efg-", n="2"),
        "square-colouring 7": dict(type="tiling", tiling_type="square-col", n=7),
        "hex-slice 7": dict(type="tiling", tiling_type="hex-slice", n=7, offset=0),
        "square-slice 7": dict(type="tiling", tiling_type="square-slice", n=7, offset=0),
        "hex-dissection 7": dict(type="tiling", tiling_type="hex-dissect", n=7, offset=0, offset_angle=0),
      },
      8: {
        "square-slice 8": dict(type="tiling", tiling_type="square-slice", n=8, offset=0),
        "plain weave abcd|efgh": dict(type="weave", weave_type="plain", strands="abcd|efgh", n="1"),
        "plain weave abcd-|efgh": dict(type="weave", weave_type="plain", strands="abcd-|efgh", n="1"),
        "plain weave abcd-|efgh-": dict(type="weave", weave_type="plain", strands="abcd-|efgh-", n="1"),
        "basket weave abcd|efgh": dict(type="weave", weave_type="basket", strands="abcd|efgh", n="4"),
        "basket weave abcd|efgh 2": dict(type="weave", weave_type="basket", strands="abcd|efgh", n="2"),
        "twill weave abcd|efgh": dict(type="weave", weave_type="twill", strands="abcd|efgh", n="4"),
        "twill weave abcd|efgh 2": dict(type="weave", weave_type="twill", strands="abcd|efgh", n="2"),
        "square-colouring 8": dict(type="tiling", tiling_type="square-col", n=8),
        "hex-slice 8": dict(type="tiling", tiling_type="hex-slice", n=8, offset=0),
        "hex-colouring 8": dict(type="tiling", tiling_type="hex-col", n=8),
      },
      9: {
        "cube weave abc|def|ghi": dict(type= "weave", weave_type="cube", strands="abc|def|ghi"),
        "hex-slice 9": dict(type="tiling", tiling_type="hex-slice", n=9, offset=0),
        "square-colouring 9": dict(type="tiling", tiling_type="square-col", n=9),
        "hex-colouring 9": dict(type="tiling", tiling_type="hex-col", n=9),
        "square-slice 9": dict(type="tiling", tiling_type="square-slice", n=9, offset=0),
        "archimedean 3.3.3.3.6": dict(type="tiling", tiling_type="archi", code="3.3.3.3.6"),
        "plain weave abcde|fghi": dict(type="weave", weave_type="plain", strands="abcde|fghi", n="1"),
        "plain weave abcde-|fghi": dict(type="weave", weave_type="plain", strands="abcde-|fghi", n="1"),
        "plain weave abcde-|fghi-": dict(type="weave", weave_type="plain", strands="abcde-|fghi-", n="1"),
        "hex-dissection 9": dict(type="tiling", tiling_type="hex-dissect", n=9, offset=0, offset_angle=0),
        "square-dissection 9": dict(type="tiling", tiling_type="square-dissect", n=9, offset=0, offset_angle=0),
      },
      10: {
        "hex-colouring 10": dict(type="tiling", tiling_type="hex-col", n=10),
        "hex-slice 10": dict(type="tiling", tiling_type="hex-slice", n=10, offset=0),
        "square-slice 10": dict(type="tiling", tiling_type="square-slice", n=10, offset=0),
        "plain weave abcde|fghij": dict(type="weave", weave_type="plain", strands="abcde|fghij", n="1"),
        "plain weave abcde-|fghij": dict(type="weave", weave_type="plain", strands="abcde-|fghij", n="1"),
        "plain weave abcde-|fghij-": dict(type="weave", weave_type="plain", strands="abcde-|fghij-", n="1"),
      },
      11: {
        "hex-colouring 11": dict(type="tiling", tiling_type="hex-col", n=11),
        "hex-slice 11": dict(type="tiling", tiling_type="hex-slice", n=11, offset=0),
        "square-slice 11": dict(type="tiling", tiling_type="square-slice", n=11, offset=0),
        "chavey E": dict(type="tiling", tiling_type="chavey", code="E"),
        "chavey G": dict(type="tiling", tiling_type="chavey", code="G"),
        "plain weave abcdef|ghijk": dict(type="weave", weave_type="plain", strands="abcdef|ghijk", n="1"),
        "plain weave abcdef|ghijk-": dict(type="weave", weave_type="plain", strands="abcdef|ghijk-", n="1"),
        "plain weave abcdef-|ghijk-": dict(type="weave", weave_type="plain", strands="abcdef-|ghijk-", n="1"),
      },
      12: {
        "hex-slice 12": dict(type="tiling", tiling_type="hex-slice", n=12, offset=0),
        "square-slice 12": dict(type="tiling", tiling_type="square-slice", n=12),
        "plain weave abcdef|ghijkl": dict(type="weave", weave_type="plain", strands="abcdef|ghijkl", n="1"),
        "plain weave abcdef-|ghijkl": dict(type="weave", weave_type="plain", strands="abcdef-|ghijkl", n="1"),
        "plain weave abcdef-|ghijkl-": dict(type="weave", weave_type="plain", strands="abcdef-|ghijkl-", n="1"),
      },
      # 13: {
      #   "chavey A": dict(type="tiling", tiling_type="chavey", code="A"),
      # },
      # 14: {
      #   "chavey B": dict(type="tiling", tiling_type="chavey", code="B"),
      # },
      # 15: {
      #   "chavey F": dict(type="tiling", tiling_type="chavey", code="F"),
      # },
      # 18: {
      #   "chavey C": dict(type="tiling", tiling_type="chavey", code="C"),
      #   "chavey D": dict(type="tiling", tiling_type="chavey", code="D"),
      # },
    }
    return (tilings_by_n,)


@app.cell(hide_code=True)
def _():
    centred = {
        "display": "flex",
        "height": "500px",
        "justify-content": "center",
        "align-items": "center",
        "text-align": "center",
    }
    return (centred,)


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
