import marimo

__generated_with = "0.11.19"
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
        mo.md(f"<span title='Weaving maps of complex data'>2025.03.13-22:15</span>").style({'background-color':'rgba(255,255,255,0.5'}),
        mo.image(src="mw.png").style(centred)
    ])
    return


@app.cell(hide_code=True)
def module_imports():
    import io
    import json
    from PIL import Image
    import matplotlib as mpl
    import math
    import pandas as pd
    import geopandas as gpd
    from shapely import is_valid
    import weavingspace as wsp
    return Image, gpd, io, is_valid, json, math, mpl, pd, wsp


@app.cell(hide_code=True)
def globals(get_colour_ramp, gpd, mo):
    dummy_data_file = "https://raw.githubusercontent.com/DOSull/weaving-space/refs/heads/main/examples/data/dummy-data.json"
    builtin_gdf = gpd.read_file(dummy_data_file, engine="fiona")
    available_palettes = [
        'Reds', 'Greys', 'Blues', 'Oranges', 'Greens', 'Purples',
        'YlGnBu', 'RdPu', 'viridis', 'summer', 'spring', 'winter']
    # make a bunch of colour ramps and save them in a dictionary
    color_ramps = {k: mo.image(get_colour_ramp(k))
                   for k in available_palettes + [p + "_r" for p in available_palettes]}
    return available_palettes, builtin_gdf, color_ramps, dummy_data_file


@app.cell(hide_code=True)
def marimo_states(available_palettes, builtin_gdf, dummy_data_file, mo):
    get_input_data, set_input_data = mo.state(dummy_data_file)
    get_gdf, set_gdf = mo.state(builtin_gdf)
    get_variables, set_variables = mo.state([])
    get_palettes, set_palettes = mo.state(available_palettes)
    get_reversed, set_reversed = mo.state([False] * 12)
    get_status_message, set_status_message = mo.state("STATUS All good!")
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
def read_gdf(
    builtin_gdf,
    get_gdf,
    get_input_data,
    get_numeric_variables,
    gpd,
    io,
    is_valid,
    set_gdf,
    set_status_message,
    set_variables,
    wsp,
):
    _done = False
    _did_repairs = False
    if type(get_input_data()) is str or len(get_input_data()) == 0:
        set_gdf(builtin_gdf)
    else:
        _old_gdf = get_gdf()
        try:
            _new_gdf = gpd.read_file(io.BytesIO(get_input_data()[0].contents).read().decode())
        except Exception as e:
            set_gdf(_old_gdf)
            set_status_message("ERROR! Exception in uploading data")
            print(e.args)
            raise
        else:
            _n = len(get_numeric_variables(_new_gdf))
            if _n < 2:
                set_status_message("WARNING! One or fewer variables, data not loaded")
                set_gdf(_old_gdf)
                _done = True
            if not _done:
                if not _new_gdf.crs.is_projected:
                    _new_gdf = _new_gdf.to_crs(3857)
                if not all(is_valid(_new_gdf.geometry)):
                    _new_gdf.geometry = wsp.tiling_utils.repair_polygon(_new_gdf.geometry)
                    if not all(is_valid(_new_gdf.geometry)):
                        set_status_message("ERROR! Geometries not valid, try another dataset")
                        set_gdf(_old_gdf)
                        _done = True
                    else:
                        _did_repairs = True
            if not _done:
                set_gdf(_new_gdf)
                set_variables(get_numeric_variables(get_gdf()))
                if _did_repairs:
                    set_status_message("WARNING! Had to repair geometry errors in your data")
                else:
                    set_status_message("STATUS All good!")
    return


@app.cell(hide_code=True)
def upload_data(mo, set_input_data, tool_tip):
    fb = mo.ui.file(filetypes=[".geojson", ".json"], 
                    on_change=set_input_data, label=f"Upload data")
    _f = tool_tip(
        fb, "Your data should be geospatial polygons - and currently only GeoJSON is readable.")
    mo.md(f"{_f}").left()
    return (fb,)


@app.cell(hide_code=True)
def _(
    centred,
    get_gdf,
    get_modded_tile_unit,
    get_selected_colour_palettes,
    get_tile_ids,
    get_variables,
    mo,
    wsp,
):
    mo.output.replace(
        mo.vstack(
            [
                mo.md("## Making tiled map..."),
                mo.md("### This might take a while"),
                mo.status.spinner(),
                mo.md(
                    "During initialisation ignore messages about missing modules"
                ),
            ]
        ).style(centred)
    )

    tiling_map = True
    _tiled_map = wsp.Tiling(get_modded_tile_unit(), get_gdf()).get_tiled_map(join_on_prototiles=False)
    _tiled_map.variables = {k: v for k, v in zip(get_tile_ids(), get_variables())}
    _tiled_map.colourmaps = {k: v for k, v in zip(get_variables(), get_selected_colour_palettes())}
    _result = _tiled_map.render(legend=False, scheme="EqualInterval")
    tiling_map = False

    _result
    return (tiling_map,)


@app.cell(hide_code=True)
def set_number_of_variables(mo, tool_tip):
    num_tiles = mo.ui.slider(steps = [x for x in range(2, 13)], value=4, debounce=True, show_value=True)
    mo.md(f"""
    ### General settings
    #### Set number of tiling elements {tool_tip(num_tiles, 'Choose the number of distinct tiles you want to use to symbolise data.')}
    """)
    return (num_tiles,)


@app.cell(hide_code=True)
def variable_palette_map_header(mo):
    mo.md("""### Variable &rarr; palette map""")
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
    if len(get_numeric_variables(get_gdf())) < num_tiles.value:
        set_status_message(f"WARNING! More tile elements ({num_tiles.value}) than variables ({len(get_numeric_variables(get_gdf()))})")

    if get_status_message() == "STATUS All good!":
        _bkgd = "lightgreen"
    elif "WARNING" in get_status_message():
        _bkgd = "pink"
    else:
        _bkgd = "red"

    _col = "white" if "ERROR" in get_status_message() else "black"

    _msg = mo.md(f"<span style='background-color:{_bkgd};color:{_col};font-face:sans-serif;padding:2px;'>{get_status_message()}</span>")

    _msg.center()
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
):
    if num_tiles.value < len(get_numeric_variables(get_gdf())):
        set_variables(get_numeric_variables(get_gdf())[:num_tiles.value])
    elif num_tiles.value > len(get_variables()):
        n_to_add = num_tiles.value - len(get_variables())
        to_add = [v for v in get_numeric_variables(get_gdf()) if v not in get_variables()][:n_to_add]
        set_variables(get_variables() + to_add)

    _chosen_palettes = get_palettes()[:num_tiles.value]
    _chosen_reversed = get_reversed()[:num_tiles.value]
    variables = mo.ui.array(
        [mo.ui.dropdown(options=get_numeric_variables(get_gdf()), value=v) for v in get_variables()], 
        label="Variables", on_change=set_variables) 
    pals = mo.ui.array(
        [mo.ui.dropdown(options=available_palettes, value=p) for p in _chosen_palettes], 
        label="Palettes", on_change=set_palettes)
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
    # mo.stop(tiling_map)
    mo.md("\n".join([
        "&nbsp;&nbsp;".join([
        f"#### Tiles `{t_id}`",
        f"{tool_tip(v, f"Variable for tiles with id {t_id}")} &rarr;",
        f"{tool_tip(p, f"Palette for variable {v.value}")}",
        f"<span style='position:relative;top:5px;'>{tool_tip(r, 'Reverse ramp')}</span>",
        f"<span style='display:inline-block;object-fit:cover;height:24px;position:relative;bottom:24px;'>{color_ramps[p.value + ("_r" if r.value else "")]}</span>",
    ]) for t_id, v, p, r in zip(get_tile_ids(), variables, pals, rev_pals)]))
    return


@app.cell(hide_code=True)
def set_spacing_limits(get_spacings, mo):
    _spacing_steps, _spacing_value = get_spacings()
    spacing = mo.ui.slider(steps=_spacing_steps, value=_spacing_value,
                          show_value=True, debounce=True)
    return (spacing,)


@app.cell
def _(mo, tile_or_weave):
    _str = "Tiling" if tile_or_weave.value == "tiling" else "Weave"
    mo.md(f"### {_str} modifiers")
    return


@app.cell(hide_code=True)
def tiling_modifier_ui_elements(mo):
    tile_rotate = mo.ui.slider(steps=[x for x in range(-90, 91, 1)], value=0, 
                               show_value=True, debounce=True)
    tile_scale_x = mo.ui.slider(steps=[x / 10 for x in range(10, 41)], value=1, 
                                show_value=True, debounce=True)
    tile_scale_y = mo.ui.slider(steps=[x / 10 for x in range(10, 41)], value=1, 
                                show_value=True, debounce=True)
    tile_skew_x = mo.ui.slider(steps=[x for x in range(-40, 41, 1)], value=0, 
                               show_value=True, debounce=True)
    tile_skew_y = mo.ui.slider(steps=[x for x in range(-40, 41, 1)], value=0, 
                               show_value=True, debounce=True)
    p_inset = mo.ui.slider(start=0, stop=10, step = 0.1, 
                           value=0, show_value=True, debounce=True)
    t_inset = mo.ui.slider(start=0, stop=5, step = 0.1, 
                           value=0, show_value=True, debounce=True)
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
    tool_tip,
):
    # mo.stop(tiling_map)
    if tile_or_weave.value == "tiling":
        _str = "\n".join([
            f"#### Spacing {tool_tip(spacing, 'In units of the map CRS, the approximate dimension of the repeating group.')}",
            f"#### Rotate by {tool_tip(tile_rotate, "Rotate tiling (degrees). Note that the tile group is rotated before any skews are applied.")}",
            f"#### Scale left-right {tool_tip(tile_scale_x, "Scale in the x direction.")}",
            f"#### Scale up-down {tool_tip(tile_scale_y, "Scale in the y direction.")}",
            f"#### Skew left-right {tool_tip(tile_skew_x, "Skew in the x direction (degrees).")}",
            f"#### Skew up-down {tool_tip(tile_skew_y, "Skew in the y direction (degrees).")}",
            f"#### Tile group inset {tool_tip(p_inset, "Inset the tile group (% spacing).")}",
            f"#### Tiles inset {tool_tip(t_inset, "Inset individual tiles (% spacing).")}"])
    else:
        _str = "\n".join([
            f"#### Spacing {tool_tip(spacing, 'In units of the map CRS, the distance between strand centre lines.')}",
            f"#### Rotate by {tool_tip(tile_rotate, "Rotate weave (degrees). Note that the weave is rotated before any skews are applied.")}",
            f"#### Scale left-right {tool_tip(tile_scale_x, "Scale in the x direction.")}",
            f"#### Scale up-down {tool_tip(tile_scale_y, "Scale in the y direction.")}",
            f"#### Skew left-right {tool_tip(tile_skew_x, "Skew in the x direction (degrees).")}",
            f"#### Skew up-down {tool_tip(tile_skew_y, "Skew in the y direction (degrees).")}",
            f"#### Strands inset {tool_tip(t_inset, "Inset strands (% width).")}"])
    mo.md(_str)
    return


@app.cell(hide_code=True)
def tiling_or_weave_chooser(mo, num_tiles, tilings_by_n, tool_tip):
    _options = list(set([v["type"] for v in tilings_by_n[num_tiles.value].values()]))
    tile_or_weave = mo.ui.dropdown(options=_options, value="tiling", label="#### Pick tiling or weave")
    mo.md(f"{tool_tip(tile_or_weave, "Choose tiling or a weave tiling")}")
    return (tile_or_weave,)


@app.cell(hide_code=True)
def tiling_type_chooser(mo, num_tiles, tile_or_weave, tilings_by_n, tool_tip):
    _type = tile_or_weave.value
    _options = [k for k, v in tilings_by_n[num_tiles.value].items() if v["type"] == _type]

    family = mo.ui.dropdown(
        options = _options, 
        value = _options[0], label = f"#### {_type.capitalize()} type") 
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
        _offset = mo.ui.slider(
            steps=[x / 100 for x in range(101)], value=0,
            label="#### Offset",
            show_value=True, debounce=True) 
    elif "hex-dissection" in family.value:
        _offset = mo.ui.number(
            start=0, stop=1,
            label="#### Offset", debounce=True)

    if "weave" in family.value:
        _aspect = mo.ui.slider(steps=[x/12 for x in range(1,13)], 
                               value=3/4, label="#### Strand width",
                               show_value=True, debounce=True)
        _over_under = mo.ui.text(value=tilings_by_n[num_tiles.value][family.value]["n"],
                                label="#### Over-under pattern")

    if tile_or_weave.value == "tiling":
        if "slice" in family.value or "dissection" in family.value:
            tile_spec = mo.ui.dictionary({"offset": _offset,})
            tooltips = ["0 starts at the base tile corners, 1 at the mid-point along segments equally dividing the base tile perimeter into the requested number of tiles."]
        else:
            tile_spec = None
    else:
        if "plain" in family.value:
            tile_spec = mo.ui.dictionary({"aspect": _aspect})
            tooltips = ["The width of the weave strands relative to strand spacing. A value of 1 will fill the map with no gaps. Progressively smaller values will leave 'holes' in the woven pattern and will eventually give the appearance of a cross hatch."]
        else:
            tile_spec = mo.ui.dictionary({
                "aspect": _aspect,
                "over_under": _over_under,
            })
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
    tool_tip,
    tooltips,
):
    # mo.stop(tiling_map)
    if tile_spec is not None:
        x = mo.md("\n".join(
            [f"### Set {tile_or_weave.value} options"] + 
            [f"#### {tool_tip(v, tt)}" for (k, v), tt in zip(tile_spec.items(), tooltips)]
        ))
    else:
        x = mo.md(f"#### No {tile_or_weave.value} options to set")
    x
    return (x,)


@app.cell(hide_code=True)
def _(mo):
    _radius = mo.ui.slider(0, 4, value=0, show_value=True)
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
    #### Show tile IDs {tool_tip(view_settings['show_ids'], 'Show the tiling element labels used to match tiles to variables in the map data.')}
    #### Show base tile {tool_tip(view_settings['show_prototile'], 'Show in fine black outline the simple tile (usually a square or hexagon) which forms the basis of the pattern.')}
    #### Show vectors {tool_tip(view_settings['show_vectors'], 'Show the translations that map repeating tiles in the pattern onto one another.')}
    #### Show &lsquo;jigsaw piece&rsquo; {tool_tip(view_settings['show_reg_prototile'], 'Show in a red outline the repeating set tile group that pieces together jigsaw-like to form the pattern.')}
    #### Show scale {tool_tip(view_settings['show_scale'], 'Give an indication of scale in map units.')}
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
    def get_base_tile_unit():
        spec = tilings_by_n[num_tiles.value][family.value]
        if tile_or_weave.value == "tiling":
            result = wsp.TileUnit(
                tiling_type=spec["tiling_type"],
                spacing=spacing.value,
                n=spec["n"] if "n" in spec else None,
                code=spec["code"] if "code" in spec else None,
                offset=tile_spec["offset"].value if tile_spec is not None else None,
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
        slice = [i in get_tile_ids() for i in result.tiles.tile_id]
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
def _():
    def get_over_under(pattern):
        if any([not c in "0123456789," for c in pattern]): return (2, 2)
        numbers = [int(s) for s in pattern.split(",")]
        length = 2 * len(numbers) // 2
        if length == 0:
            return (2, 2)
        else:
            return tuple(numbers[:length])
    return (get_over_under,)


@app.cell(hide_code=True)
def _(get_selected_colour_palettes, mpl, num_tiles, view_settings):
    def plot_tiles(tiles):
        cols = [mpl.colormaps.get(p)(2/3) for p in get_selected_colour_palettes()]
        while len(cols) < num_tiles.value:
            cols = cols + ["#ffffff00"]
        cm = mpl.colors.ListedColormap(cols)
        plot = tiles.plot(r=view_settings["radius"].value, 
                          show_vectors=view_settings["show_vectors"].value, 
                          show_prototile=view_settings["show_prototile"].value,
                          show_reg_prototile=view_settings["show_reg_prototile"].value,
                          show_ids=view_settings["show_ids"].value,
                          cmap=cm, figsize=(4.25, 4.25), r_alpha=0.5)
        if view_settings["show_scale"].value:
            plot.xaxis.set_visible(True)
            plot.yaxis.set_visible(False)
            plot.set_frame_on(False)
        else:
            plot.axis("off")
        return plot
    return (plot_tiles,)


@app.cell(hide_code=True)
def _(get_gdf, get_numeric_variables):
    def get_tile_ids():
        # return sorted(list(set(get_modded_tile_unit().tiles.tile_id)))
        # return list("abcdefghijkl")[:num_tiles.value]
        return list("abcdefghijkl")[:len(get_numeric_variables(get_gdf()))]
    return (get_tile_ids,)


@app.cell(hide_code=True)
def _(num_tiles, pals, rev_pals):
    def get_selected_colour_palettes():
        return [(p if not r else p + "_r") 
                for p, r in zip(pals.value[:num_tiles.value], 
                                rev_pals.value[:num_tiles.value])]
    return (get_selected_colour_palettes,)


@app.cell(hide_code=True)
def _(io, mpl):
    def get_colour_ramp(pal_name:str="Reds", rev:bool=False):
        fig, ax = mpl.pyplot.subplots()
        xy = [[x / 256 for x in range(257)] for i in range(2)]
        ax.imshow(xy, aspect=32, cmap=mpl.colormaps.get(pal_name + ("_r" if rev else "")))
        ax.set_axis_off()
        buf = io.BytesIO()
        mpl.pyplot.savefig(buf, dpi=24, pad_inches=0, bbox_inches="tight")
        mpl.pyplot.close(fig)
        return buf
    return (get_colour_ramp,)


@app.cell(hide_code=True)
def _(pd):
    def get_numeric_variables(_gdf):
        return [col for col in _gdf.columns if not "geom" in col 
                and pd.api.types.is_numeric_dtype(_gdf[col].dtype)]
    return (get_numeric_variables,)


@app.cell(hide_code=True)
def _(get_gdf, math, tile_or_weave):
    def get_spacings():
        _bb = get_gdf().total_bounds
        _width, _height = _bb[2] - _bb[0], _bb[3] - _bb[1]
        _max = 10 ** math.floor(math.log10(max(_width, _height))) // 5
        _mid = _max // 2
        _min = _max // 20
        _stepsize1 = _min // 10
        _stepsize2 = _min // 2
        _steps = [x for x in range(_min, _mid, _stepsize1)] + \
                 [x for x in range(_mid, _max + 1, _stepsize2)]
        if tile_or_weave.value == "weave":
            _steps = [x // 2 for x in _steps]
            _mid = _mid // 2
        return _steps, 3 * _mid // 4
    return (get_spacings,)


@app.cell(hide_code=True)
def _():
    def tool_tip(ele:str, tip:str):
        return f'<span title="{tip}">{ele}</span>'
    return (tool_tip,)


@app.cell(hide_code=True)
def setup_tilings_dictionary():
    tilings_by_n = {
      2: {
        "plain weave a|b": {"type":"weave", "weave_type": "plain", "strands": "a|b", "n": "1"},
        "twill weave a|b": {"type":"weave", "weave_type": "twill", "strands": "a|b", "n": "2"},
        "twill weave a|b-": {"type":"weave", "weave_type": "twill", "strands": "a|b-", "n": "2"},
        "twill weave a|b 3": {"type":"weave", "weave_type": "twill", "strands": "a|b", "n": "3"},
        "twill weave a|b- 3": {"type":"weave", "weave_type": "twill", "strands": "a|b-", "n": "3"},
        "twill weave a|b 4": {"type":"weave", "weave_type": "twill", "strands": "a|b", "n": "4"},
        "twill weave a|b 1,2": {"type":"weave", "weave_type": "twill", "strands": "a|b", "n": "1,2"},
        "twill weave a|b- 1,2": {"type":"weave", "weave_type": "twill", "strands": "a|b-", "n": "1,2"},
        "twill weave a|b 1,2,2,1": {"type":"weave", "weave_type": "twill", "strands": "a|b", "n": "1,2,2,1"},
        "twill weave a|b 2,3": {"type":"weave", "weave_type": "twill", "strands": "a|b", "n": "2,3"},
        "twill weave a|b 2,3,3,2": {"type":"weave", "weave_type": "twill", "strands": "a|b", "n": "2,3,3,2"},
        "archimedean 4.8.8": {"type":"tiling", "tiling_type": "archimedean", "code": "4.8.8"},
        "square-slice 2": {"type":"tiling", "tiling_type": "square-slice", "n": 2},
        "crosses 2": {"type":"tiling", "tiling_type": "crosses", "n": 2},
        "hex-colouring 2": {"type":"tiling", "tiling_type": "hex-colouring", "n": 2},
        "square-colouring 2": {"type":"tiling", "tiling_type": "square-colouring", "n": 2},
        "hex-slice 2": {"type":"tiling", "tiling_type": "hex-slice", "n": 2},
      },
      3: {
        "plain weave ab|c": {"type":"weave", "weave_type": "plain", "strands": "ab|c", "n": "1"},
        "plain weave ab-|c": {"type":"weave", "weave_type": "plain", "strands": "ab-|c", "n": "1"},
        "twill weave ab|c-": {"type":"weave", "weave_type": "twill", "strands": "ab|c-", "n": "2"},
        "basket weave ab|c-": {"type":"weave", "weave_type": "basket", "strands": "ab|c-", "n": "2"},
        "hex-slice 3": {"type":"tiling", "tiling_type": "hex-slice", "n": 3},
        # "laves 3.6.3.6": {"type":"tiling", "tiling_type": "laves", "code": "3.6.3.6"},
        "hex-colouring 3": {"type":"tiling", "tiling_type": "hex-colouring", "n": 3},
        "crosses 3": {"type":"tiling", "tiling_type": "crosses", "n": 3},
        "square-colouring 3": {"type":"tiling", "tiling_type": "square-colouring", "n": 3},
        "archimedean 3.6.3.6": {"type":"tiling", "tiling_type": "archimedean", "code": "3.6.3.6"},
        "archimedean 3.12.12": {"type":"tiling", "tiling_type": "archimedean", "code": "3.12.12"},
        "square-slice 3": {"type":"tiling", "tiling_type": "square-slice", "n": 3},
        # "twill weave ab|c-": {"type":"weave", "weave_type": "twill", "strands": "ab|c-", "n": "2"},
        # "basket weave ab|c-": {"type":"weave", "weave_type": "basket", "strands": "ab|c-", "n": "2"},
      },
      4: {
        "laves 3.3.4.3.4": {"type":"tiling", "tiling_type": "laves", "code": "3.3.4.3.4"},
        "plain weave ab|cd": {"type":"weave", "weave_type": "plain", "strands": "ab|cd", "n": "1"},
        "plain weave ab-|cd": {"type":"weave", "weave_type": "plain", "strands": "ab-|cd", "n": "1"},
        "plain weave ab-|cd-": {"type":"weave", "weave_type": "plain", "strands": "ab-|cd-", "n": "1"},
        "basket weave ab|cd": {"type":"weave", "weave_type": "basket", "strands": "ab|cd", "n": "2"},
        "basket weave ab|cd 3": {"type":"weave", "weave_type": "basket", "strands": "ab|cd", "n": "3"},
        "twill weave ab|cd": {"type":"weave", "weave_type": "twill", "strands": "ab|cd", "n": "2"},
        "twill weave ab|cd 3": {"type":"weave", "weave_type": "twill", "strands": "ab|cd", "n": "3"},
        "twill weave ab|cd 1,2": {"type":"weave", "weave_type": "twill", "strands": "ab|cd", "n": "1,2"},
        "twill weave ab|cd 1,2,2,1": {"type":"weave", "weave_type": "twill", "strands": "ab|cd", "n": "1,2,2,1"},
        "basket weave ab-|cd": {"type":"weave", "weave_type": "basket", "strands": "ab-|cd", "n": "3"},
        "twill weave ab-|cd": {"type":"weave", "weave_type": "twill", "strands": "ab-|cd", "n": "3"},
        "basket weave ab-|cd-": {"type":"weave", "weave_type": "basket", "strands": "ab-|cd-", "n": "3"},
        "twill weave ab-|cd-": {"type":"weave", "weave_type": "twill", "strands": "ab-|cd-", "n": "3"},
        # "laves 4.8.8": {"type":"tiling", "tiling_type": "laves", "code": "4.8.8"},
        "crosses 4": {"type":"tiling", "tiling_type": "crosses", "n": 4},
        "square-slice 4": {"type":"tiling", "tiling_type": "square-slice", "n": 4},
        "hex-colouring 4": {"type":"tiling", "tiling_type": "hex-colouring", "n": 4},
        # "square-colouring 4": {"type":"tiling", "tiling_type": "square-colouring", "n": 4},
        "hex-dissection 4": {"type":"tiling", "tiling_type": "hex-dissection", "n": 4},
        "hex-slice 4": {"type":"tiling", "tiling_type": "hex-slice", "n": 4},
      },
      5: {
        "square-colouring 5": {"type":"tiling", "tiling_type": "square-colouring", "n": 5},
        "crosses 5": {"type":"tiling", "tiling_type": "crosses", "n": 5},
        "plain weave abc|de": {"type":"weave", "weave_type": "plain", "strands": "abc|de", "n": "1"},
        "plain weave abc-|de": {"type":"weave", "weave_type": "plain", "strands": "abc-|de", "n": "1"},
        "plain weave abc-|de-": {"type":"weave", "weave_type": "plain", "strands": "abc-|de-", "n": "1"},
        "twill weave abc|de": {"type":"weave", "weave_type": "twill", "strands": "abc|de", "n": "3"},
        "twill weave abc|de-": {"type":"weave", "weave_type": "twill", "strands": "abc|de-", "n": "3"},
        "twill weave abc-|de-": {"type":"weave", "weave_type": "twill", "strands": "abc-|de-", "n": "3"},
        "basket weave abc|de": {"type":"weave", "weave_type": "basket", "strands": "abc|de", "n": "3"},
        "basket weave abc|de-": {"type":"weave", "weave_type": "basket", "strands": "abc|de-", "n": "3"},
        "basket weave abc-|de-": {"type":"weave", "weave_type": "basket", "strands": "abc-|de-", "n": "3"},
        "hex-colouring 5": {"type":"tiling", "tiling_type": "hex-colouring", "n": 5},
        "hex-slice 5": {"type":"tiling", "tiling_type": "hex-slice", "n": 5},
        "square-slice 5": {"type":"tiling", "tiling_type": "square-slice", "n": 5},
      },
      6: {
        "hex-slice 6": {"type":"tiling", "tiling_type": "hex-slice", "n": 6},
        "square-slice 6": {"type":"tiling", "tiling_type": "square-slice", "n": 6},
        "square-colouring 6": {"type":"tiling", "tiling_type": "square-colouring", "n": 6},
        "laves 3.3.3.3.6": {"type":"tiling", "tiling_type": "laves", "code": "3.3.3.3.6"},
        # "laves 3.4.6.4": {"type":"tiling", "tiling_type": "laves", "code": "3.4.6.4"},
        "laves 3.12.12": {"type":"tiling", "tiling_type": "laves", "code": "3.12.12"},
        "plain weave abc|def": {"type":"weave", "weave_type": "plain", "strands": "abc|def", "n": "1"},
        "plain weave abc-|def": {"type":"weave", "weave_type": "plain", "strands": "abc-|def", "n": "1"},
        "plain weave abc-|def-": {"type":"weave", "weave_type": "plain", "strands": "abc-|def-", "n": "1"},
        "basket weave abc|def": {"type":"weave", "weave_type": "basket", "strands": "abc|def", "n": "3"},
        "twill weave abc|def": {"type":"weave", "weave_type": "twill", "strands": "abc|def", "n": "3"},
        "basket weave abc-|def": {"type":"weave", "weave_type": "basket", "strands": "abc-|def", "n": "3"},
        "twill weave abc-|def": {"type":"weave", "weave_type": "twill", "strands": "abc-|def", "n": "3"},
        "basket weave abc-|def-": {"type":"weave", "weave_type": "basket", "strands": "abc-|def-", "n": "4"},
        "twill weave abc-|def-": {"type":"weave", "weave_type": "twill", "strands": "abc-|def-", "n": "4"},
        "archimedean 3.3.4.3.4": {"type":"tiling", "tiling_type": "archimedean", "code": "3.3.4.3.4"},
        "archimedean 3.4.6.4": {"type":"tiling", "tiling_type": "archimedean", "code": "3.4.6.4"},
        "archimedean 4.6.12": {"type":"tiling", "tiling_type": "archimedean", "code": "4.6.12"},
        "crosses 6": {"type":"tiling", "tiling_type": "crosses", "n": 6},
        "hex-colouring 6": {"type":"tiling", "tiling_type": "hex-colouring", "n": 6},
      },
      7: {
        "hex-colouring 7": {"type":"tiling", "tiling_type": "hex-colouring", "n": 7},
        "crosses 7": {"type":"tiling", "tiling_type": "crosses", "n": 7},
        "hex-dissection 7": {"type":"tiling", "tiling_type": "hex-dissection", "n": 7},
        "plain weave abcd|efg": {"type":"weave", "weave_type": "plain", "strands": "abcd|efg", "n": "1"},
        "plain weave abcd|efg-": {"type":"weave", "weave_type": "plain", "strands": "abcd|efg-", "n": "1"},
        "plain weave abcd-|efg-": {"type":"weave", "weave_type": "plain", "strands": "abcd-|efg-", "n": "1"},
        "basket weave abcd|efg": {"type":"weave", "weave_type": "basket", "strands": "abcd|efg", "n": "3"},
        "basket weave abcd|efg-": {"type":"weave", "weave_type": "basket", "strands": "abcd|efg-", "n": "4"},
        "basket weave abcd|efg- 2": {"type":"weave", "weave_type": "basket", "strands": "abcd|efg-", "n": "2"},
        "twill weave abcd|efg": {"type":"weave", "weave_type": "twill", "strands": "abcd|efg", "n": "3"},
        "twill weave abcd|efg-": {"type":"weave", "weave_type": "twill", "strands": "abcd|efg-", "n": "4"},
        "twill weave abcd|efg- 2": {"type":"weave", "weave_type": "twill", "strands": "abcd|efg-", "n": "2"},
        "square-colouring 7": {"type":"tiling", "tiling_type": "square-colouring", "n": 7},
        "hex-slice 7": {"type":"tiling", "tiling_type": "hex-slice", "n": 7},
        "square-slice 7": {"type":"tiling", "tiling_type": "square-slice", "n": 7},
      },
      8: {
        "square-slice 8": {"type":"tiling", "tiling_type": "square-slice", "n": 8},
        "plain weave abcd|efgh": {"type":"weave", "weave_type": "plain", "strands": "abcd|efgh", "n": "1"},
        "plain weave abcd-|efgh": {"type":"weave", "weave_type": "plain", "strands": "abcd-|efgh", "n": "1"},
        "plain weave abcd-|efgh-": {"type":"weave", "weave_type": "plain", "strands": "abcd-|efgh-", "n": "1"},
        "basket weave abcd|efgh": {"type":"weave", "weave_type": "basket", "strands": "abcd|efgh", "n": "4"},
        "basket weave abcd|efgh 2": {"type":"weave", "weave_type": "basket", "strands": "abcd|efgh", "n": "2"},
        "twill weave abcd|efgh": {"type":"weave", "weave_type": "twill", "strands": "abcd|efgh", "n": "4"},
        "twill weave abcd|efgh 2": {"type":"weave", "weave_type": "twill", "strands": "abcd|efgh", "n": "2"},
        "square-colouring 8": {"type":"tiling", "tiling_type": "square-colouring", "n": 8},
        "hex-slice 8": {"type":"tiling", "tiling_type": "hex-slice", "n": 8},
        "hex-colouring 8": {"type":"tiling", "tiling_type": "hex-colouring", "n": 8},
      },
      9: {
        # "cube weave abc|def|ghi": {"type":"weave", "weave_type": "cube", "strands": "abc|def|ghi"},
        "hex-slice 9": {"type":"tiling", "tiling_type": "hex-slice", "n": 9},
        "square-colouring 9": {"type":"tiling", "tiling_type": "square-colouring", "n": 9},
        "hex-colouring 9": {"type":"tiling", "tiling_type": "hex-colouring", "n": 9},
        "square-slice 9": {"type":"tiling", "tiling_type": "square-slice", "n": 9},
        "archimedean 3.3.3.3.6": {"type":"tiling", "tiling_type": "archimedean", "code": "3.3.3.3.6"},
        "plain weave abcde|fghi": {"type":"weave", "weave_type": "plain", "strands": "abcde|fghi", "n": "1"},
        "plain weave abcde-|fghi": {"type":"weave", "weave_type": "plain", "strands": "abcde-|fghi", "n": "1"},
        "plain weave abcde-|fghi-": {"type":"weave", "weave_type": "plain", "strands": "abcde-|fghi-", "n": "1"},
      },
      10: {
        "hex-slice 10": {"type":"tiling", "tiling_type": "hex-slice", "n": 10},
        "square-slice 10": {"type":"tiling", "tiling_type": "square-slice", "n": 10},
        "plain weave abcde|fghij": {"type":"weave", "weave_type": "plain", "strands": "abcde|fghij", "n": "1"},
        "plain weave abcde-|fghij": {"type":"weave", "weave_type": "plain", "strands": "abcde-|fghij", "n": "1"},
        "plain weave abcde-|fghij-": {"type":"weave", "weave_type": "plain", "strands": "abcde-|fghij-", "n": "1"},
      },
      11: {
        "hex-slice 11": {"type":"tiling", "tiling_type": "hex-slice", "n": 11},
        "square-slice 11": {"type":"tiling", "tiling_type": "square-slice", "n": 11},
        "plain weave abcdef|ghijk": {"type":"weave", "weave_type": "plain", "strands": "abcdef|ghijk", "n": "1"},
        "plain weave abcdef|ghijk-": {"type":"weave", "weave_type": "plain", "strands": "abcdef|ghijk-", "n": "1"},
        "plain weave abcdef-|ghijk-": {"type":"weave", "weave_type": "plain", "strands": "abcdef-|ghijk-", "n": "1"},
      },
      12: {
        "hex-slice 12": {"type":"tiling", "tiling_type": "hex-slice", "n": 12},
        # "laves 4.6.12": {"type":"tiling", "tiling_type": "laves", "code": "4.6.12"},
        "square-slice 12": {"type":"tiling", "tiling_type": "square-slice", "n": 12},
        "plain weave abcdef|ghijkl": {"type":"weave", "weave_type": "plain", "strands": "abcdef|ghijkl", "n": "1"},
        "plain weave abcdef-|ghijkl": {"type":"weave", "weave_type": "plain", "strands": "abcdef-|ghijkl", "n": "1"},
        "plain weave abcdef-|ghijkl-": {"type":"weave", "weave_type": "plain", "strands": "abcdef-|ghijkl-", "n": "1"},
      },
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
    centred_100 = {
        "display": "flex",
        "width": "100%",
        "justify-content": "center",
        "align-items": "center",
        "text-align": "center",
    }
    return centred, centred_100


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
