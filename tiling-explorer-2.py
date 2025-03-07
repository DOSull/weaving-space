import marimo

__generated_with = "0.11.13"
app = marimo.App(
    width="medium",
    app_title="MapWeaver",
    layout_file="layouts/tiling-explorer-2.grid.json",
    html_head_file="splash.html",
)


@app.cell(hide_code=True)
def _(mo):
    mo.hstack([
        mo.md(f"# MapWeaver ~ tiled maps of complex data"),
        mo.md("v2025.03.08-09:20")
    ]).center()
    return


@app.cell(hide_code=True)
def module_imports():
    import io
    from PIL import Image
    import matplotlib as mpl
    import math
    import pandas as pd
    import geopandas as gpd
    import fiona
    import weavingspace as wsp
    return Image, fiona, gpd, io, math, mpl, pd, wsp


@app.cell(hide_code=True)
def globals(gpd):
    dummy_data_file = "https://raw.githubusercontent.com/DOSull/weaving-space/refs/heads/main/examples/data/dummy-data.json"
    builtin_gdf = gpd.read_file(dummy_data_file, engine="fiona")
    available_palettes = [
        'Reds', 'Greys', 'Blues', 'Oranges', 'Greens', 'Purples',
        'YlGnBu', 'RdPu', 'viridis', 'summer', 'spring', 'winter']
    return available_palettes, builtin_gdf, dummy_data_file


@app.cell
def marimo_states(available_palettes, dummy_data_file, mo):
    get_input_data, set_input_data = mo.state(dummy_data_file)
    get_palettes, set_palettes = mo.state(available_palettes)
    get_reversed, set_reversed = mo.state([False] * 12)
    get_status_message, set_status_message = mo.state("STATUS All good!")
    return (
        get_input_data,
        get_palettes,
        get_reversed,
        get_status_message,
        set_input_data,
        set_palettes,
        set_reversed,
        set_status_message,
    )


@app.cell(hide_code=True)
def read_gdf(
    builtin_gdf,
    dummy_data_file,
    fiona,
    get_input_data,
    get_numeric_variables,
    gpd,
    io,
    set_input_data,
):
    if type(get_input_data()) is str or len(get_input_data()) == 0:
        gdf = builtin_gdf
    else:
        assert (len(get_input_data()) > 0), f"No uploaded file has been specified."
        try:
            _new_gdf = gpd.read_file(io.BytesIO(get_input_data()[0].contents), engine="fiona", mode="r")
        # except FileNotFoundError as e:
        #     print(e)
        #     raise
        # except OSError as e:
        #     print(e)
        #     raise
        # except RuntimeError as e:
        #     print(e)
        #     raise
        except fiona.errors.DriverError as e:
            print(e)
            raise
        else:
            _n = len(get_numeric_variables(_new_gdf))
            if _n < 2:
                set_input_data(dummy_data_file)
                gdf = builtin_gdf
            else:
                if not _new_gdf.crs.is_projected:
                    _new_gdf = _new_gdf.to_crs(3857)
                gdf = _new_gdf
    return (gdf,)


@app.cell(hide_code=True)
def upload_data(mo, set_input_data, tool_tip):
    # _fb = mo.ui.file_browser(multiple=False, 
    #                          on_change=set_input_data, 
    #                          label=f"### Select an input data set")
    fb = mo.ui.file(on_change=set_input_data, label=f"Upload data")
    _f = tool_tip(
        fb, "Your data should be geospatial polygons - preferably GPKG or GeoJSON, and contain a number of numerical attributes to be symbolised.")
    mo.md(f"{_f}").center()
    return (fb,)


@app.cell(hide_code=True)
def tile_the_map(
    centred,
    gdf,
    get_modded_tile_unit,
    get_selected_colour_palettes,
    get_tile_ids,
    mo,
    vars,
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
    _tiled_map = wsp.Tiling(get_modded_tile_unit(), gdf).get_tiled_map()
    _tiled_map.variables = {k: v for k, v in zip(get_tile_ids(), vars.value)}
    _tiled_map.colourmaps = {k: v for k, v in zip(vars.value, get_selected_colour_palettes())}
    tiling_map = False
    _tiled_map.render(legend=False, scheme="EqualInterval", figsize=(9, 6))
    return (tiling_map,)


@app.cell(hide_code=True)
def set_number_of_variables(mo, tool_tip):
    num_tiles = mo.ui.slider(steps = [x for x in range(2, 13) if x != 11], value=4, debounce=True, show_value=True)
    mo.md(f"""
    ### General settings
    #### Set number of tiling elements {tool_tip(num_tiles, 'Choose the number of distinct tiles you want to use to symbolise data.')}
    """)
    return (num_tiles,)


@app.cell(hide_code=True)
def build_variable_and_palette_dropdowns(
    available_palettes,
    gdf,
    get_palettes,
    get_reversed,
    mo,
    num_tiles,
    pd,
    set_palettes,
    set_reversed,
    set_status_message,
):
    available_vars = [col for col in gdf.columns if not "geom" in col 
                       and pd.api.types.is_numeric_dtype(gdf[col].dtype)]

    if len(available_vars) < num_tiles.value:
        set_status_message(f"WARNING! More tile elements ({num_tiles.value}) than variables ({len(available_vars)})")

    _chosen_vars = available_vars
    _chosen_palettes = get_palettes()[:len(available_vars)]
    _chosen_reversed = get_reversed()[:len(available_vars)]
    vars = mo.ui.array(
        [mo.ui.dropdown(options=available_vars, value=_chosen_vars[i]) 
         for i, id in enumerate(_chosen_vars)], label="Variables") 
    pals = mo.ui.array(
        [mo.ui.dropdown(options=available_palettes, value=_chosen_palettes[i]) 
         for i, id in enumerate(_chosen_palettes)], label="Palettes", on_change=set_palettes)
    rev_pals = mo.ui.array(
        [mo.ui.switch(_chosen_reversed[i]) 
         for i, id in enumerate(_chosen_reversed)], on_change=set_reversed)
    return available_vars, pals, rev_pals, vars


@app.cell
def variable_palette_map_header(mo):
    mo.md("""### Variable &rarr; palette map""")
    return


@app.cell
def status_panel(get_status_message, mo):
    _bkgd = "lightgreen" if get_status_message() == "STATUS All good!" else "pink"
    _warning = mo.md(f"<span style='background-color:{_bkgd};font-face:sans-serif;padding:2px;'>{get_status_message()}</span>")

    _warning.center()
    return


@app.cell(hide_code=True)
def build_var_palette_mapping(
    get_colour_ramp,
    get_tile_ids,
    mo,
    pals,
    rev_pals,
    tiling_map,
    tool_tip,
    vars,
):
    mo.stop(tiling_map)
    mo.md("\n".join([
        "&nbsp;&nbsp;".join([
        f"#### Tiles `{t_id}`",
        f"{tool_tip(v, f"Variable for tiles with id {t_id}")} &rarr;",
        f"{tool_tip(p, f"Palette for variable {v.value}")}",
        f"<span style='position:relative;top:5px;'>{tool_tip(r, 'Reverse ramp')}</span>",
        f"<span style='display:inline-block;object-fit:cover;height:24px;position:relative;bottom:24px;'>{mo.image(get_colour_ramp(p.value, r.value))}</span>",
    ]) for t_id, v, p, r in zip(get_tile_ids(), vars, pals, rev_pals)]))
    return


@app.cell(hide_code=True)
def draw_colour_ramps():
    # mo.stop(tiling_map)
    # _fig, _axs = mpl.pyplot.subplots(nrows = num_tiles.value + 1, 
    #                                  figsize=(1.5, 0.25 + 0.48 * num_tiles.value))
    # for ax, cm, in zip(_axs[1:], get_selected_colour_palettes()):
    #     _xy = [[x / 256 for x in range(257)], [x / 256 for x in range(257)]]
    #     ax.imshow(_xy, aspect='auto', cmap=mpl.colormaps.get(cm))
    # for ax in _axs:
    #     ax.set_axis_off()
    # ax
    return


@app.cell(hide_code=True)
def set_spacing_limits(get_spacings, mo):
    # _spacings = get_spacings()
    # if tile_or_weave.value == "tiling":
    #     spacing = mo.ui.slider(
    #         start=_spacings['min spacing'], 
    #         stop=_spacings['max spacing'], 
    #         step=_spacings['step spacing'], 
    #         value=_spacings['max spacing'] // 2,
    #         show_value=True, debounce=True)
    # else:
    #     spacing = mo.ui.slider(
    #         start=_spacings['min strand width'], 
    #         stop=_spacings['max strand width'], 
    #         step=_spacings['step strand width'], 
    #         value=_spacings['max strand width'] // 2,
    #         show_value=True, debounce=True)
    _spacing_steps, _spacing_value = get_spacings()
    spacing = mo.ui.slider(steps=_spacing_steps, value=_spacing_value,
                          show_value=True, debounce=True)
    return (spacing,)


@app.cell(hide_code=True)
def _(mo):
    tile_rotate = mo.ui.slider(start=-90, stop=90, step=5, value=0, show_value=True, debounce=True)
    tile_skew_x = mo.ui.slider(start=-40, stop=40, step=5, value=0, show_value=True, debounce=True)
    tile_skew_y = mo.ui.slider(start=-40, stop=40, step=5, value=0, show_value=True, debounce=True)
    return tile_rotate, tile_skew_x, tile_skew_y


@app.cell(hide_code=True)
def _(mo):
    p_inset = mo.ui.slider(start=0, stop=10, step = 0.1, 
                           value=0, show_value=True, debounce=True)
    return (p_inset,)


@app.cell(hide_code=True)
def _(mo, p_inset, spacing, tile_or_weave):
    _max_t_inset = p_inset.value / 3 if tile_or_weave.value == "tiling" else spacing.value / 10
    _value_t_inset = _max_t_inset / 3 if tile_or_weave.value == "tiling" else 0
    t_inset = mo.ui.slider(start=0, stop=2.5, step = 0.1, 
                           value=_value_t_inset, show_value=True, debounce=True)
    return (t_inset,)


@app.cell(hide_code=True)
def setup_tiling_modifiers(
    mo,
    p_inset,
    spacing,
    t_inset,
    tile_or_weave,
    tile_rotate,
    tile_skew_x,
    tile_skew_y,
    tiling_map,
    tool_tip,
):
    mo.stop(tiling_map)
    if tile_or_weave.value == "tiling":
        _str = "\n".join([
            f"### Tiling modifiers",
            f"#### Spacing {tool_tip(spacing, 'In units of the map CRS, the approximate dimension of the repeating group')}",
            f"#### Rotate by {tool_tip(tile_rotate, "Rotate tiling (degrees). Note that the tile group is rotated before any skews are applied.")}",
            f"#### Skew left-right {tool_tip(tile_skew_x, "Skew in the x direction (degrees)")}",
            f"#### Skew up-down {tool_tip(tile_skew_y, "Skew in the y direction (degrees)")}",
            f"#### Tile group inset {tool_tip(p_inset, "Inset the tile group (% spacing)")}",
            f"#### Tiles inset {tool_tip(t_inset, "Inset individual tiles (% spacing)")}"])
    else:
        _str = "\n".join([
            f"### Weave modifiers",
            f"#### Spacing {tool_tip(spacing, 'In units of the map CRS, the distance between strand centre lines.')}",
            f"#### Rotate by {tool_tip(tile_rotate, "Rotate weave (degrees). Note that the weave is rotated before any skews are applied.")}",
            f"#### Skew left-right {tool_tip(tile_skew_x, "Skew in the x direction (degrees)")}",
            f"#### Skew up-down {tool_tip(tile_skew_y, "Skew in the y direction (degrees)")}",
            f"#### Strands inset {tool_tip(t_inset, "Inset strands (% spacing)")}"])
    mo.md(_str)
    return


@app.cell(hide_code=True)
def _(mo, num_tiles, tilings_by_n, tool_tip):
    _options = list(set([v["type"] for v in tilings_by_n[num_tiles.value].values()]))
    tile_or_weave = mo.ui.dropdown(options=_options, value="tiling", label="#### Pick tiling or weave")
    mo.md(f"{tool_tip(tile_or_weave, "Choose tiling or a weave tiling")}")
    return (tile_or_weave,)


@app.cell(hide_code=True)
def _(mo, num_tiles, tile_or_weave, tilings_by_n, tool_tip):
    _type = tile_or_weave.value
    _options = [k for k, v in tilings_by_n[num_tiles.value].items() if v["type"] == _type]

    family = mo.ui.dropdown(
        options = _options, 
        value = _options[0], label = f"#### {_type.capitalize()} type") 
    mo.md(f"{tool_tip(family, "Choose tiling family")}")
    return (family,)


@app.cell(hide_code=True)
def _(get_modded_tile_unit, plot_tiles):
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
                               value=1, label="#### Strand width",
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
def _(mo, tile_or_weave, tile_spec, tiling_map, tool_tip, tooltips):
    mo.stop(tiling_map)
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
def _(mo, tool_tip, view_settings):
    mo.md(f"""
    ### Design view options
    #### Tile group 'shells' to show {tool_tip(view_settings['radius'], 'The number of &lsquo;shells&rsquo; of the tiling to show around the base tile group.')}
    #### Show tile IDs {tool_tip(view_settings['show_ids'], 'Show the tiling element labels used to match tiles to variables in the map data.')}
    #### Show &lsquo;jigsaw piece&rsquo; {tool_tip(view_settings['show_reg_prototile'], 'Show in a red outline the repeating set tile group that pieces together jigsaw-like to form the pattern.')}
    #### Show vectors {tool_tip(view_settings['show_vectors'], 'Show the translations that map repeating tiles in the pattern onto one another.')}
    #### Show base tile {tool_tip(view_settings['show_prototile'], 'Show in fine black outline the simple tile (usually a square or hexagon) which forms the basis of the pattern.')}
    #### Show scale {tool_tip(view_settings['show_scale'], 'Give an indication of scale in map units.')}
    """)
    return


@app.cell(hide_code=True)
def _(
    family,
    gdf,
    get_over_under,
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
            return wsp.TileUnit(
                tiling_type=spec["tiling_type"],
                spacing=spacing.value,
                n=spec["n"] if "n" in spec else None,
                code=spec["code"] if "code" in spec else None,
                offset=tile_spec["offset"].value if tile_spec is not None else None,
                crs=3857 if gdf is None else gdf.crs)
        else:
            return wsp.WeaveUnit(
                weave_type=spec["weave_type"],
                spacing=spacing.value,
                strands=spec["strands"],
                n=get_over_under(tile_spec["over_under"].value) \
                    if spec["weave_type"] in ["twill", "basket"] else 1,
                aspect=tile_spec["aspect"].value,
                crs=3857 if gdf is None else gdf.crs)
    return (get_base_tile_unit,)


@app.cell(hide_code=True)
def _(
    get_base_tile_unit,
    p_inset,
    spacing,
    t_inset,
    tile_or_weave,
    tile_rotate,
    tile_skew_x,
    tile_skew_y,
):
    def get_modded_tile_unit():
        if tile_or_weave.value == "tiling":
            return get_base_tile_unit() \
               .transform_rotate(tile_rotate.value) \
               .transform_skew(tile_skew_x.value, tile_skew_y.value) \
               .inset_tiles(t_inset.value * spacing.value / 100) \
               .inset_prototile(p_inset.value * spacing.value / 100)
        else:
            return get_base_tile_unit() \
               .transform_rotate(tile_rotate.value) \
               .transform_skew(tile_skew_x.value, tile_skew_y.value) \
               .inset_tiles(t_inset.value * spacing.value / 100)
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
            cols = cols + ["#00000000"]
        cm = mpl.colors.ListedColormap(cols)
        plot = tiles.plot(r=view_settings["radius"].value, 
                          show_vectors=view_settings["show_vectors"].value, 
                          show_prototile=view_settings["show_prototile"].value,
                          show_reg_prototile=view_settings["show_reg_prototile"].value,
                          show_ids=view_settings["show_ids"].value,
                          cmap=cm, figsize=(4, 4))
        if view_settings["show_scale"].value:
            plot.xaxis.set_visible(True)
            plot.yaxis.set_visible(False)
            plot.set_frame_on(False)
        else:
            plot.axis("off")
        return plot
    return (plot_tiles,)


@app.cell(hide_code=True)
def _(get_modded_tile_unit, num_tiles):
    def get_tile_ids():
        return sorted(list(set(get_modded_tile_unit().tiles.tile_id)))
        return list("abcdefghijkl")[:num_tiles.value]
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
        return buf
    return (get_colour_ramp,)


@app.cell(hide_code=True)
def _(pd):
    def get_numeric_variables(_gdf):
        return [col for col in _gdf.columns if not "geom" in col 
                and pd.api.types.is_numeric_dtype(_gdf[col].dtype)]
    return (get_numeric_variables,)


@app.cell(hide_code=True)
def _(gdf, math, tile_or_weave):
    def get_spacings():
        _bb = gdf.total_bounds
        _width, _height = _bb[2] - _bb[0], _bb[3] - _bb[1]
        print((_width, _height))
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
        "basket weave a|b": {"type":"weave", "weave_type": "basket", "strands": "a|b", "n": "2"},
        "archimedean 4.8.8": {"type":"tiling", "tiling_type": "archimedean", "code": "4.8.8"},
        "square-slice 2": {"type":"tiling", "tiling_type": "square-slice", "n": 2},
        "crosses 2": {"type":"tiling", "tiling_type": "crosses", "n": 2},
        "hex-colouring 2": {"type":"tiling", "tiling_type": "hex-colouring", "n": 2},
        "square-colouring 2": {"type":"tiling", "tiling_type": "square-colouring", "n": 2},
        "hex-slice 2": {"type":"tiling", "tiling_type": "hex-slice", "n": 2},
      },
      3: {
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
        "basket weave ab|cd": {"type":"weave", "weave_type": "basket", "strands": "ab|cd", "n": "2"},
        "twill weave ab|cd": {"type":"weave", "weave_type": "twill", "strands": "ab|cd", "n": "2"},
        "basket weave ab|cd-": {"type":"weave", "weave_type": "basket", "strands": "ab|cd-", "n": "3"},
        "twill weave ab|cd-": {"type":"weave", "weave_type": "twill", "strands": "ab|cd-", "n": "3"},
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
        "twill weave abc|de": {"type":"weave", "weave_type": "twill", "strands": "abc|de", "n": "3"},
        "twill weave abc|de-": {"type":"weave", "weave_type": "twill", "strands": "abc|de-", "n": "3"},
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
        "basket weave abc|def": {"type":"weave", "weave_type": "basket", "strands": "abc|def", "n": "3"},
        "twill weave abc|def": {"type":"weave", "weave_type": "twill", "strands": "abc|def", "n": "3"},
        "basket weave abc|def-": {"type":"weave", "weave_type": "basket", "strands": "abc|def-", "n": "4"},
        "twill weave abc|def-": {"type":"weave", "weave_type": "twill", "strands": "abc|def-", "n": "4"},
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
        "twill weave abcd|efg": {"type":"weave", "weave_type": "twill", "strands": "abcd|defg", "n": "4"},
        "twill weave abcd|efg-": {"type":"weave", "weave_type": "twill", "strands": "abcd|defg-", "n": "4"},
        "square-colouring 7": {"type":"tiling", "tiling_type": "square-colouring", "n": 7},
        "hex-slice 7": {"type":"tiling", "tiling_type": "hex-slice", "n": 7},
        "square-slice 7": {"type":"tiling", "tiling_type": "square-slice", "n": 7},
      },
      8: {
        "square-slice 8": {"type":"tiling", "tiling_type": "square-slice", "n": 8},
        "basket weave abcd|efgh": {"type":"weave", "weave_type": "basket", "strands": "abcd|efgh", "n": "4"},
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
      },
      10: {
        "hex-slice 10": {"type":"tiling", "tiling_type": "hex-slice", "n": 10},
        "square-slice 10": {"type":"tiling", "tiling_type": "square-slice", "n": 10},
      },
      11: {
        "hex-slice 11": {"type":"tiling", "tiling_type": "hex-slice", "n": 11},
        "square-slice 11": {"type":"tiling", "tiling_type": "square-slice", "n": 11},
      },
      12: {
        "hex-slice 12": {"type":"tiling", "tiling_type": "hex-slice", "n": 12},
        # "laves 4.6.12": {"type":"tiling", "tiling_type": "laves", "code": "4.6.12"},
        "square-slice 12": {"type":"tiling", "tiling_type": "square-slice", "n": 12},
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
    return (centred,)


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
