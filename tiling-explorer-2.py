import marimo

__generated_with = "0.11.9"
app = marimo.App(
    width="medium",
    layout_file="layouts/tiling-explorer-2.grid.json",
)


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(f"# MapWeaver ~ tiled maps of multivariate data")
    return


@app.cell(hide_code=True)
def _():
    import matplotlib as mpl
    import geopandas as gpd
    from pandas.api.types import is_numeric_dtype
    import weavingspace as wsp
    return gpd, is_numeric_dtype, mpl, wsp


@app.cell(hide_code=True)
def _():
    # gdf = gpd.read_file("https://raw.githubusercontent.com/DOSull/weaving-space/refs/heads/main/examples/data/dummy-data.json", engine="fiona")
    return


@app.cell(hide_code=True)
def _():
    def tool_tip(ele:str, tip:str):
        return f'<span title="{tip}">{ele}</span>'
    return (tool_tip,)


@app.cell(hide_code=True)
def _(gdf, is_numeric_dtype, mo, num_tiles):
    _tile_ids = list("abcdefghijkl")[:num_tiles.value]
    _var_names = [_ for _ in gdf.columns if not "geom" in _ and is_numeric_dtype(gdf[_].dtype)]
    _pal_names = ['Reds', 'Greens', 'Greys', 'Blues', 'Oranges', 'Purples', 
                  'YlGnBu', 'RdPu', 'viridis', 'summer', 'spring', 'winter']

    vars = mo.ui.array([mo.ui.dropdown(options=_var_names, value=_var_names[i]) 
                        for i, id in enumerate(_tile_ids)], label="Variables") 
    pals = mo.ui.array([mo.ui.dropdown(options=_pal_names, value=_pal_names[i]) 
                        for i in range(len(_tile_ids))], label="Palettes")
    rev_pals = mo.ui.array([mo.ui.switch(False) for i in range(len(_tile_ids))])
    return pals, rev_pals, vars


@app.cell(hide_code=True)
def _(mo):
    get_input_data, set_input_data = mo.state("https://raw.githubusercontent.com/DOSull/weaving-space/refs/heads/main/examples/data/dummy-data.json")
    return get_input_data, set_input_data


@app.cell(hide_code=True)
def _(get_input_data, gpd):
    def get_gdf():
        if type(get_input_data()) is str:
            _gdf = gpd.read_file(get_input_data(), engine="fiona")
        else:
            _gdf = gpd.read_file(get_input_data()[0].path, engine="fiona")
        if not _gdf.crs.is_projected:
            _gdf = _gdf.to_crs(3857)
        return _gdf

    gdf = get_gdf()
    _bb = gdf.total_bounds
    _w, _h = _bb[2] - _bb[0], _bb[3] - _bb[1]
    _a = _w * _h
    _min_spacing = max(_w, _h) // 200
    spacing_lims = {
        "min spacing": _min_spacing,
        "max spacing": _min_spacing * 10,
        "step spacing": (_min_spacing * 9) // 100,
        "min strand width": _min_spacing // 3,
        "max strand width": _min_spacing // 3  * 10,
        "step strand width": (_min_spacing * 3) // 100,
    }
    return gdf, get_gdf, spacing_lims


@app.cell(hide_code=True)
def _(mo, set_input_data, tool_tip):
    f = tool_tip(
        mo.ui.file_browser(multiple=False, on_change=set_input_data, label=f"### Select an input data set"),
        "Your data should be geospatial polygons - preferably GPKG or GeoJSON, and contain a number of numerical attributes to be symbolised.")
    mo.md(f"{f}")
    return (f,)


@app.cell
def _(mo):
    mo.md(f"### General settings")
    return


@app.cell(hide_code=True)
def _(mo, tool_tip):
    num_tiles = mo.ui.dropdown(options = {str(x): x for x in range(2, 13)}, value="6")
    mo.md(f"#### Set number of variables {tool_tip(num_tiles, "Choose the number of distinct tiles you want to use to symbolise data.")}")
    return (num_tiles,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""### Variable to palette mapping""")
    return


@app.cell(hide_code=True)
def _(get_palettes, mo, mpl, pals, tiling_map):
    mo.stop(tiling_map)
    _n = len(pals)
    _fig, _axs = mpl.pyplot.subplots(nrows = _n + 1, figsize=(1.5, 0.25 + 0.48 * _n))
    for ax, cm, in zip(_axs[1:], get_palettes()):
        _xy = [[x / 256 for x in range(257)], [x / 256 for x in range(257)]]
        ax.imshow(_xy, aspect='auto', cmap=mpl.colormaps.get(cm))
    for ax in _axs:
        ax.set_axis_off()
    ax
    return ax, cm


@app.cell(hide_code=True)
def _(mo, pals, rev_pals, tile, tiling_map, tool_tip, vars):
    mo.stop(tiling_map)
    _tile_ids = sorted(list(set((tile.tiles.tile_id))))
    mo.md("\n".join([
        "&nbsp;&nbsp;".join([
            f"#### Tiles `{t_id}`",
            f"{tool_tip(v, f"Variable for tiles with id {t_id}")} &rarr;",
            f"{tool_tip(p, f"Palette for variable {v.value}")}",
            f"<span style='position:relative;top:5px;'>{tool_tip(r, 'Reverse ramp')}</span>",
        ])
        for t_id, v, p, r in zip(_tile_ids, vars, pals, rev_pals)
    ]))
    return


@app.cell(hide_code=True)
def _():
    # tile_map_button = mo.ui.run_button(label="Tile map!")
    # mo.md(f"{tool_tip(tile_map_button, "Click to tile the map data")}").center().style({"padding": "5px"})
    return


@app.cell(hide_code=True)
def _(final_tile, gdf, get_palettes, mo, tile, vars, wsp):
    _centred = {"display": "flex", 
                "height": "500px", 
                "justify-content": "center", 
                "align-items": "center", 
                "text-align": "center"}

    mo.output.replace(
        mo.vstack([
            mo.md("## Making tiled map..."),
            mo.md("### This might take a while"),
            mo.status.spinner()
        ]).style(_centred)
    )

    # _btn_text = f"Tile map!"
    # _btn = f"<span style='background-color:#efefef;font-family:monospace;padding:3px;box-shadow:3px 3px #888;'>{_btn_text}</span>"

    # _msg = "\n".join([
    #     f"## Click {_btn} above left", 
    #     f"## to start, and also to remake the", 
    #     f"## map after design changes",
    # ])
    # mo.stop(not tile_map_button.value, mo.md(_msg).style(_centred))

    # _tile_ids = sorted(list(set((tile.tiles.tile_id))))
    # tiled_map = wsp.Tiling(final_tile, gdf).get_tiled_map()
    # tiled_map.variables = {k: v for k, v in zip(_tile_ids, vars.value)}
    # tiled_map.colourmaps = {k: v for k, v in zip(vars.value, get_palettes())}
    # tiled_map.render(legend=False)_tile_ids = sorted(list(set((tile.tiles.tile_id))))

    tiling_map = True
    _tile_ids = sorted(list(set((tile.tiles.tile_id))))
    _tiled_map = wsp.Tiling(final_tile, gdf).get_tiled_map()
    _tiled_map.variables = {k: v for k, v in zip(_tile_ids, vars.value)}
    _tiled_map.colourmaps = {k: v for k, v in zip(vars.value, get_palettes())}
    tiling_map = False
    _tiled_map.render(legend=False, scheme="EqualInterval")
    return (tiling_map,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""### Tiling modifiers""")
    return


@app.cell(hide_code=True)
def _(mo, spacing_lims, tile_or_weave):
    if tile_or_weave.value == "tiling":
        spacing = mo.ui.slider(
            start=spacing_lims['min spacing'], 
            stop=spacing_lims['max spacing'], 
            step=spacing_lims['step spacing'], 
            value=spacing_lims['max spacing'] // 2,
            show_value=True, debounce=True)
    else:
        spacing = mo.ui.slider(
            start=spacing_lims['min strand width'], 
            stop=spacing_lims['max strand width'], 
            step=spacing_lims['step strand width'], 
            value=spacing_lims['max strand width'] // 2,
            show_value=True, debounce=True)
    # mo.md(f"#### Spacing {tool_tip(spacing, 'In units of the map CRS. For tiles, size of the repeat unit; for weaves, width of the ribbons.')}")
    return (spacing,)


@app.cell(hide_code=True)
def _(mo):
    tile_rotate = mo.ui.slider(start=-90, stop=90, step=5, value=0, show_value=True, debounce=True)
    return (tile_rotate,)


@app.cell(hide_code=True)
def _(mo):
    p_inset = mo.ui.slider(start=0, stop=10, step = 0.1, 
                           value=0, show_value=True, debounce=True)
    return (p_inset,)


@app.cell(hide_code=True)
def _(mo, p_inset, tile, tile_or_weave):
    _max_t_inset = p_inset.value / 3 if tile_or_weave.value == "tiling" else tile.spacing / 10
    _value_t_inset = _max_t_inset / 3 if tile_or_weave.value == "tiling" else 0
    t_inset = mo.ui.slider(start=0, stop=2.5, step = 0.1, 
                           value=_value_t_inset, show_value=True, debounce=True)
    return (t_inset,)


@app.cell(hide_code=True)
def _(
    mo,
    p_inset,
    spacing,
    t_inset,
    tile_or_weave,
    tile_rotate,
    tiling_map,
    tool_tip,
):
    mo.stop(tiling_map)
    if tile_or_weave.value == "tiling":
        _str = "\n".join([
            f"#### Spacing {tool_tip(spacing, 'In units of the map CRS. For tiles, size of the repeat unit; for weaves, width of the ribbons.')}",
            f"#### Rotate by {tool_tip(tile_rotate, "Rotate tile unit (degrees)")}",
            f"#### Prototile inset {tool_tip(p_inset, "Inset group of tiles (% spacing)")}",
            f"#### Tile inset {tool_tip(t_inset, "Inset tiles (% spacing)")}"])
    else:
        _str = "\n".join([
            f"#### Spacing {tool_tip(spacing, 'In units of the map CRS. For tiles, size of the repeat unit; for weaves, width of the ribbons.')}",
            f"#### Rotate by {tool_tip(tile_rotate, "Rotate tile unit (degrees)")}",
            f"#### Tile inset {tool_tip(t_inset, "Inset tiles (% spacing)")}"])
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
def _(family, mo, num_tiles, tile_or_weave, tilings_by_n):
    _chosen_family = family.value

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
        _strands = mo.ui.text(value=tilings_by_n[num_tiles.value][family.value]["strands"], 
                              label="#### Strands spec") 
        _aspect = mo.ui.slider(steps=[x / 6 for x in range(1,7)], value=5/6, label="#### Strand width",
                               show_value=True, debounce=True)
        _over_under = mo.ui.text(value="2,2" if "twill" in family.value else "2", 
                                 label="#### Over-under pattern")

    if tile_or_weave.value == "tiling":
        if "slice" in family.value or "dissection" in family.value:
            tile_spec = mo.ui.dictionary({
                "offset": _offset,
            })
            tooltips = ["0 starts at the base tile corners, 1 at the mid-point along segments equally dividing the base tile perimeter."]
        else:
            tile_spec = None
    else:
        if "plain" in family.value:
            tile_spec = mo.ui.dictionary({
                "strands": _strands,
                "aspect": _aspect,
            })
            tooltips = ["A code expressing the sequence of distinct weave elements as a series of letters in the warp and weft directions, separated by a pipe | symbol.", "The width of the weave ribbons relative to spacing."]
        else:
            tile_spec = mo.ui.dictionary({
                "strands": _strands,
                "aspect": _aspect,
                "over_under": _over_under,
            })
            tooltips = ["A code expressing the sequence of distinct weave elements as a series of letters in the warp and weft directions, separated by a pipe | symbol.", "The width of the weave ribbons relative to spacing.", "A comma-separated sequence of how many strands in the other direction each ribbon should go over, then under."]
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
    _spec = tilings_by_n[num_tiles.value][family.value]
    if tile_or_weave.value == "tiling":
        tile = wsp.TileUnit(
            tiling_type=_spec["tiling_type"],
            spacing=spacing.value,
            n=_spec["n"] if "n" in _spec else None,
            code=_spec["code"] if "code" in _spec else None,
            offset=tile_spec["offset"].value if tile_spec is not None else None,
            crs=3857 if gdf is None else gdf.crs)
    else:
        tile = wsp.WeaveUnit(
            weave_type=_spec["weave_type"],
            spacing=spacing.value,
            strands=tile_spec["strands"].value,
            n=1 if "plain" in family.value else get_over_under(tile_spec["over_under"].value),
            aspect=tile_spec["aspect"].value,
            crs=3857 if gdf is None else gdf.crs)
    return (tile,)


@app.cell(hide_code=True)
def _(p_inset, spacing, t_inset, tile, tile_or_weave, tile_rotate):
    if tile_or_weave.value == "tiling":
        final_tile = tile \
           .transform_rotate(tile_rotate.value) \
           .inset_tiles(t_inset.value * spacing.value / 100) \
           .inset_prototile(p_inset.value * spacing.value / 100)
    else:
        final_tile = tile \
           .transform_rotate(tile_rotate.value) \
           .inset_tiles(t_inset.value * spacing.value / 100)
    return (final_tile,)


@app.cell(hide_code=True)
def _(final_tile, plot_tiles):
    plot_tiles(final_tile)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(f"### Tile design view settings")
    return


@app.cell(hide_code=True)
def _(mo):
    _radius = mo.ui.slider(0, 4, value=0, show_value=True)
    _show_prototile = mo.ui.switch(value=False)
    _show_reg_prototile = mo.ui.switch(value=True)
    _show_vectors = mo.ui.switch(value=False)
    _show_ids = mo.ui.switch(value=False)

    view_settings = mo.ui.dictionary({
        "radius": _radius,
        "show_prototile": _show_prototile,
        "show_reg_prototile": _show_reg_prototile,
        "show_vectors": _show_vectors,
        "show_ids": _show_ids
    })
    return (view_settings,)


@app.cell(hide_code=True)
def _(mo, tool_tip, view_settings):
    mo.md("\n".join([
       f"#### Repeats to show {tool_tip(view_settings['radius'], 'The number of &lsquo;shells&rsquo; of the tiling to show away from the base tile unit.')}",
       f"#### Show tile IDs {tool_tip(view_settings['show_ids'], 'Show the tile labels used to match tiling elements to variables in the map data.')}",
       f"#### Show &lsquo;jigsaw piece&rsquo; {tool_tip(view_settings['show_reg_prototile'], 'Show in a red outline the repeating set of tiles that piece together jigsaw-like to form the pattern.')}",
       f"#### Show vectors {tool_tip(view_settings['show_vectors'], 'Show the translations that map repeating tiles in the pattern onto one another.')}",
       f"#### Show base tile {tool_tip(view_settings['show_prototile'], 'Show in fine black outline the repeating base tile (usually a square or hexagon) which forms the basis of the pattern.')}",
    ]))
    return


@app.cell(hide_code=True)
def _(get_palettes, mpl, view_settings):
    def plot_tiles(tiles):
        cm = mpl.colors.ListedColormap([mpl.colormaps.get(p)(2/3) for p in get_palettes()])
        plot = tiles.plot(r=view_settings["radius"].value, 
                          show_vectors=view_settings["show_vectors"].value, 
                          show_prototile=view_settings["show_prototile"].value,
                          show_reg_prototile=view_settings["show_reg_prototile"].value,
                          show_ids=view_settings["show_ids"].value,
                          cmap=cm)
        plot.axis("off")
        return plot
    return (plot_tiles,)


@app.cell(hide_code=True)
def _(pals, rev_pals):
    def get_palettes():
        return [(p if not r else p + "_r") for p, r in zip(pals.value, rev_pals.value)]
    return (get_palettes,)


@app.cell(hide_code=True)
def _():
    tilings_by_n = {
      2: {
        "plain weave a|b": {"type":"weave", "weave_type": "plain", "strands": "a|b"},
        "twill weave a|b": {"type":"weave", "weave_type": "twill", "strands": "a|b"},
        "basket weave a|b": {"type":"weave", "weave_type": "basket", "strands": "a|b"},
        "archimedean 4.8.8": {"type":"tiling", "tiling_type": "archimedean", "code": "4.8.8"},
        "square-slice 2": {"type":"tiling", "tiling_type": "square-slice", "n": 2},
        "crosses 2": {"type":"tiling", "tiling_type": "crosses", "n": 2},
        "hex-colouring 2": {"type":"tiling", "tiling_type": "hex-colouring", "n": 2},
        "square-colouring 2": {"type":"tiling", "tiling_type": "square-colouring", "n": 2},
        "hex-slice 2": {"type":"tiling", "tiling_type": "hex-slice", "n": 2},
      },
      3: {
        "hex-slice 3": {"type":"tiling", "tiling_type": "hex-slice", "n": 3},
        "laves 3.6.3.6": {"type":"tiling", "tiling_type": "laves", "code": "3.6.3.6"},
        "hex-colouring 3": {"type":"tiling", "tiling_type": "hex-colouring", "n": 3},
        "crosses 3": {"type":"tiling", "tiling_type": "crosses", "n": 3},
        "square-colouring 3": {"type":"tiling", "tiling_type": "square-colouring", "n": 3},
        "archimedean 3.6.3.6": {"type":"tiling", "tiling_type": "archimedean", "code": "3.6.3.6"},
        "archimedean 3.12.12": {"type":"tiling", "tiling_type": "archimedean", "code": "3.12.12"},
        "square-slice 3": {"type":"tiling", "tiling_type": "square-slice", "n": 3},
      },
      4: {
        "laves 3.3.4.3.4": {"type":"tiling", "tiling_type": "laves", "code": "3.3.4.3.4"},
        "basket weave ab|cd": {"type":"weave", "weave_type": "basket", "strands": "ab|cd"},
        "laves 4.8.8": {"type":"tiling", "tiling_type": "laves", "code": "4.8.8"},
        "crosses 4": {"type":"tiling", "tiling_type": "crosses", "n": 4},
        "square-slice 4": {"type":"tiling", "tiling_type": "square-slice", "n": 4},
        "hex-colouring 4": {"type":"tiling", "tiling_type": "hex-colouring", "n": 4},
        "square-colouring 4": {"type":"tiling", "tiling_type": "square-colouring", "n": 4},
        "hex-dissection 4": {"type":"tiling", "tiling_type": "hex-dissection", "n": 4},
        "hex-slice 4": {"type":"tiling", "tiling_type": "hex-slice", "n": 4},
      },
      5: {
        "square-colouring 5": {"type":"tiling", "tiling_type": "square-colouring", "n": 5},
        "crosses 5": {"type":"tiling", "tiling_type": "crosses", "n": 5},
        "twill weave abc|cd-": {"type":"weave", "weave_type": "twill", "strands": "abc|cd-"},
        "hex-colouring 5": {"type":"tiling", "tiling_type": "hex-colouring", "n": 5},
        "hex-slice 5": {"type":"tiling", "tiling_type": "hex-slice", "n": 5},
        "square-slice 5": {"type":"tiling", "tiling_type": "square-slice", "n": 5},
      },
      6: {
        "hex-slice 6": {"type":"tiling", "tiling_type": "hex-slice", "n": 6},
        "square-slice 6": {"type":"tiling", "tiling_type": "square-slice", "n": 6},
        "square-colouring 6": {"type":"tiling", "tiling_type": "square-colouring", "n": 6},
        "laves 3.3.3.6": {"type":"tiling", "tiling_type": "laves", "code": "3.3.3.6"},
        "laves 3.4.6.4": {"type":"tiling", "tiling_type": "laves", "code": "3.4.6.4"},
        "laves 3.12.12": {"type":"tiling", "tiling_type": "laves", "code": "3.12.12"},
        "basket weave abc|def": {"type":"weave", "weave_type": "basket", "strands": "abc|def"},
        "archimedean 3.3.4.3.4": {"type":"tiling", "tiling_type": "archimedean", "code": "3.3.4.3.3"},
        "archimedean 3.4.6.4": {"type":"tiling", "tiling_type": "archimedean", "code": "3.4.6.4"},
        "archimedean 4.6.12": {"type":"tiling", "tiling_type": "archimedean", "code": "4.6.12"},
        "crosses 6": {"type":"tiling", "tiling_type": "crosses", "n": 6},
        "hex-colouring 6": {"type":"tiling", "tiling_type": "hex-colouring", "n": 6},
      },
      7: {
        "hex-colouring 7": {"type":"tiling", "tiling_type": "hex-colouring", "n": 7},
        "crosses 7": {"type":"tiling", "tiling_type": "crosses", "n": 7},
        "hex-dissection 7": {"type":"tiling", "tiling_type": "hex-dissection", "n": 7},
        "twill weave abc|def-": {"type":"weave", "weave_type": "twill", "strands": "abc|def-"},
        "square-colouring 7": {"type":"tiling", "tiling_type": "square-colouring", "n": 7},
        "hex-slice 7": {"type":"tiling", "tiling_type": "hex-slice", "n": 7},
        "square-slice 7": {"type":"tiling", "tiling_type": "square-slice", "n": 7},
      },
      8: {
        "square-slice 8": {"type":"tiling", "tiling_type": "square-slice", "n": 8},
        "basket weave abcd|efgh": {"type":"weave", "weave_type": "basket", "strands": "abcd|efgh"},
        "square-colouring 8": {"type":"tiling", "tiling_type": "square-colouring", "n": 8},
        "hex-slice 8": {"type":"tiling", "tiling_type": "hex-slice", "n": 8},
        "hex-colouring 8": {"type":"tiling", "tiling_type": "hex-colouring", "n": 8},
      },
      9: {
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
        "laves 4.6.12": {"type":"tiling", "tiling_type": "laves", "code": "4.6.12"},
        "square-slice 12": {"type":"tiling", "tiling_type": "square-slice", "n": 12},
      },
    }
    return (tilings_by_n,)


if __name__ == "__main__":
    app.run()
