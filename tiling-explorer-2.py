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
    import weavingspace as wsp
    return gpd, mpl, wsp


@app.cell(hide_code=True)
def _(gpd):
    gdf = gpd.read_file("https://raw.githubusercontent.com/DOSull/weaving-space/refs/heads/main/examples/data/dummy-data.json", engine="fiona")
    return (gdf,)


@app.cell(hide_code=True)
def _():
    def tool_tip(ele:str, tip:str):
        return f'<span title="{tip}">{ele}</span>'
    return (tool_tip,)


@app.cell(hide_code=True)
def _(gdf, mo, tile):
    _tile_ids = sorted(list(set(tile.tiles.tile_id)))
    _var_names = [_ for _ in gdf.columns if not "geom" in _]
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
    _tiled_map.render(legend=False)
    return (tiling_map,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""### Tiling modifiers""")
    return


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
    _max_t_inset = p_inset.value / 3 if tile_or_weave.value == "tiles" else tile.spacing / 10
    _value_t_inset = _max_t_inset / 3 if tile_or_weave.value == "tiles" else 0
    t_inset = mo.ui.slider(start=0, stop=2.5, step = 0.1, 
                           value=_value_t_inset, show_value=True, debounce=True)
    return (t_inset,)


@app.cell(hide_code=True)
def _(mo, p_inset, t_inset, tile_or_weave, tile_rotate, tiling_map, tool_tip):
    mo.stop(tiling_map)
    if tile_or_weave.value == "tiles":
        _str = "\n".join([f"#### Rotate by {tool_tip(tile_rotate, "Rotate tile unit (degrees)")}",
                          f"#### Prototile inset {tool_tip(p_inset, "Inset group of tiles (% spacing)")}",
                          f"#### Tile inset {tool_tip(t_inset, "Inset tiles (% spacing)")}"])
    else:
        _str = "\n".join([f"#### Rotate by {tool_tip(tile_rotate, "Rotate tile unit (degrees)")}",
                          f"#### Tile inset {tool_tip(t_inset, "Inset tiles (% spacing)")}"])
    mo.md(_str)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(f"### General tiling specification")
    return


@app.cell(hide_code=True)
def _(mo, tool_tip):
    tile_or_weave = mo.ui.dropdown(options=["tiles", "weave"], value="tiles", label="#### Pick tile or weave")
    mo.md(f"{tool_tip(tile_or_weave, "Choose tiles or a weave tiling")}")
    return (tile_or_weave,)


@app.cell(hide_code=True)
def _(mo, tile_or_weave, tool_tip):
    if tile_or_weave.value == "tiles":
        family = mo.ui.dropdown(
            options = ["hex-slice", "square-slice", "hex-dissection", "cross",
                       "laves", "archimedean", "square-colouring", "hex-colouring"], 
            value = "laves", label = "#### Tiling type") 
    else:
        family = mo.ui.dropdown(
            options = ["plain", "twill", "basket"], #, "cube"], 
            value = "plain", label = "#### Weave type")

    mo.md(f"{tool_tip(family, "Choose tiling family")}")
    return (family,)


@app.cell(hide_code=True)
def _(mo, tile_or_weave, tool_tip):
    if tile_or_weave.value == "tiles":
        spacing = mo.ui.slider(start=50, stop=5000, step=50, value=750,
                               show_value=True, debounce=True)
    else:
        spacing = mo.ui.slider(start=100, stop=1000, step=10, value=250,
                               show_value=True, debounce=True)

    mo.md(f"#### Spacing {tool_tip(spacing, 'In units of the map CRS. For tiles, size of the repeat unit; for weaves, width of the ribbons.')}")
    return (spacing,)


@app.cell(hide_code=True)
def _(family, mo, tile_or_weave):
    if "slice" in family.value:
        _n_parts = mo.ui.slider(
            steps = [x for x in range(2, 13)], value = 6,
            label=f"#### Number of slices", 
            show_value=True, debounce=True)
    elif family.value == "hex-dissection":
        _n_parts = mo.ui.slider(
            steps = [4, 7, 9], value = 7,
            label=f"#### Number of pieces", 
            show_value=True, debounce=True)
    elif family.value == "cross":
        _n_parts = mo.ui.slider(
            steps = [x for x in range(2, 8)], value = 5,
            label=f"#### Number of crosses", 
            show_value=True, debounce=True)
    elif "colour" in family.value:
        _n_parts = mo.ui.slider(
            steps = [x for x in range(2, 10)], value = 7,
            label=f"#### Number of colours",
            show_value=True, debounce=True)

    if "slice" in family.value:
        _offset = mo.ui.slider(
            steps=[x / 100 for x in range(101)], value=0,
            label="#### Offset",
            show_value=True, debounce=True) 
    elif family.value == "hex-dissection":
        _offset = mo.ui.number(
            start=0, stop=1,
            label="#### Offset", debounce=True) 

    _code = mo.ui.dropdown(
        options = ["3.3.3.3.6", "3.3.4.3.4", "3.4.6.4", "3.6.3.6", 
                  "3.12.12", "4.6.12", "4.8.8"],
        value = "3.3.4.3.4",
        label = "Code",
    )

    if family.value == "twill" or family.value == "basket":
        _strands = mo.ui.text(value="ab|cd", label="#### Strands spec") 
    elif family.value == "plain":
        _strands = mo.ui.text(value="a|b", label="#### Strands spec")
    else:
        _strands = mo.ui.text(value="abc|def|ghi", label="#### Strands spec")

    _aspect = mo.ui.slider(steps=[x / 6 for x in range(1,7)], value=5/6, label="#### Strand width",
                           show_value=True, debounce=True)

    _over_under = mo.ui.text(value="2,2" if family.value == "twill" else "2", 
                             label="#### Over-under pattern")

    if tile_or_weave.value == "tiles":
        if family.value in ["laves", "archimedean"]:
            tile_spec = mo.ui.dictionary({
                "code": _code,
            })
            tooltips = ["The sequence of vertex degrees around each polygon in the tiling (Laves), or the number of sides of each polygon around each vertex (Archimedean)."]
        elif "slice" in family.value or family.value == "hex-dissection":
            tile_spec = mo.ui.dictionary({
                "n": _n_parts,
                "offset": _offset,
            })
            tooltips = ["The number of pie-slices of the base tile.", "0 starts at the base tile corners, 1 at the mid-point along segments equally dividing the base tile perimeter."]
        else:
            tile_spec = mo.ui.dictionary({
                "n": _n_parts,
            })
            tooltips = ["How many parts to dissect the hexagon into."]
    else:
        if family.value == "plain":
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
def _(mo, tile_spec, tiling_map, tool_tip, tooltips):
    mo.stop(tiling_map)
    mo.md("\n".join([
        f"#### {tool_tip(v, tt)}" for (k, v), tt in zip(tile_spec.items(), tooltips)
    ]))
    return


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
def _(family, gdf, get_over_under, spacing, tile_or_weave, tile_spec, wsp):
    if tile_or_weave.value == "tiles":
        tile = wsp.TileUnit(
            tiling_type=family.value,
            spacing=spacing.value,
            n=tile_spec["n"].value if "n" in tile_spec else None,
            code=tile_spec["code"].value if "code" in tile_spec else None,
            offset=tile_spec["offset"].value if "offset" in tile_spec else None,
            crs=3857 if gdf is None else gdf.crs)
    else:
        tile = wsp.WeaveUnit(
            weave_type=family.value,
            spacing=spacing.value,
            strands=tile_spec["strands"].value,
            n=1 if family.value == "plain" else get_over_under(tile_spec["over_under"].value),
            aspect=tile_spec["aspect"].value,
            crs=3857 if gdf is None else gdf.crs)
    return (tile,)


@app.cell(hide_code=True)
def _(p_inset, spacing, t_inset, tile, tile_or_weave, tile_rotate):
    if tile_or_weave.value == "tiles":
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


@app.cell
def _(pals, rev_pals):
    def get_palettes():
        return [(p if not r else p + "_r") for p, r in zip(pals.value, rev_pals.value)]
    return (get_palettes,)


if __name__ == "__main__":
    app.run()
