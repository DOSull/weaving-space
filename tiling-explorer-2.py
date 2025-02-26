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
def _():
    import asyncio
    import math
    return asyncio, math


@app.cell(hide_code=True)
def _():
    import matplotlib as mpl
    import geopandas as gpd
    from weavingspace import TileUnit
    from weavingspace import WeaveUnit
    from weavingspace import Tiling
    return TileUnit, Tiling, WeaveUnit, gpd, mpl


@app.cell(hide_code=True)
def _(gpd):
    gdf = gpd.read_file("https://raw.githubusercontent.com/DOSull/weaving-space/refs/heads/main/examples/data/dummy-data.json", engine="fiona")
    return (gdf,)


@app.cell(hide_code=True)
def _(Tiling, gdf, tile):
    tiled_map = Tiling(tile, gdf).get_tiled_map()
    return (tiled_map,)


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
    return pals, vars


@app.cell(hide_code=True)
def _(mo, pals, tile, vars):
    mo.md("\n".join([f"### Variable to palette mapping"] +
        [f"#### Tiles `{t_id}` {v} &rarr; {p}"
         for t_id, v, p in zip(sorted(list(set((tile.tiles.tile_id)))), vars, pals)]))
    return


@app.cell(hide_code=True)
def _(mpl, pals):
    _n = len(pals)
    _fig, _axs = mpl.pyplot.subplots(nrows = _n + 1, figsize=(2, 0.55 + 0.39 * _n))
    for ax, cm in zip(_axs[1:], pals.value):
        _xy = [[x / 256 for x in range(257)], [x / 256 for x in range(257)]]
        ax.imshow(_xy, aspect='auto', cmap=mpl.colormaps.get(cm))
    for ax in _axs:
        ax.set_axis_off()
    ax
    return ax, cm


@app.cell(hide_code=True)
def _(mo):
    tile_map_button = mo.ui.run_button(label="Tile map!")
    tile_map_button
    return (tile_map_button,)


@app.cell(hide_code=True)
async def _(asyncio, mo, pals, tile, tile_map_button, tiled_map, vars):
    _centred = {"display": "flex", "height": "500px", 
                "justify-content": "center", "align-items": "center", 
                "text-align": "center"}

    mo.output.replace(
        mo.vstack([
            mo.md("## Making tiled map..."),
            mo.md("### This might take a while"),
            mo.status.spinner()
        ]).style(_centred)
    )
    await asyncio.sleep(1)

    _msg = "\n".join([
        "## Click the **Tile map!** button to start", 
        "## or to refresh the map",
        "## after design changes",
    ])
    mo.stop(not tile_map_button.value, mo.md(_msg).style(_centred))
    tiled_map.variables = {k: v for k, v in zip(list(set(tile.tiles.tile_id)), vars.value)}
    tiled_map.colourmaps = {k: v for k, v in zip(vars.value, pals.value)}
    tiled_map.render(legend=False)
    return


@app.cell(hide_code=True)
def _(mo, tile_or_weave):
    _tile_rotate = mo.ui.slider(start=-90, stop=90, step=5, value=0, show_value=True, debounce=True)
    if tile_or_weave.value == "tiles":
        _spacing = mo.ui.slider(start=50, stop=5000, step=50, value=750, show_value=True, debounce=True)
    else:
        _spacing = mo.ui.slider(start=100, stop=1000, step=10, value=250, show_value=True, debounce=True)
    spacing_rotation = mo.ui.dictionary({
      "tile_rotate": _tile_rotate,
      "spacing": _spacing
    })
    return (spacing_rotation,)


@app.cell(hide_code=True)
def _(mo, spacing_rotation):
    mo.md("\n".join([
        f"### Tiling settings",
        f"#### Spacing {spacing_rotation["spacing"]}",
        f"#### Rotate by {spacing_rotation["tile_rotate"]}"
    ]))
    return


@app.cell(hide_code=True)
def _(mo, spacing_rotation, tile_or_weave):
    _max_prototile_spacing = spacing_rotation["spacing"].value // 6
    p_inset = mo.ui.slider(start=0, stop=_max_prototile_spacing, step = 1, 
                           value=0, show_value=True, debounce=True)
    _str = f"#### Prototile inset {p_inset}" if tile_or_weave.value == "tiles" else ""
    mo.md(_str)
    return (p_inset,)


@app.cell(hide_code=True)
def _(mo, p_inset, tile_or_weave):
    _max_tile_spacing = p_inset.value // 3 if tile_or_weave.value == "tiles" else 10
    t_inset = mo.ui.slider(start=0, stop=_max_tile_spacing, step = 1, 
                           value=0, show_value=True, debounce=True)
    mo.md(f"#### Tile inset {t_inset}")
    return (t_inset,)


@app.cell(hide_code=True)
def _(mo):
    tile_or_weave = mo.ui.dropdown(options=["tiles", "weave"], value="tiles")
    mo.md(f"#### Pick tile or weave {tile_or_weave}")
    return (tile_or_weave,)


@app.cell(hide_code=True)
def _(mo, tile_or_weave):
    if tile_or_weave.value == "tiles":
        type = mo.ui.dropdown(
            options = ["hex-slice", "square-slice", "hex-dissection", "cross",
                       "laves", "archimedean", "square-colouring", "hex-colouring"], 
            value = "laves") 
        _str = "Tiling type"
    else:
        type = mo.ui.dropdown(
            options = ["plain", "twill", "basket", "cube"], 
            value = "plain")
        _str = "Weave type"

    mo.md(f"#### {_str} {type}")
    return (type,)


@app.cell(hide_code=True)
def _(mo, tile_or_weave, type):
    if "slice" in type.value:
        _n = [x for x in range(2, 13)]
        _v = 6
        _prefix = "Number of slices: "
    elif type.value == "hex-dissection":
        _n = [4, 7, 9]
        _v = 7
        _prefix = "Number of pieces: "
    elif type.value == "cross":
        _n = [x for x in range(2, 8)]
        _v = 5
        _prefix = "Number of crosses: "
    elif "colour" in type.value:
        _n = [x for x in range(2, 10)]
        _v = 7
        _prefix = "Number of colours: "
    else:
        _n = [0, 1]
        _prefix = "NA"
        _v = 0

    n = mo.ui.slider(steps = _n, value = _v, show_value=True, debounce=True)

    if "slice" in type.value:
        _offset = [x / 100 for x in range(101)]
        _v_offset = 0
    else:
        _offset = [0, 1]
        _v_offset = 0

    offset = mo.ui.slider(steps = _offset, value = _v_offset, show_value=True)

    code = mo.ui.dropdown(
        options = ["3.3.3.3.6", "3.3.4.3.4", "3.4.6.4", "3.6.3.6", 
                  "3.12.12", "4.6.12", "4.8.8"],
        value = "3.3.4.3.4"     
    )

    if type.value == "twill" or type.value == "basket":
        _strands = "ab|cd"
    elif type.value == "plain":
        _strands = "a|b"
    else:
        _strands = "abc|def|ghi"

    strands = mo.ui.text(value=_strands)
    aspect = mo.ui.slider(steps=[1/3,0.4,0.5,2/3,3/4,4/5,9/10,1], value=3/4, show_value=True)

    if tile_or_weave.value == "tiles":
        if type.value in ["laves", "archimedean"]:
            _str = f"#### Code: {code}"
        else:
            _str = f"#### {_prefix} (`n`) {n}"
            if "slice" in type.value or type.value == "hex-dissection":
                _str = "\n".join([_str, f"#### Offset: {offset}"])
    else:
        _str = "\n".join([
            f"#### Strands spec: {strands}",
            f"#### Aspect: {aspect}"
        ])

    mo.md(_str)
    return aspect, code, n, offset, strands


@app.cell(hide_code=True)
def _(
    TileUnit,
    WeaveUnit,
    aspect,
    code,
    gdf,
    n,
    offset,
    spacing_rotation,
    strands,
    tile_or_weave,
    type,
):
    if tile_or_weave.value == "tiles":
        tile = TileUnit(tiling_type=type.value, 
                        n = n.value, 
                        code = code.value, 
                        offset = offset.value,
                        spacing=spacing_rotation["spacing"].value, 
                        crs=gdf.crs)
    else:
        tile = WeaveUnit(weave_type=type.value, 
                         strands=strands.value, 
                         spacing=spacing_rotation["spacing"].value, 
                         aspect=aspect.value, 
                         crs=gdf.crs)
    return (tile,)


@app.cell(hide_code=True)
def _(p_inset, plot_tiles, spacing_rotation, t_inset, tile, tile_or_weave):
    if tile_or_weave.value == "tiles":
        _plot_tile = tile \
           .transform_rotate(spacing_rotation["tile_rotate"].value) \
           .inset_tiles(t_inset.value) \
           .inset_prototile(p_inset.value)
    else:
        _plot_tile = tile \
           .transform_rotate(spacing_rotation["tile_rotate"].value) \
           .inset_tiles(t_inset.value)

    plot_tiles(_plot_tile)
    return


@app.cell(hide_code=True)
def _(mo):
    _radius = mo.ui.slider(0, 4, value=1, show_value=True)
    _show_prototile = mo.ui.switch(value=False)
    _show_reg_prototile = mo.ui.switch(value=False)
    _show_vectors = mo.ui.switch(value=False)
    _show_ids = mo.ui.switch(value=True)

    view_settings = mo.ui.dictionary({
        "radius": _radius,
        "show_prototile": _show_prototile,
        "show_reg_prototile": _show_reg_prototile,
        "show_vectors": _show_vectors,
        "show_ids": _show_ids
    })
    return (view_settings,)


@app.cell(hide_code=True)
def _(mo, view_settings):
    mo.md("\n".join(["### Tile design view settings",
      f"#### Set radius {view_settings["radius"]}",
      f"#### Show vectors {view_settings["show_vectors"]}",
      f"#### Show base tile {view_settings["show_prototile"]}",
      f"#### Show prototile {view_settings["show_reg_prototile"]}",
      f"#### Show tile IDs {view_settings["show_ids"]}"]))
    return


@app.cell(hide_code=True)
def _(mpl, pals, view_settings):
    def plot_tiles(tiles):
        cm = mpl.colors.ListedColormap([mpl.colormaps.get(p)(0.75) for p in pals.value])
        plot = tiles.plot(r=view_settings["radius"].value, 
                          show_vectors=view_settings["show_vectors"].value, 
                          show_prototile=view_settings["show_prototile"].value,
                          show_reg_prototile=view_settings["show_reg_prototile"].value,
                          show_ids=view_settings["show_ids"].value,
                          cmap=cm)
        plot.axis("off")
        return plot
    return (plot_tiles,)


if __name__ == "__main__":
    app.run()
