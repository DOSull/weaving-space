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
    from weavingspace import Tiling
    return TileUnit, Tiling, gpd, mpl


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
    _vars = [_ for _ in gdf.columns if not "geom" in _]
    _pals = ['Reds', 'summer', 'Oranges', 'viridis',
             'spring', 'Greens', 'YlGnBu', 'Blues',
             'winter', 'Purples', 'RdPu', 'Greys']
    vars = mo.ui.array([mo.ui.dropdown(options=_vars, value=_vars[i], label=f"Tiles '{id}'") 
                        for i, id in enumerate(tile.tiles.tile_id)], label="Variables") 
    pals = mo.ui.array([mo.ui.dropdown(options=_pals, value=_pals[i]) for i, id in enumerate(tile.tiles.tile_id)], 
                      label="Palettes")
    mo.md(f"#### {mo.hstack([vars, pals], align="center")}")
    return pals, vars


@app.cell(hide_code=True)
def _(mpl, pals):
    _n = len(pals)
    _fig, _axs = mpl.pyplot.subplots(nrows = _n, figsize=(2, 0.2 + 0.34 * _n))
    for ax, cm in zip(_axs, pals.value):
        _xy = [[x / 256 for x in range(257)], [x / 256 for x in range(257)]]
        ax.imshow(_xy, aspect='auto', cmap=mpl.colormaps.get(cm))
    for ax in _axs:
        ax.set_axis_off()
    ax
    return ax, cm


@app.cell
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
    tiled_map.variables = {k: v for k, v in zip(tile.tiles.tile_id, vars.value)}
    tiled_map.colourmaps = {k: v for k, v in zip(vars.value, pals.value)}
    tiled_map.render(legend=False)
    return


@app.cell(hide_code=True)
def _(mo):
    tile_rotate = mo.ui.slider(start=-90, stop=90, step=5, value=0, show_value=True, debounce=True)
    spacing = mo.ui.slider(start=50, stop=5000, step=50, value=750, debounce=True)

    mo.md("\n".join([
        "### Tiling settings",
        f"#### Spacing {spacing}",
        f"#### Rotate by {tile_rotate}"
    ]))
    return spacing, tile_rotate


@app.cell
def _(mo, spacing):
    _max_prototile_spacing = spacing.value // 6
    p_inset = mo.ui.slider(start=0, stop=_max_prototile_spacing, step = 1, 
                           value=0, show_value=True, debounce=True)
    mo.md(f"#### Prototile inset {p_inset}")
    return (p_inset,)


@app.cell
def _(mo, p_inset):
    _max_tile_spacing = p_inset.value // 5
    t_inset = mo.ui.slider(start=0, stop=_max_tile_spacing, step = 1, 
                           value=0, show_value=True, debounce=True)
    mo.md(f"#### Tile inset {t_inset}")
    return (t_inset,)


@app.cell(hide_code=True)
def _(mo):
    tiling_type = mo.ui.dropdown(
        options = ["hex-slice", "square-slice", "hex-dissection", "cross",
                   "laves", "archimedean", "square-colouring", "hex-colouring"],
        value = "laves"
    )
    mo.md(f"#### Type (`tiling_type`): {tiling_type}")
    return (tiling_type,)


@app.cell(hide_code=True)
def _(mo, tiling_type):
    if "slice" in tiling_type.value:
        _n = [x for x in range(2, 13)]
        _v = 6
        _prefix = "Number of slices: "
    elif tiling_type.value == "hex-dissection":
        _n = [4, 7, 9]
        _v = 7
        _prefix = "Number of pieces: "
    elif tiling_type.value == "cross":
        _n = [x for x in range(2, 8)]
        _v = 5
        _prefix = "Number of crosses: "
    elif "colour" in tiling_type.value:
        _n = [x for x in range(2, 10)]
        _v = 7
        _prefix = "Number of colours: "
    else:
        _n = [0, 1]
        _prefix = "NA"
        _v = 0

    n = mo.ui.slider(steps = _n, value = _v, show_value=True)

    if "slice" in tiling_type.value:
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

    if tiling_type.value in ["laves", "archimedean"]:
        _str = f"#### Code (`code`): {code}"
    else:
        _str = f"#### {_prefix} (`n`) {n}"
        if "slice" in tiling_type.value or tiling_type.value == "hex-dissection":
            _str = "\n".join([_str, f"#### Offset (`offset`): {offset}"])

    mo.md(_str)
    return code, n, offset


@app.cell(hide_code=True)
def _(
    TileUnit,
    code,
    gdf,
    n,
    offset,
    p_inset,
    spacing,
    t_inset,
    tile_rotate,
    tiling_type,
):
    tile = TileUnit(tiling_type=tiling_type.value, n = n.value, code = code.value, offset = offset.value,
                   spacing=spacing.value, crs=gdf.crs) \
          .transform_rotate(tile_rotate.value) \
          .inset_tiles(t_inset.value) \
          .inset_prototile(p_inset.value)
    return (tile,)


@app.cell
def _(plot_tiles, tile):
    plot_tiles(tile)
    return


@app.cell(hide_code=True)
def _(
    code,
    gdf,
    mo,
    n,
    offset,
    p_inset,
    spacing,
    t_inset,
    tile_rotate,
    tiling_type,
):
    _str = f'### `TileUnit(tiling_type="{tiling_type.value}"'
    if tiling_type.value in ["laves", "archimedean"]:
        _str = _str + f', code="{code.value}"'
    else:
        _str = _str + f', n={n.value}'
        if tiling_type.value in ["hex-slice", "square-slice", "hex-dissection"]:
            _str = _str + f', offset={offset.value}'
    _str = _str + f', spacing={spacing.value}, crs={gdf.crs})'
    if tile_rotate.value != 0:
        _str = _str + f'.transform_rotate({tile_rotate.value})'
    if t_inset.value != 0:
        _str = _str + f'.inset_tiles({t_inset.value})'
    if p_inset.value != 0:
        _str = _str + f'.inset_prototile({p_inset.value})'
    _str = _str + f'`'

    mo.md(_str)
    return


@app.cell(hide_code=True)
def _(mo):
    radius = mo.ui.slider(0, 4, value=1, show_value=True)
    show_prototile = mo.ui.switch(value=False)
    show_reg_prototile = mo.ui.switch(value=False)
    show_vectors = mo.ui.switch(value=False)
    show_ids = mo.ui.switch(value=True)
    palette = mo.ui.dropdown(options=["Spectral", "tab10", "tab20"], value="Spectral")

    mo.md("\n".join(["### Tile design view settings",
      f"#### Set radius {radius}",
      f"#### Show vectors {show_vectors}",
      f"#### Show base tile {show_prototile}",
      f"#### Show prototile {show_reg_prototile}",
      f"#### Show tile IDs {show_ids}",
      f"#### Palette {palette}"]))
    return (
        palette,
        radius,
        show_ids,
        show_prototile,
        show_reg_prototile,
        show_vectors,
    )


@app.cell(hide_code=True)
def _(
    mpl,
    pals,
    radius,
    show_ids,
    show_prototile,
    show_reg_prototile,
    show_vectors,
):
    def plot_tiles(tiles):
        cm = mpl.colors.ListedColormap([mpl.colormaps.get(p)(0.75) for p in pals.value])
        plot = tiles.plot(r=radius.value, 
                          show_vectors=show_vectors.value, 
                          show_prototile=show_prototile.value,
                          show_reg_prototile=show_reg_prototile.value,
                          show_ids=show_ids.value,
                          cmap=cm)
        plot.axis("off")
        return plot
    return (plot_tiles,)


if __name__ == "__main__":
    app.run()
