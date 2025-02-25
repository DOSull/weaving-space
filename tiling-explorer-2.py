import marimo

__generated_with = "0.11.0"
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
    import math
    return (math,)


@app.cell(hide_code=True)
def _():
    import geopandas as gpd
    import weavingspace.tile_unit as tu
    import weavingspace.tile_map as tm
    # from weavingspace.tile_unit import TileUnit
    # from weavingspace.tile_map import Tiling
    return (gpd,) #TileUnit, Tiling, gpd


@app.cell(hide_code=True)
def _(gpd):
    gdf = gpd.read_file("https://dosull.github.io/weaving-space/examples/data/dummy-data.json", engine="fiona")
    return (gdf,)


@app.cell(hide_code=True)
def _(Tiling, gdf, tile):
    tiled_map = tm.Tiling(tile, gdf).get_tiled_map()
    return (tiled_map,)


@app.cell(hide_code=True)
def _(gdf, mo, tile):
    _vars = [_ for _ in gdf.columns if not "geom" in _]
    _pals = ['viridis', 'magma', 'Greys', 'Purples', 'Blues', 
             'Greens', 'Oranges', 'Reds', 'YlOrBr', 'YlOrRd', 
             'OrRd', 'PuRd', 'RdPu', 'BuPu', 'GnBu', 
             'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    vars = mo.ui.array([mo.ui.dropdown(options=_vars, value=_vars[i], label=f"Tiles '{id}'") 
                        for i, id in enumerate(tile.tiles.tile_id)], label="Variables") 
    pals = mo.ui.array([mo.ui.dropdown(options=_pals, value="Reds") for i, id in enumerate(tile.tiles.tile_id)], 
                      label="Palettes")
    mo.md(f"#### {mo.hstack([vars, pals])}")
    return pals, vars


@app.cell(hide_code=True)
def _(pals, tile, tiled_map, vars):
    tiled_map.variables = {k: v for k, v in zip(tile.tiles.tile_id, vars.value)}
    tiled_map.colourmaps = {k: v for k, v in zip(vars.value, pals.value)}
    tiled_map.render(legend=False)
    return


@app.cell(hide_code=True)
def _(mo):
    t_inset = mo.ui.slider(steps=[0,1,2,3,5,10,20,30,50], value=0, show_value=True)
    p_inset = mo.ui.slider(steps=[0,1,2,3,5,10,20,30,50,100,250,500], value=0, show_value=True)
    tile_rotate = mo.ui.slider(start=-90, stop=90, step=5, value=0, show_value=True)
    aspect = mo.ui.slider(start=0.3, stop=3, step=0.01, value=1, show_value=True)
    crs = mo.ui.dropdown(options={"3857":3857, "4326":4326, "2193":2193}, value="3857")
    spacing = mo.ui.slider(start=50, stop=5000, step=50, value=750)

    mo.md("\n".join(["### Tiling settings",
      f"#### Spacing {spacing}&nbsp;&nbsp;CRS {crs}",               
      f"#### Rotate by {tile_rotate}&nbsp;&nbsp;Width/height {aspect}",               
      f"#### Tile inset {t_inset}&nbsp;&nbsp;Prototile inset {p_inset}"]))
    return aspect, crs, p_inset, spacing, t_inset, tile_rotate


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
        _n = [x for x in range(2, 21)]
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
    aspect,
    code,
    crs,
    math,
    n,
    offset,
    p_inset,
    spacing,
    t_inset,
    tile_rotate,
    tiling_type,
):
    tile = tu.TileUnit(tiling_type=tiling_type.value, n = n.value, code = code.value, offset = offset.value,
                   spacing=spacing.value, crs=crs.value) \
          .transform_rotate(tile_rotate.value) \
          .transform_scale(math.sqrt(aspect.value), 1/math.sqrt(aspect.value)) \
          .inset_tiles(t_inset.value) \
          .inset_prototile(p_inset.value)
    return (tile,)


@app.cell
def _(plot_tiles, tile):
    plot_tiles(tile)
    return


@app.cell(hide_code=True)
def _(
    aspect,
    code,
    crs,
    math,
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
    _str = _str + f', spacing={spacing.value}, crs={crs.value})'
    if tile_rotate.value != 0:
        _str = _str + f'.transform_rotate({tile_rotate.value})'
    if aspect.value != 1:
        _str = _str + f'.transform_scale({math.sqrt(aspect.value):.3f}, {1/math.sqrt(aspect.value):.3f})'
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
    palette,
    radius,
    show_ids,
    show_prototile,
    show_reg_prototile,
    show_vectors,
):
    def plot_tiles(tiles):
        plot = tiles.plot(r=radius.value, 
                          show_vectors=show_vectors.value, 
                          show_prototile=show_prototile.value,
                          show_reg_prototile=show_reg_prototile.value,
                          show_ids=show_ids.value,
                          cmap=palette.value)
        plot.axis("off")
        return plot
    return (plot_tiles,)


if __name__ == "__main__":
    app.run()
