import marimo

__generated_with = "0.11.0"
app = marimo.App(width="medium")


@app.cell
def _():
    from weavingspace.tile_unit import TileUnit
    return (TileUnit,)


@app.cell(hide_code=True)
def _(mo):
    radius = mo.ui.slider(0, 4, value=1, show_value=True)
    t_inset = mo.ui.slider(0, 20, 1, value=0, show_value=True)
    p_inset = mo.ui.slider(0, 50, 5, value=10, show_value=True)
    show_prototile = mo.ui.switch(value=False)
    show_reg_prototile = mo.ui.switch(value=False)
    show_vectors = mo.ui.switch(value=False)
    tile_rotate = mo.ui.slider(steps = [x for x in range(-90, 91)], value=0, show_value=True)
    aspect = mo.ui.slider(start=0.3, stop=3, step=0.01, value=1, show_value=True)
    palette = mo.ui.dropdown(options=["Spectral", "tab10", "tab20"], value="Spectral")

    mo.md("\n".join(["### General settings",
      f"#### Rotate tile unit by {tile_rotate}&nbsp;&nbsp;Width/height {aspect}&nbsp;&nbsp;Set radius {radius}",               
      f"#### Tile inset {t_inset}&nbsp;&nbsp;Prototile inset {p_inset}",
      f"#### Show vectors {show_vectors}&nbsp;&nbsp;Show base tile {show_prototile}&nbsp;&nbsp;Show repeat unit {show_reg_prototile}&nbsp;&nbsp;Palette {palette}"]))
    return (
        aspect,
        p_inset,
        palette,
        radius,
        show_prototile,
        show_reg_prototile,
        show_vectors,
        t_inset,
        tile_rotate,
    )


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
        _v = 12
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
    math,
    n,
    offset,
    p_inset,
    plot_tiles,
    t_inset,
    tile_rotate,
    tiling_type,
):
    tile = TileUnit(tiling_type=tiling_type.value, n = n.value, code = code.value, offset = offset.value) \
          .transform_rotate(tile_rotate.value) \
          .transform_scale(math.sqrt(aspect.value), 1/math.sqrt(aspect.value)) \
          .inset_tiles(t_inset.value) \
          .inset_prototile(p_inset.value)
    plot_tiles(tile)
    return (tile,)


@app.cell(hide_code=True)
def _(palette, radius, show_prototile, show_reg_prototile, show_vectors):
    def plot_tiles(tiles):
        plot = tiles.plot(r=radius.value, 
                          show_vectors=show_vectors.value, 
                          show_prototile=show_prototile.value,
                          show_reg_prototile=show_reg_prototile.value,
                          show_ids=False,
                          cmap=palette.value)
        plot.axis("off")
        return plot
    return (plot_tiles,)


@app.cell
def _():
    import math
    return (math,)


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
