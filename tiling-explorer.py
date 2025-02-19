import marimo

__generated_with = "0.11.0"
app = marimo.App(
    width="medium",
    layout_file="layouts/tiling-explorer.grid.json",
)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Tilings explorer""")
    return


@app.cell(hide_code=True)
def _():
    from weavingspace.tile_unit import TileUnit
    from weavingspace.tile_unit import TileShape
    return TileShape, TileUnit


@app.cell(hide_code=True)
def _(mo):
    radius = mo.ui.slider(0, 4, value=1)
    t_inset = mo.ui.slider(0, 20, 1, value=0)
    p_inset = mo.ui.slider(0, 50, 5, value=10)
    show_prototile = mo.ui.switch(value=False)
    show_reg_prototile = mo.ui.switch(value=False)
    tile_rotate = mo.ui.slider(start=0, stop=90, step=5, value=0)
    aspect = mo.ui.slider(start=0.3, stop=3, step=0.01, value=1)
    palette = mo.ui.dropdown(options=["Spectral", "tab10", "tab20"], value="Spectral")

    mo.md("\n".join(["### General settings",
      f"#### Rotate tile unit by {tile_rotate}&nbsp;&nbsp;Width/height {aspect}",               
      f"#### Tile inset {t_inset}&nbsp;&nbsp;Prototile inset {p_inset}",
      f"#### Set radius {radius}&nbsp;&nbsp;Show base tile {show_prototile}&nbsp;&nbsp;Show repeat unit {show_reg_prototile}&nbsp;Palette {palette}"]))
    return (
        aspect,
        p_inset,
        palette,
        radius,
        show_prototile,
        show_reg_prototile,
        t_inset,
        tile_rotate,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Regular tiling colourings
        Regular colourings of hexagonal or square tilings, with any number from 2 to 9 (inclusive) of colours allowed. Potentially more complex square colourings can be achieved with weave patterns.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    hex_or_square = mo.ui.dropdown(
        options=["hex-colouring", "square-colouring"], 
        value="hex-colouring")
    mo.md(f"#### Hex or square? {hex_or_square}")
    return (hex_or_square,)


@app.cell(hide_code=True)
def _(list_to_dict, mo):
    n_cols = mo.ui.dropdown(
        options=list_to_dict(range(2, 10)),
        value="7"
    )
    mo.md(f"#### Number of colours {n_cols}")
    return (n_cols,)


@app.cell(hide_code=True)
def _(
    TileUnit,
    aspect,
    hex_or_square,
    math,
    n_cols,
    p_inset,
    plot_tiles,
    t_inset,
    tile_rotate,
):
    colourings = TileUnit(tiling_type=hex_or_square.value, n = n_cols.value) \
      .transform_rotate(tile_rotate.value) \
      .transform_scale(math.sqrt(aspect.value), 1/math.sqrt(aspect.value)) \
      .inset_tiles(t_inset.value) \
      .inset_prototile(p_inset.value)
    plot_tiles(colourings)
    return (colourings,)


@app.cell(hide_code=True)
def _(colourings, hex_or_square, mo, n_cols):
    _download_button = mo.download(data=colourings.tiles.to_json().encode('utf-8'), 
                                  filename=f'{hex_or_square.value}-{n_cols.value}_tile_unit.json', 
                                  mimetype='text/plain', 
                                  label='Download')
    mo.md(f'Click to download **colourings** tile unit {_download_button}')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Hexagon and square slices
        Pie slices of hexagons or squares, with an offset specifying how far along from start to midpoint of the first segment slices should start.
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    hex_or_square_slice = mo.ui.dropdown(
        options=["hex-slice", "square-slice"],
        value="hex-slice"
    )
    mo.md(f"#### Hex or square? {hex_or_square_slice}")
    return (hex_or_square_slice,)


@app.cell(hide_code=True)
def _(mo):
    n_slices = mo.ui.slider(start=2, stop=24, step=1, value=6)
    mo.md(f"#### Set number of slices {n_slices}")
    return (n_slices,)


@app.cell
def _(mo):
    offset_slices = mo.ui.slider(steps=[0, 1/4, 1/3, 1/2, 2/3, 3/4, 1], 
                                 value=0, show_value=True)
    mo.md(f"#### Slice offset {offset_slices}")
    return (offset_slices,)


@app.cell
def _(
    TileUnit,
    aspect,
    hex_or_square_slice,
    math,
    n_slices,
    offset_slices,
    p_inset,
    plot_tiles,
    t_inset,
    tile_rotate,
):
    slices = TileUnit(tiling_type=hex_or_square_slice.value, 
                   n=n_slices.value, offset=offset_slices.value) \
      .transform_rotate(tile_rotate.value) \
      .transform_scale(math.sqrt(aspect.value), 1/math.sqrt(aspect.value)) \
      .inset_tiles(t_inset.value) \
      .inset_prototile(p_inset.value)
    plot_tiles(slices)
    return (slices,)


@app.cell
def _(hex_or_square_slice, mo, n_slices, offset_slices, slices):
    _download_button = mo.download(data=slices.tiles.to_json().encode('utf-8'), 
                filename=f'{hex_or_square_slice.value}-{n_slices.value}-{offset_slices.value:.3f}_tile_unit.json', 
                mimetype='text/plain', label='Download')
    mo.md(f'Click to download **square/hex slice** tile unit {_download_button}')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        ## Hex dissections
        4, 7 and 9 dissections of the hexagon in two different varieties.
        """
    )
    return


@app.cell(hide_code=True)
def _(list_to_dict, mo):
    n_pieces = mo.ui.dropdown(
        options = list_to_dict([4, 7, 9]),
        value = "7"
    )
    mo.md(f"#### Number of pieces&nbsp;&nbsp;{n_pieces}")
    return (n_pieces,)


@app.cell(hide_code=True)
def _(mo):
    dissection_offset = mo.ui.switch(value = True)
    mo.md(f"#### Start from corners?&nbsp;&nbsp;{dissection_offset}")
    return (dissection_offset,)


@app.cell(hide_code=True)
def _(
    TileUnit,
    aspect,
    dissection_offset,
    math,
    n_pieces,
    p_inset,
    plot_tiles,
    t_inset,
    tile_rotate,
):
    pieces = TileUnit(tiling_type="hex-dissection", 
                   n=n_pieces.value, offset=0 if dissection_offset.value else 1) \
      .transform_rotate(tile_rotate.value) \
      .transform_scale(math.sqrt(aspect.value), 1/math.sqrt(aspect.value)) \
      .inset_tiles(t_inset.value) \
      .inset_prototile(p_inset.value)
    plot_tiles(pieces)
    return (pieces,)


@app.cell
def _(dissection_offset, mo, n_pieces, pieces):
    _download_button = mo.download(data=pieces.tiles.to_json().encode('utf-8'), 
                filename=f'hex-dissection-{n_pieces.value}-{dissection_offset.value}_tile_unit.json', 
                mimetype='text/plain', label='Download')
    mo.md(f'Click to download **Hex dissection** tile unit {_download_button}')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Laves tilings and their Archimedean duals
        The Laves tilings are monohedral. They are dual to the Archimedean tilings by regular polygons. We support 7 of the 11 possible tilings.
        """
    )
    return


@app.cell
def _(mo):
    laves_or_arch = mo.ui.dropdown(options=["laves", "archimedean"], value="laves")
    mo.md(f"#### Laves or Archimedean? {laves_or_arch}")
    return (laves_or_arch,)


@app.cell(hide_code=True)
def _(mo):
    supported_tilings = ["3.3.3.3.6", "3.3.4.3.4", "3.4.6.4", "3.6.3.6", "3.12.12", "4.6.12", "4.8.8"]
    tiling_code = mo.ui.dropdown(options=supported_tilings, value="3.3.4.3.4")
    mo.md(f"#### Tiling code {tiling_code}")
    return supported_tilings, tiling_code


@app.cell(hide_code=True)
def _(
    TileUnit,
    aspect,
    laves_or_arch,
    math,
    p_inset,
    plot_tiles,
    t_inset,
    tile_rotate,
    tiling_code,
):
    laves_or_arch_tiles = TileUnit(tiling_type=laves_or_arch.value,
                             code=tiling_code.value) \
      .transform_rotate(tile_rotate.value) \
      .transform_scale(math.sqrt(aspect.value), 1/math.sqrt(aspect.value)) \
      .inset_tiles(t_inset.value) \
      .inset_prototile(p_inset.value)
    plot_tiles(laves_or_arch_tiles)
    return (laves_or_arch_tiles,)


@app.cell
def _(laves_or_arch, laves_or_arch_tiles, mo, tiling_code):
    _download_button = mo.download(data=laves_or_arch_tiles.tiles.to_json().encode('utf-8'), 
                filename=f'{laves_or_arch.value}-{tiling_code.value}_tile_unit.json', 
                mimetype='text/plain', label='Download')
    mo.md(f'Click to download **Laves/Archimedean** tile unit {_download_button}')
    return


@app.cell(hide_code=True)
def _(palette, radius, show_prototile, show_reg_prototile):
    def list_to_dict(lst):
        return dict(zip([str(v) for v in lst], lst))

    def plot_tiles(tiles):
        plot = tiles.plot(r=radius.value, 
                          show_vectors=False, 
                          show_prototile=show_prototile.value,
                          show_reg_prototile=show_reg_prototile.value,
                          show_ids=False,
                          cmap=palette.value)
        plot.axis("off")
        return plot
    return list_to_dict, plot_tiles


@app.cell(hide_code=True)
def _():
    import math
    return (math,)


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
