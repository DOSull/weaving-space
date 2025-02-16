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
    palette = mo.ui.dropdown(options=["tab10", "tab20", "tab20b"], value="tab10")

    mo.md("\n".join(["### General settings",
      f"#### Set radius {radius}&nbsp;&nbsp;Tile inset {t_inset}&nbsp;&nbsp;Prototile inset {p_inset}",
      f"#### Show base tile {show_prototile}&nbsp;&nbsp;Show repeat unit {show_reg_prototile}&nbsp;Palette {palette}"]))
    return (
        p_inset,
        palette,
        radius,
        show_prototile,
        show_reg_prototile,
        t_inset,
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
def _(TileUnit, hex_or_square, n_cols, p_inset, plot_tiles, t_inset):
    colourings = TileUnit(tiling_type=hex_or_square.value, n = n_cols.value).inset_tiles(t_inset.value).inset_prototile(p_inset.value)
    plot_tiles(colourings)
    return (colourings,)


@app.cell
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
        ## Hexagon slices and dissections
        We can subdivide hexagons in many different ways. The most easily specified are 'pie slices' which, given the 6 sides of a hexagon yield regular tilings when 2, 3, 4, 6 or 12 slices are made. Dissections are more arbitrary (there are many more possibilities than implemented here).
        """
    )
    return


@app.cell(hide_code=True)
def _(mo):
    slice_or_dissect = mo.ui.dropdown(
        options=["hex-slice", "hex-dissection"],
        value="hex-slice"
    )
    mo.md(f"#### Slice or dissection? {slice_or_dissect}")
    return (slice_or_dissect,)


@app.cell(hide_code=True)
def _(mo, slice_or_dissect):
    n_pieces = mo.ui.dropdown(
        options=({"2": 2, "3": 3, "4": 4, "6": 6, "12": 12}
                 if slice_or_dissect.value == "hex-slice"
                 else {"4": 4, "7": 7, "9": 9}),
        value="6" if slice_or_dissect.value == "hex-slice" else "7"
    )
    mo.md(f"#### Choose number of pieces {n_pieces}")
    return (n_pieces,)


@app.cell(hide_code=True)
def _(mo):
    offset_cuts = mo.ui.switch(value=True)
    mo.md(f"#### Centre cuts on edges {offset_cuts}")
    return (offset_cuts,)


@app.cell
def _(
    TileUnit,
    n_pieces,
    offset_cuts,
    p_inset,
    plot_tiles,
    slice_or_dissect,
    t_inset,
):
    hex = TileUnit(tiling_type=slice_or_dissect.value, 
                   n=n_pieces.value, offset=1 if offset_cuts.value else 0) \
          .inset_tiles(t_inset.value) \
          .inset_prototile(p_inset.value)
    plot_tiles(hex)
    return (hex,)


@app.cell
def _(hex, mo, n_pieces, offset_cuts, slice_or_dissect):
    _download_button = mo.download(data=hex.tiles.to_json().encode('utf-8'), 
                filename=f'{slice_or_dissect.value}-{n_pieces.value}-{offset_cuts.value}_tile_unit.json', 
                mimetype='text/plain', label='Download')
    mo.md(f'Click to download **hex slice/dissection** tile unit {_download_button}')
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Laves tilings and their Archimedean duals""")
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
def _(TileUnit, laves_or_arch, p_inset, plot_tiles, t_inset, tiling_code):
    laves_or_arch_tiles = TileUnit(tiling_type=laves_or_arch.value,
                             code=tiling_code.value) \
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
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
