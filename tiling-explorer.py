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
    radius = mo.ui.slider(0, 4, value = 1)
    t_inset = mo.ui.slider(0, 20, 1, value = 0)
    p_inset = mo.ui.slider(0, 50, 5, value = 10)
    show_prototile = mo.ui.switch(value=False)
    show_reg_prototile = mo.ui.switch(value=False)

    mo.md("\n".join(["### General settings",
                     f"#### Set radius {radius}&nbsp;&nbsp;Tile inset {t_inset}&nbsp;&nbsp;Prototile inset {p_inset}",
                     f"#### Show base tile {show_prototile}&nbsp;&nbsp;Show repeat unit {show_reg_prototile}"]))
    return p_inset, radius, show_prototile, show_reg_prototile, t_inset


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Regular tiling colourings
        These are simple regular colourings of hexagonal or square tilings. 3, 4, and 7 colourings of hexagonal tiles are available, but only a 5-colouring of squares. More complex rectangular colourings can be achieved with weave units.
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
def _(hex_or_square, mo):
    n_cols = mo.ui.dropdown(
        options={"3":3, "4":4, "7":7} if hex_or_square.value == "hex-colouring" else {"5":5},
        value="3" if hex_or_square.value == "hex-colouring" else "5"
    )
    mo.md(f"#### Number of colours {n_cols}")
    return (n_cols,)


@app.cell(hide_code=True)
def _(TileUnit, hex_or_square, n_cols, p_inset, plot_tiles, t_inset):
    colourings = TileUnit(tiling_type=hex_or_square.value, n = n_cols.value).inset_tiles(t_inset.value).inset_prototile(p_inset.value)
    plot_tiles(colourings)
    return (colourings,)


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
    plot_tiles(TileUnit(tiling_type=slice_or_dissect.value, n=n_pieces.value, offset=1 if offset_cuts.value else 0).inset_tiles(t_inset.value).inset_prototile(p_inset.value))
    return


@app.cell(hide_code=True)
def _(radius, show_prototile, show_reg_prototile):
    def plot_tiles(tiles):
        plot = tiles.plot(r=radius.value, 
                          show_vectors=False, 
                          show_prototile=show_prototile.value,
                          show_reg_prototile=show_reg_prototile.value,
                          show_ids=False)
        plot.axis("off")
        return plot
    return (plot_tiles,)


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


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
    plot_tiles(TileUnit(tiling_type=laves_or_arch.value, code=tiling_code.value).inset_tiles(t_inset.value).inset_prototile(p_inset.value))
    return


if __name__ == "__main__":
    app.run()
