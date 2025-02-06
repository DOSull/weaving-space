import marimo

__generated_with = "0.11.0"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""# Tilings explorer""")
    return


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _():
    from weavingspace.tile_unit import TileUnit
    from weavingspace.tile_unit import TileShape
    return TileShape, TileUnit


@app.cell(hide_code=True)
def _(mo):
    radius = mo.ui.slider(0, 4, label = "Choose radius `r`, for tiling displays", value = 1)
    radius
    return (radius,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Hexagon colourings""")
    return


@app.cell(hide_code=True)
def _(mo):
    n_hex_cols = mo.ui.dropdown(
        options = {"3": 3, "4": 4, "7": 7},
        label = "Choose number of colours `n`",
        value = "3"
    )
    n_hex_cols
    return (n_hex_cols,)


@app.cell(hide_code=True)
def _(TileUnit, n_hex_cols, radius):
    _plot = TileUnit(tiling_type="hex-colouring", n = n_hex_cols.value).plot(
        r=radius.value, show_vectors=True, show_reg_prototile=False)
    _plot.axis("off")
    _plot
    return


if __name__ == "__main__":
    app.run()
