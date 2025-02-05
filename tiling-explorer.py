import marimo

__generated_with = "0.10.19"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    from weavingspace import TileShape
    from weavingspace import TileUnit

    mo.md("# All the tiles!")
    return TileShape, TileUnit, mo


@app.cell
def _(mo):
    mo.md(r"""## Base tile units""")
    return


@app.cell(hide_code=True)
def _(TileShape, mo):
    base_shapes = mo.ui.dropdown(
        options=dict(zip(
            ["RECTANGLE", "HEXAGON", "TRIANGLE", "DIAMOND"],
            list(TileShape)
        )),
        label="Choose base tiling shape",
        value="RECTANGLE"
    )
    base_shapes
    return (base_shapes,)


@app.cell
def _(TileUnit, base_shapes):
    _unit = TileUnit(base_shape=base_shapes.value)
    _plot = _unit.plot(r=1, show_vectors=True, show_reg_prototile=False)
    _plot.axis('off')
    _plot
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## Hexagon colourings""")
    return


@app.cell
def _(mo):
    _n = [3, 4, 7]
    _n_str = [str(x) for x in _n]
    n_colours = mo.ui.dropdown(
        options=dict(zip(_n_str, _n)),
        label="Choose how many colours",
        value="3"
    )
    n_colours
    return (n_colours,)


@app.cell
def _(TileUnit, n_colours):
    _unit = TileUnit(tiling_type="hex-colouring", n=n_colours.value)
    _plot = _unit.plot(r=1, show_vectors=True, show_reg_prototile=False)
    _plot.axis("off")
    _plot
    return


if __name__ == "__main__":
    app.run()
