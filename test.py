import marimo

__generated_with = "0.11.0"
app = marimo.App(width="medium")


@app.cell
def _():
    return


@app.cell
def _():
    import geopandas as gpd
    return (gpd,)


@app.cell
def _(gpd):
    gpd.read_file("https://dosull.github.io/weaving-space/examples/data/dummy-data.json").shape
    return


if __name__ == "__main__":
    app.run()
