import marimo

__generated_with = "0.11.9"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    import io
    import json
    import geopandas
    import shapely.geometry as geom
    return geom, geopandas, io, json


@app.cell
def _(mo):
    f = mo.ui.file()
    f
    return (f,)


@app.cell
def _(geom):
    def get_shape(coords, shape_type):
        if shape_type == "Polygon":
            return geom.Polygon(coords)
        elif shape_type == "MultiPolygon":
            return geom.MultiPolygon(coords)
        return None
    return (get_shape,)


@app.cell
def _(f, geopandas, get_shape, io, json):
    _data = None if len(f.value) == 0 else f.value[0].contents
    _geojson_dict = json.loads(io.BytesIO(_data).read())

    _crs = _geojson_dict["crs"]["properties"]["name"]
    _features = _geojson_dict["features"]
    _props = [p["properties"] for p in _features]
    _geoms = [get_shape(p["geometry"]["coordinates"], p["geometry"]["type"]) for p in _features]

    gdf = geopandas.GeoDataFrame(data=_props, geometry=_geoms, crs=_crs).to_crs(3857)
    return (gdf,)


@app.cell
def _(gdf):
    gdf.plot()
    return


if __name__ == "__main__":
    app.run()
