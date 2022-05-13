root = "/home/osullid3/Documents/code/weaving-space/"
code_path = root + "weavingspace/"
data_path = root + "data/"

import sys
sys.path.append(code_path)

import geopandas

from weave_units import WeaveUnit
from tile_map import Tiling

ak = geopandas.read_file(data_path + "/imd-auckland-2018.gpkg")
w1 = WeaveUnit(weave_type = "twill", spacing = 200,  aspect = 0.75,
                    margin = .03, strands = "abc|defg", crs = 2193)

weave1 = Tiling(w1, ak, id_var = "DZ2018")
textile1 = weave1.get_tiled_map(rotation = 30,  prioritise_tiles = True)
textile1.to_file(data_path + "output.gpkg")

vl = QgsVectorLayer(data_path + "/output.gpkg")
QgsProject.instance().addMapLayer(vl)