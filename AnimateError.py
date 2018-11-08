## \file AnimateError.py
# This is a script which runs the animateDnSerror function from visualization_manager_DJ.py

## \file example_input.json
# This is an example JSON windfarm input file for the Floris model.

from floris.floris import Floris
from copy import deepcopy
from visualization_manager_DJ import VisualizationManager
import json

with open("example_input.json") as WFJSON:
    ## a JSON windfarm object read from a file
    WF = json.load(WFJSON) # a JSON windfarm object read from a file
## grid resolution: [x_resolution, y_resolution, z_resolution]
grid_resolution = [25, 25, 15]  # grid resolution: [x_resolution, y_resolution, z_resolution]
## a new VisualizationManager object
vman = VisualizationManager(WF, grid_resolution)  # a new VisualizationManager object

# these are the parameters used by the animateDnSerror() function
## average wind speed
vman.params.speed = WF['farm']['properties']['wind_speed']
## {min speed, max speed, step}
vman.params.Srange = [vman.params.speed-5.0, vman.params.speed+5.0, 1]
## {min dir error, max dir error, step}
vman.params.Drange = [-10.0, 10.0, 1]

# call the animateDnSerror function
vman.animateDnSerror()
