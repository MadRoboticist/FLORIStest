## \file AnimateError.py
# This is a script which runs the animateDnSerror function from visualization_manager_DJ.py

## \file example_input.json
# This is an example JSON windfarm input file for the Floris model.

from FLORIS_iface import FLORIS_sub as Floris
from VIZMAN_new import VisualizationManager
import numpy as np

grid_resolution = [48, 15, 15]
FL = Floris("twoTurb_input_new.json", grid_resolution)
FL.set_incoming(8.0,270.0)



vman = VisualizationManager(FL, grid_resolution)  # a new VisualizationManager object

# these are the parameters used by the animateDnSerror() function
## {min speed, max speed, step}
vman.params.Srange = [FL.wind_speed-1.0, FL.wind_speed+1.0, 1]
## {min dir error, max dir error, step}
vman.params.Drange = [FL.wind_direction-1.0, FL.wind_direction+1.0, 1]

# call the animateDnSerror function
vman.animateDnSerror()
