## \file ReducedSM.py
# This is a script which runs the reducedSM function from visualization_manager_DJ.py
#

from floris.floris import Floris
from copy import deepcopy
from visualization_manager_DJ import VisualizationManager
import json
import numpy as np
import matplotlib.pyplot as plt

with open("example_input.json") as WFJSON:
    ## a JSON windfarm object read from a file
    WF = json.load(WFJSON) # a JSON windfarm object read from a file
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
grid_resolution = [25, 25, 15] # [x_res, y_res, z_res]

## VisualizationManager object
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager

# the following are the parameters used by the reducedSM() function
vman.params.v0 = 8.0 # initial speed estimate
vman.params.d0 = np.deg2rad(0.0) # initial direction estimate
vman.params.epSpeed = 0.001  # speed epsilon (ev)
vman.params.epDir = 0.0001  # direction epsilon (ed)
vman.params.spErrMax = 0.1  # speed error threshold
vman.params.dirErrMax = 0.01  # direction error threshold
vman.params.iterMax = 100  # iteration threshold
vman.params.vBar = 7.5 # actual wind speed
vman.params.dBar = np.deg2rad(0.5) # actual wind direction
vman.params.mask_thresh = 0.7 # masking threshold for normalized sensitivity matrix

# run the reducedSM function
vman.reducedSM()