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
## The initial speed estimate
vman.params.v0 = 8.0 # initial speed estimate
## The initial direction estimate
vman.params.d0 = np.deg2rad(0.0) # initial direction estimate
## speed epsilon (ev)
vman.params.epSpeed = 0.001  # speed epsilon (ev)
## direction epsilon (ed)
vman.params.epDir = 0.0001  # direction epsilon (ed)
## speed error minimum threshold
vman.params.spErrMax = 0.1  # speed error threshold
## direction error minimum threshold
vman.params.dirErrMax = 0.01  # direction error threshold
## maximum number of iterations to be completed
vman.params.iterMax = 100  # iteration threshold
## The actual wind speed
vman.params.vBar = 7.5 # actual wind speed
## The actual wind direction
vman.params.dBar = np.deg2rad(0.5) # actual wind direction
## A masking threshold for reducing the sensitivity matrix
vman.params.mask_thresh = 0.7 # masking threshold for reducing normalized sensitivity matrix

# run the reducedSM function
vman.reducedSM()