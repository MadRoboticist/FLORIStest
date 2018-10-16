## \file PlotSensitivityMatrix.py
# This is a script which runs the plotSensitivityMatrix function from visualization_manager_DJ.py

from floris.floris import Floris
from visualization_manager_DJ import VisualizationManager
import json
import numpy as np

with open("example_input.json") as WFJSON:
    ## a JSON windfarm object read from a file
    WF = json.load(WFJSON) # a JSON windfarm object read from a file
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
grid_resolution = [25, 25, 15] # [x_res, y_res, z_res]
## VisualizationManager object
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager
# These are the parameters used by the plotSensitivityMatrix() function
vman.params.v0 = 8.0 # the initial wind estimate
vman.params.d0 = np.deg2rad(0.0) # the initial direction estimate
vman.params.epSpeed = 0.001  # speed epsilon (ev)
vman.params.epDir = 0.0001  # direction epsilon (ed)
# plot the sensitivity matrix
vman.plotSensitivityMatrix()
