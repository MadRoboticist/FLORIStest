## \file PlotSensitivityMatrix.py
# This is a script which runs the plotSensitivityMatrix function from visualization_manager_DJ.py

from old.visualization_manager_DJ import VisualizationManager
import json
import numpy as np

with open("example_input.json") as WFJSON:
    ## a JSON windfarm object read from a file
    WF = json.load(WFJSON) # a JSON windfarm object read from a file
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
grid_resolution = [48, 15, 15] # [x_res, y_res, z_res]
## VisualizationManager object
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager
# These are the parameters used by the plotSensitivityMatrix() function
## The initial wind estimate
vman.params.v0 = 13.0 # the initial wind estimate
## The initial direction estimate
vman.params.d0 = np.deg2rad(10.0) # the initial direction estimate
## speed epsilon (ev)
vman.params.epSpeed = 0.001  # speed epsilon (ev)
## direction epsilon (ed)
vman.params.epDir = 0.0001  # direction epsilon (ed)
# plot the sensitivity matrix
vman.plotSensitivityMatrix()
