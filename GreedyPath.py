## \file GreedyPath.py
# This is a script which utilizes the greedyPath function in pathPlan.py
#

from visualization_manager_DJ import VisualizationManager
from pathPlan import UAV, dir, PathPlanner
import json
import numpy as np

with open("example_input.json") as WFJSON:
    ## a JSON windfarm object read from a file
    WF = json.load(WFJSON) # a JSON windfarm object read from a file
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
grid_resolution = [25, 25, 15] # [x_res, y_res, z_res]

## VisualizationManager object
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager
vman.params.v0 = 8.0 # the initial wind estimate
vman.params.d0 = np.deg2rad(0.0) # the initial direction estimate
vman.params.epSpeed = 0.001  # speed epsilon (ev)
vman.params.epDir = 0.0001  # direction epsilon (ed)
vman.params.spErrMax = 0.1  # speed error threshold
vman.params.dirErrMax = 0.01  # direction error threshold
vman.params.iterMax = 100  # iteration threshold
vman.params.vBar = 7.5 # actual wind speed
vman.params.dBar = np.deg2rad(0.5) # actual wind direction
vman.params.mask_thresh = 0.7 # masking threshold for reducing normalized sensitivity matrix

planner = PathPlanner(vman)
UAV1 = UAV()
UAV1.GPS = [200,-250]
planner.greedyPath(UAV1, 100)
planner.plotScoreMapUAV(UAV1)
planner.show()
#UAV1.GPS = [200,-250]
#planner.pathGreedy(UAV1, 150)
#planner.plotScoreMapUAV(UAV1)
#UAV1.GPS = [200,-250]
#planner.pathGreedy(UAV1, 200)
#planner.plotScoreMapUAV(UAV1)


