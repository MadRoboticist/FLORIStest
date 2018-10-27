## \file GreedyPathOnline.py
# This is a script which utilizes the greedyPath function in pathPlan.py
# while doing online updates of the sensitivity matrix and associated estimates

from visualization_manager_DJ import VisualizationManager
from pathPlan import PathPlanner
from UAV import UAV
import json
import numpy as np
import matplotlib.pyplot as plt

with open("example_input.json") as WFJSON:
    ## a JSON windfarm object read from a file
    WF = json.load(WFJSON) # a JSON windfarm object read from a file
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
grid_resolution = [25, 25, 15] # [x_res, y_res, z_res]

## @var vman
# @brief VisualizationManager object
#
# see definition at visualization_manager_DJ.VisualizationManager
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager

## The initial wind speed estimate
#
# see definition at visualization_manager_DJ.VisualizationManager.params.v0
vman.params.v0 = 8.0 # the initial wind speed estimate
## The initial wind direction estimate
#
# see definition at visualization_manager_DJ.VisualizationManager.params.d0
vman.params.d0 = np.deg2rad(0.0) # the initial direction estimate
## Speed epsilon (ev) which is used to calculate the sensitivity matrix
#
# see definition at visualization_manager_DJ.VisualizationManager.params.epSpeed
vman.params.epSpeed = 0.001  # speed epsilon (ev)
## Direction epsilon (ed) which is used to calculate the sensitivity matrix
#
# see definition at visualization_manager_DJ.VisualizationManager.params.epDir
vman.params.epDir = 0.0001  # direction epsilon (ed)
## Speed error minimum threshold
#
# see definition at visualization_manager_DJ.VisualizationManager.params.spErrMax
vman.params.spErrMax = 0.1  # speed error threshold
## Direction error minimum threshold
#
# see definition at visualization_manager_DJ.VisualizationManager.params.dirErrMax
vman.params.dirErrMax = 0.01  # direction error threshold
## Maximum number of iterations to complete
#
# see definition at visualization_manager_DJ.VisualizationManager.params.iterMax
vman.params.iterMax = 100  # iteration threshold
## Actual wind speed
#
# see definition at visualization_manager_DJ.VisualizationManager.params.vBar
vman.params.vBar = 7.5 # actual wind speed
## Actual wind direction in radians
#
# see definition at visualization_manager_DJ.VisualizationManager.params.dBar
vman.params.dBar = np.deg2rad(0.5) # actual wind direction

## A PathPlanner object
#
# see definition at pathPlan.PathPlanner
planner = PathPlanner(vman)

## A UAV object
#
# see definition at pathPlan.UAV
UAV1 = UAV(planner)
## Holds a UAV's "GPS" point
# starting point for this script is 200,-250
#
# see definition at pathPlan.UAV.GPS
UAV1.GPS = [200,-250]
#run the first iteration
UAV1.plan_horizon = 30
UAV1.moves2recalc = 10
UAV1.patrolMax = 125
planner.greedyPath(UAV1)

for i in range(50):
    # run the greedy path algorithm from pathPlan

    # plot the score map with the UAV's path on it
    UAV1.move()
    planner.updateEstimates(UAV1)
    planner.greedyPath(UAV1)
filename = "UAV"
planner.plotHistory(UAV1,filename)




