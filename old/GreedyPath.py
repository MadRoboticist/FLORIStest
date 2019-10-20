## \file GreedyPath.py
# This is a script which utilizes the greedyPath function in pathPlan.py
#

from old.visualization_manager_DJ import VisualizationManager
from old.pathPlan import PathPlanner
from UAV import UAV
import json
import numpy as np

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

## A PathPlanner object
#
# see definition at pathPlan.PathPlanner
planner = PathPlanner(vman)

## A UAV object
#
# see definition at UAV.UAV
UAV1 = UAV(planner)
## Holds a UAV's GPS point
#
# see definition at UAV.UAV.GPS
UAV1.GPS = [200,-250]
## Holds the UAV's planning horizon
#  (the number of steps to plan ahead)
#
# see definition at UAV.UAV.plan_horizon
UAV1.plan_horizon = 150
# run the greedy path algorithm from pathPlan for 150 iterations
planner.greedyPath(UAV1)
# plot the score map with the UAV's path on it
planner.plotScoreMapUAV(UAV1)
# show the plot
planner.show()


