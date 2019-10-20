## \file COSPlanOnline.py
# This is a script which utilizes the greedyPath function in pathPlan.py
# while doing online updates of the sensitivity matrix and associated estimates

from visualization_manager_DJ import VisualizationManager
from pathPlan import PathPlanner
from readVTK import VTKreader
from UAV import UAV
import json
import numpy as np
import scipy.io as scio
from copy import deepcopy

with open("twoTurb_input.json") as WFJSON:
    ## a JSON windfarm object read from a file
    WF = json.load(WFJSON) # a JSON windfarm object read from a file
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
## grid_resolution = [60, 52, 15] # [x_res, y_res, z_res] ## 36 turbines
grid_resolution = [48, 15, 15] ## 2 turbines
## grid_resolution = [30, 15, 15] ## 1 turbine
grid_size = grid_resolution[0]*grid_resolution[1]
coverage = 1.0
iterations = 10

print("Initializing UAV")
WF['farm']['properties']['wind_speed']=13.0 # actual wind speed
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager
vman.params.vBar = 13.0
vman.params.dBar = np.deg2rad(10.0)
vman.params.spErrMax = -0.1  # speed error threshold
print("speed error threshold: "+str(vman.params.dirErrMax))
## direction error minimum threshold
vman.params.dirErrMax = np.deg2rad(0.1)  # direction error threshold
print("direction error threshold: "+str(vman.params.dirErrMax))

## Initialize UAV
planner = PathPlanner(vman)
uav = UAV(planner)
uav.GPS = [2000,800]
uav.minX = 0
uav.maxX = grid_resolution[0]
uav.minY = 0
uav.maxY = grid_resolution[1]
uav.d0 = np.deg2rad(0.0)
planner.error[0][1]=planner.params.dBar-uav.d0
uav.v0 =8.0
planner.error[0][0]=planner.params.vBar-uav.v0

## A parameter which decides how many moves the UAV will hold in memory.
# this is basically the number of nodes in the pseudoinverse mask
#
# see definition at pathPlan.UAV.patrolMax
uav.patrolMax = int(grid_size*coverage)

#print("patrol max: " + str(UAV1.patrolMax))
## A parameter which decides how far ahead the planner will work
# for the first iteration.
#
# see definition at pathPlan.UAV.plan_horizon
uav.init_plan_horizon = uav.patrolMax

## A parameter which decides how many moves the UAV will take on
# the first planned path before the first recalculation
#
# see definition at pathPlan.UAV.moves2recalc
uav.init_mask_size = int(uav.patrolMax*0.9)

## A parameter which decides how many moves the UAV will take on
# the planned path before it recalculates the estimates and the plan
#
# see definition at pathPlan.UAV.moves2recalc
#uav.moves2recalc = int(np.amin([uav.plan_horizon*0.5,15]))
uav.moves2recalc = int(uav.patrolMax*0.15)

## A parameter which decides how far ahead the planner will work
# after the first iteration.
#
# see definition at pathPlan.UAV.plan_horizon
#uav.plan_horizon = int(np.amin([uav.patrolMax*0.9,30]))
uav.plan_horizon = int(uav.patrolMax*0.3)






planner.percent_plan = 1 #percentage of planned steps that will be taken
# run it for the indicated number of recalculations
i=0

#while dirErr[i]>planner1.params.dirErrMax or spdErr[i]>planner.params.spErrMax:
for i in range(iterations):
    print('\n\n\niteration ' + str(i))
    #uav.planner.COSPlan(uav)  # recalculate the plan
    uav.planner.greedyPath(uav)
    uav.move()  # move the UAV

    if (len(uav.GPSpath) >= uav.init_mask_size):
        uav.planner.updateEstimatesRM(uav)  # recalculate the estimates


        uav.planner.error.append([uav.planner.params.vBar - uav.v0,
                           uav.planner.params.dBar - uav.d0])
    #dirErr1.append(abs(planner1.params.dBar - planner1.params.d0))
    #spdErr1.append(abs(planner1.params.vBar - planner1.params.v0))
    #i=i+1

MAT = True

if(MAT):
    scio.savemat('mat/data/ACCsub/2turbs_1UAV_'+str(coverage)+'.mat',mdict={'UAV': uav.planner.error})
print('max coverage: '+str(100*np.amax(uav.mask_size)/grid_size)+'%')
print('min coverage: '+str(100*np.amin(uav.mask_size)/grid_size)+'%')
print('average coverage: ' + str(100*np.average(uav.mask_size)/grid_size) + '%')
filename = "UAV_2turb_15perc" # set the file name
# animate what happened

planner.plotHistory([uav], None, False) # don't save the animation, show the plot
# planner.plotHistory(uav, filename, False) # save the animation as an mp4, don't show the plot