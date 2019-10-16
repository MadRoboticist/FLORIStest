## \file LearnPathOnline.py
# This is a script which utilizes the learnPath function in pathPlan.py
#

from old.visualization_manager_DJ import VisualizationManager
from pathPlan import PathPlanner
from UAV import UAV
import json
import numpy as np
import scipy.io as scio

with open("twoTurb_input.json") as WFJSON:
    ## a JSON windfarm object read from a file
    WF = json.load(WFJSON) # a JSON windfarm object read from a file
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
#grid_resolution = [60, 52, 15] # [x_res, y_res, z_res] 36_turb
grid_resolution = [48, 15, 15] # [x_res, y_res, z_res] twoTurb

grid_size = grid_resolution[0]*grid_resolution[1]
coverage = 0.15
## @var vman
# @brief VisualizationManager object
#
# see definition at visualization_manager_DJ.VisualizationManager
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager

## The initial wind speed estimate
#
# see definition at visualization_manager_DJ.VisualizationManager.params.v0
vman.params.v0 = 13.0 # the initial wind speed estimate
## The initial wind direction estimate
#
# see definition at visualization_manager_DJ.VisualizationManager.params.d0
vman.params.d0 = np.deg2rad(10.0) # the initial direction estimate
## Speed epsilon (ev) which is used to calculate the sensitivity matrix
#
# see definition at visualization_manager_DJ.VisualizationManager.params.epSpeed
vman.params.epSpeed = 0.01  # speed epsilon (ev)
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
vman.params.dirErrMax = np.deg2rad(0.1)  # direction error threshold
## Actual wind speed
#
# see definition at visualization_manager_DJ.VisualizationManager.params.vBar
vman.params.vBar = 8.0 # actual wind speed
## Actual wind direction in radians
#
# see definition at visualization_manager_DJ.VisualizationManager.params.dBar
vman.params.dBar = np.deg2rad(0.0) # actual wind direction

## A PathPlanner object
#
# see definition at pathPlan.PathPlanner
planner = PathPlanner(vman)

VTKpath = "L:\\SOWFAdata\\36Turb_lowTI\\postProcessing\\sliceDataInstantaneous\\21600"

VTKfilename = "U_slice_horizontal_1.vtk"
#readVTK(VTKpath+'\\'+VTKfilename,False)
#VTK = VTKreader(VTKpath, VTKfilename, vman)
#planner.Xbar = deepcopy(VTK.u_field)
## the number of learning iterations per planned path
planner.learn_length = 2000 #number of learning iterations
planner.learn_alpha = 0.9
planner.learn_gamma = 0.7
planner.learn_thresh = 0.5 #probability that the UAV will take a random step while learning
planner.max_score_box = 6# mask radius around a score maximum
planner.percent_plan = 0.90 #percentage of planned steps that will be taken
planner.learn_path_length = int(grid_size*coverage*2)
planner.MAX_BONUS = 1
## A UAV object
#
# see definition at UAV.UAV
UAV1 = UAV(planner)
## @var maskSUB
#   sets a subtractive penalty value to be applied to a node's
#   score mask each time it is visited
#   (set to 0 to apply no subtractive penalty)
UAV1.maskSUB = 100
## @var maskMUL
#   sets a multiplicative penalty value to be applied to a node's
#   score mask each time it is visited
#   (set to 1 to apply no multiplicative penalty)
UAV1.maskMUL = 1.5
## @var scoreWt
#   sets the weight of the sensitivity score in path calculations
UAV1.scoreWt = 0.0
## @var dSwt
#   sets the weight of the derivative of the sensitivity score
#   wrt transition in path calculations
UAV1.dSwt = 0
## @var d2Swt
#   sets the weight of the 2nd derivative of the sensitivity score
#   wrt transition in path calculations
UAV1.d2Swt = 0
## @var windWt
#   sets the weight of the wind direction vs. transition heading
#   in path calculations
UAV1.windWt = 0
## @var headWt
#   sets the weight of the UAV heading vs. transition heading
#   in path calculations
UAV1.headWt = 0
## @var waveWt
#   sets the weight of the wave map in path calculations
UAV1.waveWt = 0
## @var timeWt
#   sets the weight of time/iterations elapsed in path calculations
UAV1.timeWt = 0
UAV1.GPS = [2000,800]
## A parameter which decides how many moves the UAV will hold in memory.
# this is basically the number of nodes in the pseudoinverse mask
#
# see definition at pathPlan.UAV.patrolMax
UAV1.patrolMax = int(grid_size*coverage)
#print("patrol max: " + str(UAV1.patrolMax))
## A parameter which decides how far ahead the planner will work
# for the first iteration.
#
# see definition at pathPlan.UAV.plan_horizon
UAV1.init_plan_horizon = UAV1.patrolMax
## A parameter which decides how many moves the UAV will take on
# the first planned path before the first recalculation
#
# see definition at pathPlan.UAV.moves2recalc
UAV1.init_mask_size = int(UAV1.patrolMax*0.9)
## A parameter which decides how far ahead the planner will work
# after the first iteration.
#
# see definition at pathPlan.UAV.plan_horizon
UAV1.plan_horizon = int(grid_size*coverage*0.5)
## A parameter which decides how many moves the UAV will take on
# the planned path before it recalculates the estimates and the plan
#
# see definition at pathPlan.UAV.moves2recalc
UAV1.moves2recalc = int(UAV1.plan_horizon*0.5)
dirErr = list()
dirErr.append(abs(planner.params.dBar - planner.params.d0))
spdErr = list()
spdErr.append(abs(planner.params.vBar - planner.params.v0))
# run it for the indicated number of recalculations
i=0
while dirErr[i]>planner.params.dirErrMax or spdErr[i]>planner.params.spErrMax:
#for i in range(10):

    planner.learnPath(UAV1) # recalculate the plan
    UAV1.move() # move the UAV
    if(len(UAV1.GPSpath)>=UAV1.init_mask_size):
        print('iteration ' + str(i))
        i = i + 1
        planner.updateEstimates(UAV1) # recalculate the estimates
    dirErr.append(abs(planner.params.dBar - planner.params.d0))
    spdErr.append(abs(planner.params.vBar - planner.params.v0))


MAT = False

if(MAT):
    scio.savemat('mat/2turbs_'+str(vman.params.vBar) + '_' + \
         str(np.rad2deg(vman.params.dBar)) + '_' + str(int(100*coverage)) + \
         '_spdERR_UAV.mat',mdict={'spdErrUAV'+str(int(100*coverage)): spdErr})
    scio.savemat('mat/2turbs_' + str(vman.params.vBar) + '_' + \
         str(np.rad2deg(vman.params.dBar)) + '_' + str(int(100*coverage)) + \
         '_dirERR_UAV.mat', mdict={'dirErrUAV'+str(int(100*coverage)): dirErr})

# planner.printAVGs()
filename = "UAV_2turb_Qlearn3" # set the file name
# animate what happened
#planner.plotHistory(UAV1, None, True) # don't save the animation, show the plot
planner.plotHistory(UAV1, filename, False) # save the animation as an mp4, don't show the plot