## \file MultiUAV.py
# This is a script which utilizes the COSplan function in pathPlan.py
# while doing online updates of the sensitivity matrix and associated estimates
# using a swarm of UAVs with independent consensus-based nearest-neighbor estimates
# and truth data from a static FLORIS model

from VIZMAN_new import VisualizationManager
from pathPlan_new import PathPlanner
from UAV import UAV
from FLORIS_iface import FLORIS_sub as Floris
import numpy as np
import scipy.io as scio

# 'actual' values
vBar = 8.0 # meters per second
dBar = 280.0 # degrees counterclockwise from North

############### XBAR FROM SOWFA MODEL #####################
# read *.u file into list
'''
print("reading in xBar from SOWFA")
SOWFAfile = "36_Turb_lowTI.u"
xBar = []
with open(SOWFAfile) as f:
    lines = []  # list to collect lines
    while 1:
        aline = f.readline()
        if aline.strip():
            lines.append(aline)     # nonempty line
        else:              # empty line
            if len(lines)==0: break
            xBar.append(np.loadtxt(lines, dtype=int))
            lines = []

xBar = np.array(xBar) # convert to array
grid_resolution = [xBar.shape[1], xBar.shape[2], 15] # [x_res, y_res, z_res]
'''
########## RESOLUTION FOR FLORIS V. FLORIS ####################
########### UNCOMMENT THE NEXT TWO LINES TO USE FLORIS XBAR ######################
grid_resolution = [60, 52, 15] ## 36 turbines
xBar = None

grid_size = grid_resolution[0]*grid_resolution[1]
coverage = 0.15

### set up FLORIS model and vman
print("starting FLORIS")
JSONfile = "36_turb_input_new.json"
WF = Floris(JSONfile, grid_resolution)
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager
vman.params.vBar = vBar
vman.params.dBar = np.deg2rad(dBar)

## Initialize UAVs
vman.params.d0 = np.deg2rad(270.0+45.0)
vman.params.v0 = 5.5
print("Initializing UAV1")
planner1 = PathPlanner(vman, xBar)
UAV1 = UAV(planner1)
UAV1.GPS = [350,1700]
UAV1.minX = 0
UAV1.maxX = 20
UAV1.minY = 0
UAV1.maxY = 18
UAV1.d0 = planner1.params.d0
planner1.error[0][1]=planner1.params.dBar-UAV1.d0
UAV1.v0 = planner1.params.v0
planner1.error[0][0]=planner1.params.vBar-UAV1.v0

print("Initializing UAV2")
vman.params.d0 = np.deg2rad(270.0+33.75)
vman.params.v0 = 6.875
planner2 = PathPlanner(vman, xBar)
UAV2 = UAV(planner2)
UAV2.GPS = [2200,350]
UAV2.minX = 20
UAV2.maxX = 40
UAV2.minY = 0
UAV2.maxY = 18
UAV2.d0 = planner2.params.d0
planner2.error[0][1]=planner2.params.dBar-UAV2.d0
UAV2.v0 = planner2.params.v0
planner2.error[0][0]=planner2.params.vBar-UAV2.v0

print("Initializing UAV3")
vman.params.d0 = np.deg2rad(270.0+22.5)
vman.params.v0 = 8.75
planner3 = PathPlanner(vman, xBar)
UAV3 = UAV(planner3)
UAV3.GPS = [4000,1700]
UAV3.minX = 40
UAV3.maxX = 60
UAV3.minY = 0
UAV3.maxY = 18
UAV3.d0 = planner3.params.d0
planner3.error[0][1]=planner3.params.dBar-UAV3.d0
UAV3.v0 = planner3.params.v0
planner3.error[0][0]=planner3.params.vBar-UAV3.v0

print("Initializing UAV4")
vman.params.d0 = np.deg2rad(270.0+11.25)
vman.params.v0 = 10.625
planner4 = PathPlanner(vman, xBar)
UAV4 = UAV(planner4)
UAV4.GPS = [350,3200]
UAV4.minX = 0
UAV4.maxX = 20
UAV4.minY = 18
UAV4.maxY = 36
UAV4.d0 = planner4.params.d0
planner4.error[0][1]=planner4.params.dBar-UAV4.d0
UAV4.v0 = planner4.params.v0
planner4.error[0][0]=planner4.params.vBar-UAV4.v0

print("Initializing UAV5")
vman.params.d0 = np.deg2rad(270.0)
vman.params.v0 = 12.5
planner5 = PathPlanner(vman, xBar)
UAV5 = UAV(planner5)
UAV5.GPS = [3500,1800]
UAV5.minX = 20
UAV5.maxX = 40
UAV5.minY = 18
UAV5.maxY = 36
UAV5.d0 = planner5.params.d0
planner5.error[0][1]=planner5.params.dBar-UAV5.d0
UAV5.v0 = planner5.params.v0
planner5.error[0][0]=planner5.params.vBar-UAV5.v0

print("Initializing UAV6")
vman.params.d0 = np.deg2rad(270.0-11.25)
vman.params.v0 = 14.375
planner6 = PathPlanner(vman, xBar)
UAV6 = UAV(planner6)
UAV6.GPS = [4000,3200]
UAV6.minX = 40
UAV6.maxX = 60
UAV6.minY = 18
UAV6.maxY = 36
UAV6.d0 = planner6.params.d0
planner6.error[0][1]=planner6.params.dBar-UAV6.d0
UAV6.v0 = planner6.params.v0
planner6.error[0][0]=planner6.params.vBar-UAV6.v0

print("Initializing UAV7")
vman.params.d0 = np.deg2rad(270.0-22.5)
vman.params.v0 = 16.25
planner7 = PathPlanner(vman, xBar)
UAV7 = UAV(planner7)
UAV7.GPS = [350,4500]
UAV7.minX = 0
UAV7.maxX = 20
UAV7.minY = 36
UAV7.maxY = 52
UAV7.d0 = planner7.params.d0
planner7.error[0][1]=planner7.params.dBar-UAV7.d0
UAV7.v0 = planner7.params.v0
planner7.error[0][0]=planner7.params.vBar-UAV7.v0

print("Initializing UAV8")
vman.params.d0 = np.deg2rad(270.0-33.75)
vman.params.v0 = 18.125
planner8 = PathPlanner(vman, xBar)
UAV8 = UAV(planner8)
UAV8.GPS = [2500,3600]
UAV8.minX = 20
UAV8.maxX = 40
UAV8.minY = 36
UAV8.maxY = 52
UAV8.d0 = planner8.params.d0
planner8.error[0][1]=planner8.params.dBar-UAV8.d0
UAV8.v0 = planner8.params.v0
planner8.error[0][0]=planner8.params.vBar-UAV8.v0

print("Initializing UAV9")
vman.params.d0 = np.deg2rad(270.0-45.0)
vman.params.v0 = 20.0
planner9 = PathPlanner(vman, xBar)
UAV9 = UAV(planner9)
UAV9.GPS = [4000,4500]
UAV9.minX = 40
UAV9.maxX = 60
UAV9.minY = 36
UAV9.maxY = 52
UAV9.d0 = planner9.params.d0
planner9.error[0][1]=planner9.params.dBar-UAV9.d0
UAV9.v0 = planner9.params.v0
planner9.error[0][0]=planner9.params.vBar-UAV9.v0

UAVs=[UAV1, UAV2, UAV3, UAV4, UAV5, UAV6, UAV7, UAV8, UAV9]

for uav in UAVs:
## A parameter which decides how many moves the UAV will hold in memory.
# this is basically the number of nodes in the pseudoinverse mask
#
# see definition at pathPlan.UAV.patrolMax
    uav.patrolMax = int(grid_size*coverage/len(UAVs))

## A parameter which decides how far ahead the planner will work
# for the first iteration.
#
# see definition at pathPlan.UAV.plan_horizon
    uav.init_plan_horizon = UAV1.patrolMax

## A parameter which decides how many moves the UAV will take on
# the first planned path before the first recalculation
#
# see definition at pathPlan.UAV.moves2recalc
    uav.init_mask_size = int(UAV1.patrolMax*0.9)

## A parameter which decides how far ahead the planner will work
# after the first iteration.
#
# see definition at pathPlan.UAV.plan_horizon
    uav.plan_horizon = int(grid_size*coverage*0.3/len(UAVs))

## A parameter which decides how many moves the UAV will take on
# the planned path before it recalculates the estimates and the plan
#
# see definition at pathPlan.UAV.moves2recalc
    uav.moves2recalc = int(UAV1.plan_horizon*0.9)

################## COSPlan params ########################
## @var maskSUB
#   sets a subtractive penalty value to be applied to a node's
#   score mask each time it is visited
#   (set to 0 to apply no subtractive penalty)
    uav.maskSUB = 1.0
## @var maskMUL
#   sets a multiplicative penalty value to be applied to a node's
#   score mask each time it is visited
#   (set to 1 to apply no multiplicative penalty)
    uav.maskMUL = 0.5
## @var maskMULthenSUB
#   boolean which decides the order maskSUB and maskMUL are carried out
    uav.maskMULthenSUB = True
## @var scoreWt
#   sets the weight of the sensitivity score in path calculations
    uav.scoreWt = 0
## @var dSwt
#   sets the weight of the derivative of the sensitivity score
#   wrt transition in path calculations
    uav.dSwt = 0
## @var d2Swt
#   sets the weight of the 2nd derivative of the sensitivity score
#   wrt transition in path calculations
    uav.d2Swt = 0
## @var windWt
#   sets the weight of the wind direction vs. transition heading
#   in path calculations
    uav.windWt = 0
## @var headWt
#   sets the weight of the UAV heading vs. transition heading
#   in path calculations
    uav.headWt = 0
## @var waveWt
#   sets the weight of the wave map in path calculations
    uav.waveWt = 0
## @var timeWt
#   sets the weight of time/iterations elapsed in path calculations
    uav.timeWt = 0

planner1.percent_plan = 1 #percentage of planned steps that will be taken
# run it for the indicated number of recalculations
i=0

NNradius = 2200
alpha = 0.75

for i in range(30):
    print('\n\n\niteration ' + str(i))
    for ind, uav in enumerate(UAVs):
        print("UAV"+str(ind+1))
        uav.planner.COSPlan(uav)  # recalculate the plan
        #uav.move()  # move the UAV, static xBar/FLORIS v. FLORIS
        uav.moveTV() # move the UAV, time varying xBar/FLORIS v. SOWFA
        if (len(uav.GPSpath) >= uav.init_mask_size):
            uav.planner.updateEstimates(uav)  # recalculate the estimates

            Knn = 0
            Vnn = 0
            Dnn = 0
            for indx, nUAV in enumerate(UAVs):
                if ind != indx:
                    if uav.planner.euclidDist(uav.GPS,nUAV.GPS)<NNradius:
                        Knn = Knn + 1
                        Vnn = Vnn + nUAV.v0
                        Dnn = Dnn + nUAV.d0
            if Knn>0:
                print(str(Knn)+" nearest neighbors.")
                uav.v0 = (uav.v0 * alpha) + ((Vnn * (1 - alpha)) / Knn)
                uav.d0 = (uav.d0 * alpha) + ((Dnn * (1 - alpha)) / Knn)
                print("Adjusted Estimate: " + str([uav.v0, uav.d0]))
                print("Adjusted Error: ["+str(uav.planner.params.vBar-uav.v0)+", "+\
                      str(uav.planner.params.dBar-uav.d0)+"]")
            uav.planner.error.append([uav.planner.params.vBar - uav.v0,
                               uav.planner.params.dBar - uav.d0])

# set true to output error to .mat files
MAT = True

if(MAT):
    for ind, uav in enumerate(UAVs):
        scio.savemat('mat/36turbs_UAV' + str(ind) + '.mat',mdict={'UAV'+str(ind): uav.planner.error})

#planner.printAVGs()
filename = "UAV_36turb_swarm_75_varinit_15perc_NEW_FL" # set the file name
# animate what happened
planner = PathPlanner(vman)

planner.plotHistory(UAVs, None, False) # don't save the animation, show the plot
#planner.plotHistory(UAVs, fiilename, False) # save the animation as an mp4, don't show the plot

