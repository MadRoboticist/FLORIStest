## \file SingleUAV.py
# This is a script which utilizes functions in pathPlan_new.py
# while doing online updates of the sensitivity matrix and associated estimates

from VIZMAN_new import VisualizationManager
from pathPlan_new import PathPlanner
from UAV import UAV
import numpy as np
import scipy.io as scio
from FLORIS_iface import FLORIS_sub as Floris

# 'actual' values
vBar = 8.0 # meters per second
yBar = 0.0 # yaw in degrees
# estimated values
v0 = 13.0 # meters per second
y0 = 10.0 # yaw in degrees
# number of steps in memory / maximum map coverage
coverage = 0.15
# number of iterations to run the simulation
iterations = 5
# starting 'GPS' location of UAV
GPS = [2000,800]
# set to true to write the error data out to .mat files
MAT = False
# set to true to write the animation to mp4 file
MP4 = False
mp4filename = "UAV_2turb_15perc" # set the file name
# set to True to show the plot when finished
PLOT = True

############# JSON input for FLORIS model #####################
########### UNCOMMENT ONE #####################################
#JSONfile = "oneTurb_input_new.json" # one turbine case
JSONfile = "twoTurb_input_new.json" # two turbine case
#JSONfile = "36_turb_input_new.json" # 36 turbine case

########## RESOLUTION FOR FLORIS V. FLORIS ####################
########### UNCOMMENT ONE #####################################
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
# grid_resolution = [30, 15, 15] ## 1 turbine
grid_resolution = [48, 15, 15] ## 2 turbines
# grid_resolution = [60, 52, 15] ## 36 turbines

# set true to read xBar from SOWFAfile
SOWFA = False
SOWFAfile = "twoTurb_lowTI.u"

############### XBAR FROM SOWFA MODEL #####################
# read *.u file into list
if SOWFA:
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
else:
    xBar = None

## number of nodes on map
grid_size = grid_resolution[0]*grid_resolution[1]

## instatiate the Floris (interface)
WF = Floris(JSONfile, grid_resolution) # a JSON windfarm object read from a file
WF.set_vy(vBar, yBar)

print("Initializing UAV")
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager
vman.params.vBar = vBar
vman.params.yBar = np.deg2rad(yBar)
vman.params.v0 = v0
vman.params.y0 = np.deg2rad(y0)
## speed error threshold, set to negative to continue running after convergence
vman.params.spErrMax = -0.1
print("speed error threshold: "+str(vman.params.dirErrMax))
## direction error threshold, set negative to continue running after convergence
vman.params.yawErrMax = np.deg2rad(-0.1)
print("yaw error threshold: "+str(vman.params.yawErrMax))

## Initialize UAV
planner = PathPlanner(vman, xBar, True) # instantiate planner, Floris v. Floris, estimating yaw
uav = UAV(planner) # instantiate UAV object
uav.GPS = GPS # set start 'GPS'
# set UAV patrol boundaries
uav.minX = 0
uav.maxX = grid_resolution[0]
uav.minY = 0
uav.maxY = grid_resolution[1]
# set UAV estimates
uav.v0 = v0
uav.y0 = np.deg2rad(y0)
# calculate initial error
planner.error[0][1]=planner.params.yBar-uav.y0
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
uav.plan_horizon = int(uav.patrolMax*0.3)

planner.percent_plan = 1 #percentage of planned steps that will be taken


# run it for the indicated number of recalculations
#while dirErr[i]>planner1.params.dirErrMax or spdErr[i]>planner.params.spErrMax:
for i in range(iterations):
    print('\n\n\niteration ' + str(i))
    #uav.planner.COSPlan(uav)  # recalculate the plan
    uav.planner.greedyPath(uav)
    if not SOWFA:
        uav.move()  # move the UAV, static xBar/FLORIS v. FLORIS
    else:
        uav.moveTV() # move the UAV, time variant xBar/FLORIS v. SOWFA

    if (len(uav.GPSpath) >= uav.init_mask_size):
        uav.planner.updateEstimates_vy(uav)  # recalculate the estimates
        uav.planner.error.append([uav.planner.params.vBar - uav.v0,
                           uav.planner.params.yBar - uav.y0])

if(MAT):
    scio.savemat('mat/data/ACCsub/2turbs_1UAV_'+str(coverage)+'.mat',mdict={'UAV': uav.planner.error})
print('max coverage: '+str(100*np.amax(uav.mask_size)/grid_size)+'%')
print('min coverage: '+str(100*np.amin(uav.mask_size)/grid_size)+'%')
print('average coverage: ' + str(100*np.average(uav.mask_size)/grid_size) + '%')
# animate what happened
if not MP4:
    mp4filename = None

planner.plotHistory_vy([uav], mp4filename, PLOT)
