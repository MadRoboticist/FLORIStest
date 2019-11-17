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
from copy import deepcopy

# 'actual' values
vBar = 8.0  # meters per second
yBar = 0.0  # yaw in degrees
iterations = 10
# initial estimates
v0 = [vBar-2.0, vBar-1.0, vBar, vBar+1.0, vBar+2.0, vBar+3.0, vBar+4.0, vBar+5.0, vBar+6.0]
y0 = [yBar-40.0, yBar-30.0, yBar-20.0, yBar-10.0, yBar, yBar+10.0, yBar+20.0, yBar+30.0, yBar+40.0]
# UAV index boundaries
xMin = [0, 20, 40, 0, 20, 40, 0, 20, 40]
xMax = [20, 40, 60, 20, 40, 60, 20, 40, 60]
yMin = [0, 0, 0, 18, 18, 18, 36, 36, 36]
yMax = [18, 18, 18, 36, 36, 36, 52, 52, 52]
# UAV starting 'GPS' positions
GPS = [[350, 1700],
       [2200, 350],
       [4000, 1700],
       [350, 3200],
       [3500, 1800],
       [4000, 3200],
       [350, 4500],
       [2500, 3600],
       [4000, 4500]]
# Nearest Neighbor radius
NNradius = 2200
# weight of UAV's own estimate in NN average
alpha = 0.75
# set total amount of coverage
coverage = 0.15
# set true to output errors to .mat files
MAT = True
# set True to output animation to .mp4 file (will not show plot when finished)
MP4 = False
mp4filename = "UAV_36turb_swarm_yaw"  # set the file name
# set True to show plot when finished
PLOT = True

# ######### RESOLUTION FOR FLORIS V. FLORIS ####################
JSONfile = "36_turb_input_new.json"
grid_resolution = [60, 52, 15]  # 36 turbines
xBar = None  # xBar is None for FLORIS V. FLORIS
SOWFA = False  # set to True to use SOWFA input from SOWFAfile
SOWFAfile = "36_Turb_lowTI.u"

# ############## XBAR FROM SOWFA MODEL #####################
# read *.u file into list
if SOWFA:
    print("reading in xBar from SOWFA")
    xBar = []
    with open(SOWFAfile) as f:
        lines = []  # list to collect lines
        while 1:
            aline = f.readline()
            if aline.strip():
                lines.append(aline)     # nonempty line
            else:              # empty line
                if len(lines) == 0:
                    break
                xBar.append(np.loadtxt(lines, dtype=int))
                lines = []
    xBar = np.array(xBar)  # convert to array
    grid_resolution = [xBar.shape[1], xBar.shape[2], 15]  # [x_res, y_res, z_res]

# Number of nodes on the grid:
grid_size = grid_resolution[0]*grid_resolution[1]

# set up FLORIS model and vman
print("starting FLORIS")
WF = Floris(JSONfile, grid_resolution)
WF.set_vy(vBar, yBar)
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager
vman.params.vBar = vBar
vman.params.yBar = np.deg2rad(yBar)

# set up UAVs
UAVs = list()
for i in range(len(GPS)):
    print("Initializing UAV", i+1)
    vman.params.y0 = np.deg2rad(y0[i])
    vman.params.v0 = v0[i]
    UAVs.append(UAV(PathPlanner(vman, xBar, True)))
    UAVs[i].GPS = GPS[i]
    UAVs[i].minX = xMin[i]
    UAVs[i].maxX = xMax[i]
    UAVs[i].minY = yMin[i]
    UAVs[i].maxY = yMax[i]
    UAVs[i].y0 = deepcopy(UAVs[i].planner.params.y0)
    UAVs[i].planner.error[0][1] = UAVs[i].planner.params.yBar-UAVs[i].y0
    UAVs[i].v0 = deepcopy(UAVs[i].planner.params.v0)
    UAVs[i].planner.error[0][0] = UAVs[i].planner.params.vBar-UAVs[i].v0

    ## A parameter which decides how many moves the UAV will hold in memory.
    # this is basically the number of nodes in the pseudoinverse mask
    #
    # see definition at pathPlan.UAV.patrolMax
    UAVs[i].patrolMax = int(grid_size*coverage/len(GPS))

    ## A parameter which decides how far ahead the planner will work
    # for the first iteration.
    #
    # see definition at pathPlan.UAV.plan_horizon
    UAVs[i].init_plan_horizon = UAVs[i].patrolMax

    ## A parameter which decides how many moves the UAV will take on
    # the first planned path before the first recalculation
    #
    # see definition at pathPlan.UAV.moves2recalc
    UAVs[i].init_mask_size = int(UAVs[i].patrolMax*0.9)

    ## A parameter which decides how far ahead the planner will work
    # after the first iteration.
    #
    # see definition at pathPlan.UAV.plan_horizon
    UAVs[i].plan_horizon = int(grid_size*coverage*0.3/len(GPS))

    ## A parameter which decides how many moves the UAV will take on
    # the planned path before it recalculates the estimates and the plan
    #
    # see definition at pathPlan.UAV.moves2recalc
    UAVs[i].moves2recalc = int(UAVs[i].plan_horizon*0.9)

    # ################# COSPlan params ########################
    ## @var maskSUB
    #   sets a subtractive penalty value to be applied to a node's
    #   score mask each time it is visited
    #   (set to 0 to apply no subtractive penalty)
    UAVs[i].maskSUB = 1.0
    ## @var maskMUL
    #   sets a multiplicative penalty value to be applied to a node's
    #   score mask each time it is visited
    #   (set to 1 to apply no multiplicative penalty)
    UAVs[i].maskMUL = 0.5
    ## @var scoreWt
    #   sets the weight of the sensitivity score in path calculations
    UAVs[i].scoreWt = 0
    ## @var dSwt
    #   sets the weight of the derivative of the sensitivity score
    #   wrt transition in path calculations
    UAVs[i].dSwt = 0
    ## @var d2Swt
    #   sets the weight of the 2nd derivative of the sensitivity score
    #   wrt transition in path calculations
    UAVs[i].d2Swt = 0
    ## @var windWt
    #   sets the weight of the wind direction vs. transition heading
    #   in path calculations
    UAVs[i].windWt = 0
    ## @var headWt
    #   sets the weight of the UAV heading vs. transition heading
    #   in path calculations
    UAVs[i].headWt = 0
    ## @var waveWt
    #   sets the weight of the wave map in path calculations
    UAVs[i].waveWt = 0
    ## @var timeWt
    #   sets the weight of time/iterations elapsed in path calculations
    UAVs[i].timeWt = 0
    # percentage of planned steps that will be taken (x100)
    UAVs[i].planner.percent_plan = 1

# run it for the indicated number of recalculations
for i in range(iterations):
    print('\n\n\niteration ' + str(i))
    for ind, uav in enumerate(UAVs):
        print("UAV"+str(ind+1))
        uav.planner.COSPlan(uav)  # recalculate the plan
        if SOWFA:
            uav.moveTV()  # move the UAV, time varying xBar/FLORIS v. SOWFA
        else:
            uav.move()  # move the UAV, static xBar/FLORIS v. FLORIS
        if len(uav.GPSpath) >= uav.init_mask_size:
            uav.planner.updateEstimates(uav)  # recalculate the estimates
            Knn = 0
            Vnn = 0
            Ynn = 0
            for indx, nUAV in enumerate(UAVs):
                if ind != indx:
                    if uav.planner.euclidDist(uav.GPS, nUAV.GPS) < NNradius:
                        Knn = Knn + 1
                        Vnn = Vnn + nUAV.v0
                        Ynn = Ynn + nUAV.y0
            if Knn > 0:
                print(str(Knn)+" nearest neighbors.")
                uav.v0 = (uav.v0 * alpha) + ((Vnn * (1 - alpha)) / Knn)
                uav.y0 = (uav.y0 * alpha) + ((Ynn * (1 - alpha)) / Knn)
                print("Adjusted Estimate: " + str([uav.v0, uav.y0]))
                print("Adjusted Error: ["+str(uav.planner.params.vBar-uav.v0)+", " +
                      str(uav.planner.params.yBar-uav.y0)+"]")
            uav.planner.error.append([uav.planner.params.vBar - uav.v0,
                                      uav.planner.params.yBar - uav.y0])

if MAT:
    for ind, uav in enumerate(UAVs):
        scio.savemat('mat/36turbs_UAV' + str(ind) + '.mat', mdict={'UAV'+str(ind): uav.planner.error})
if not MP4:
    mp4filename = None

# animate what happened
planner = PathPlanner(vman, YAW=True)
planner.plotHistory(UAVs, mp4filename, PLOT)  # don't save the animation, show the plot
