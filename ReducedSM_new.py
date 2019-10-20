## \file ReducedSM.py
# This is a script which runs the reducedSM function from visualization_manager_DJ.py
#

from VIZMAN_new import VisualizationManager
from FLORIS_iface import FLORIS_sub as Floris
import numpy as np

json_file = "twoTurb_input_new.json"
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
grid_resolution = [48, 15, 15] # [x_res, y_res, z_res]
FL = Floris(json_file, grid_resolution)
FL.set_incoming(8.0,270.0)
## VisualizationManager object
vman = VisualizationManager(FL, grid_resolution) # set up the visualization manager

# the following are the parameters used by the reducedSM() function
## The initial speed estimate
vman.params.v0 = 8.0 # initial speed estimate
## The initial direction estimate
vman.params.d0 = np.deg2rad(270.0) # initial direction estimate
## speed epsilon (ev)
vman.params.epSpeed = 0.01  # speed epsilon (ev)
## direction epsilon (ed)
vman.params.epDir = 0.001  # direction epsilon (ed)
## speed error minimum threshold
vman.params.spErrMax = -0.1  # speed error threshold
print("speed error threshold: "+str(vman.params.dirErrMax))
## direction error minimum threshold
vman.params.dirErrMax = -np.deg2rad(0.1)  # direction error threshold
print("direction error threshold: "+str(vman.params.dirErrMax))
## maximum number of iterations to be completed
vman.params.iterMax = 10 # iteration threshold
## The actual wind speed
vman.params.vBar = 13.0 # actual wind speed
## The actual wind direction
vman.params.dBar = np.deg2rad(280.0) # actual wind direction
## A masking threshold for reducing the sensitivity matrix
vman.params.mask_thresh = 0.5 # masking threshold for reducing normalized sensitivity matrix
## A dampening factor for the algorithm
vman.params.damping = 1.0

vman.reducedSM(False, False, None, None) # output to .mat file, no animation, no mp4, xBar from FLORIS