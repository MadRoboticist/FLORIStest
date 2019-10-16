## \file ReducedSM.py
# This is a script which runs the reducedSM function from visualization_manager_DJ.py
#

from old.visualization_manager_DJ import VisualizationManager
import json
import numpy as np

with open("twoTurb_input.json") as WFJSON:
    ## a JSON windfarm object read from a file
    WF = json.load(WFJSON) # a JSON windfarm object read from a file
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
grid_resolution = [48, 15, 15] # [x_res, y_res, z_res]

## VisualizationManager object
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager

# the following are the parameters used by the reducedSM() function
## The initial speed estimate
vman.params.v0 = 8.0 # initial speed estimate
## The initial direction estimate
vman.params.d0 = np.deg2rad(0.0) # initial direction estimate
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
vman.params.iterMax = 6 # iteration threshold
## The actual wind speed
vman.params.vBar = 13.0 # actual wind speed
## The actual wind direction
vman.params.dBar = np.deg2rad(10.0) # actual wind direction
## A masking threshold for reducing the sensitivity matrix
vman.params.mask_thresh = 0.0 # masking threshold for reducing normalized sensitivity matrix
## A dampening factor for the algorithm
vman.params.damping = 1.0
VTKpath = "L:\\SOWFAdata\\oneTurb_lowTI\\postProcessing\\slicedatainstantaneous\\22000"

VTKfilename = "T_slice_horizontal_1.vtk"
# VTK = VTKreader(VTKpath, VTKfilename, vman)
# xBar = deepcopy(VTK.u_field)
# run the reducedSM function
vman.reducedSM(True, False, None, None) # output to .mat file, no animation, no mp4, xBar from FLORIS