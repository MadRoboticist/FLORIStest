## \file TV_SM_SOWFA.py
# This is a script which runs the reducedSM function from visualization_manager_DJ.py
# using time-varying SOWFA data from a *.u file

from VIZMAN_new import VisualizationManager
from FLORIS_iface import FLORIS_sub as Floris
import numpy as np

# read *.u file into list
SOWFAfile = "twoTurb_lowTI.u"
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
#print(xBar.shape)
grid_resolution = [xBar.shape[1], xBar.shape[2], 15] # [x_res, y_res, z_res]
#print("grid_resolution: ",grid_resolution)
JSONfile = "twoTurb_input_new.json"
WF = Floris(JSONfile, grid_resolution)## grab the X,Y dimensions from xBar

## VisualizationManager object
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager

# the following are the parameters used by the reducedSM() function
## The initial speed estimate
vman.params.v0 = 8.0 # initial speed estimate
## The initial direction estimate
vman.params.d0 = np.deg2rad(270.0) # initial direction estimate
## speed epsilon (ev)
vman.params.epSpeed = 0.001  # speed epsilon (ev)
## direction epsilon (ed)
vman.params.epDir = 0.0001  # direction epsilon (ed)
## speed error minimum threshold
vman.params.spErrMax = -0.01  # speed error threshold
## direction error minimum threshold
vman.params.dirErrMax = -0.001  # direction error threshold
## The actual wind speed
vman.params.vBar = 8.0 # actual wind speed
## The actual wind direction
vman.params.dBar = np.deg2rad(270.0) # actual wind direction
## A masking threshold for reducing the sensitivity matrix
vman.params.mask_thresh = 0.0 # masking threshold for reducing normalized sensitivity matrix
# 0 for no mask
vman.params.iterMax = len(xBar) # iteration threshold, set to number of SOWFA frames
#vman.params.iterMax = 10
mp4file = "twoTurbSOWFA_reducedSM_0init-err_no-mask_new"
# run the function, save .mat file, animate plot, save mp4, use xBar from SOWFA, do not show plot.
vman.reducedSM(True, False, None, xBar, False)
