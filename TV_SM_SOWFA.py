## \file TV_SM_SOWFA.py
# This is a script which runs the reducedSM function from visualization_manager_DJ.py
# using time-varying SOWFA data from a *.u file

from visualization_manager_DJ import VisualizationManager
import json
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
print(xBar.shape) # print the shape

JSONfile = "twoTurb_input.json"
with open(JSONfile) as WFJSON:
    ## a JSON windfarm object read from a file
    WF = json.load(WFJSON) # a JSON windfarm object read from a file
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
## grab the X,Y dimensions from xBar
grid_resolution = [xBar.shape[1], xBar.shape[2], 15] # [x_res, y_res, z_res]

## VisualizationManager object
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager

# the following are the parameters used by the reducedSM() function
## The initial speed estimate
vman.params.v0 = 13.0 # initial speed estimate
## The initial direction estimate
vman.params.d0 = np.deg2rad(10.0) # initial direction estimate
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
vman.params.dBar = np.deg2rad(0.0) # actual wind direction
## A masking threshold for reducing the sensitivity matrix
vman.params.mask_thresh = 0.0 # masking threshold for reducing normalized sensitivity matrix
# 0 for no mask
vman.params.iterMax = len(xBar) # iteration threshold, set to number of SOWFA frames

mp4file = "twoTurbSOWFA_reducedSM_non0init-err_no-mask"
# run the function, save .mat file, animate plot, save mp4, use xBar from SOWFA, do not show plot.
vman.reducedSM(True, True, mp4file, xBar, False)
