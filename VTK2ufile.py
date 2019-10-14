import VTKreader
from copy import deepcopy
from visualization_manager_DJ import VisualizationManager
import json
import numpy as np

###### CHOOSE ONE ########
JSONfile = "oneTurb_input.json"
#JSONfile = "twoTurb_input.json"
#JSONfile = "36_turb_input.json"

with open(JSONfile) as WFJSON:
    ## a JSON windfarm object read from a file
    WF = json.load(WFJSON) # a JSON windfarm object read from a file

######### CHOOSE ONE #############
## flowfield resolution given as [x_resolution, y_resolution, z_resolution]
grid_resolution = [30, 15, 15] # one turbine
#grid_resolution = [48, 15, 15] # two turbines
#grid_resolution = [60, 52, 15] # 36 turbines

## VisualizationManager object
vman = VisualizationManager(WF, grid_resolution) # set up the visualization manager

############ CHOOSE ONE ###############
VTKpath = "L:\\SOWFAdata\\oneTurb_lowTI\\postProcessing\\sliceDataInstantaneous" # one turbine
#VTKpath = "L:\\SOWFAdata\\twoTurb_lowTI\\postProcessing\\sliceDataInstantaneous" # two turbines
#VTKpath = "L:\\SOWFAdata\\WF_lowTI\\postProcessing\\sliceDataInstantaneous" # 36 turbines

############## CHOOSE ONE #############
#VTKrange = [20145, 21875] #two Turbines, cut initialization
#VTKrange = [20005, 21875] #two Turbines, leave in initialization
#VTKrange = [20005, 22000] #one or 36 Turbines, leave in initialization
VTKrange = [20145, 22000] #one or 36 Turbines, cut initialization

increment = 5
VTKfilename = "U_slice_horizontal_1.vtk"
outputFile = "oneTurb_lowTI.u"
# run the reducedSM function

VTKreader.readVTKtoFile(VTKpath, VTKfilename, VTKrange, increment, outputFile, vman)
