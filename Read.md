/*
@mainpage

The FLORIS model is being used to make online predictions about the best locations in a windfarm 
to take measurements. From this, UAV trajectories will be generated. 
The files included in this documentation have been created 
(or altered from the Floris source code) in order to carry out this research.

visualization_manager_DJ.py

This package contains functions which handle visualizations of Floris data.

AnimateError.py

This is a test script which runs the animateDnSerror function from visualization_manager_DJ.py

PlotSensitivityMatrix.py

This is a test script which runs the plotSensitivityMatrix function from visualization_manager_DJ.py

ReducedSM.py

This is a script which runs the reducedSM function from visualization_manager_DJ.py

example_input.json

This is an example JSON object which is used as the input parameters to create a Floris windfarm simulation.

input_reader.py

This file was altered from the original to allow a Floris object to be created from a JSON object rather than just a file.