# Copyright 2019 NREL

# Licensed under the Apache License, Version 2.0 (the "License"); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.

# See read the https://floris.readthedocs.io for documentation
import numpy as np
import matplotlib.pyplot as plt
import floris.tools as wfct
import floris.simulation as sim
import json

class FLORIS_sub:

    def __init__(self, json_file, grid_resolution):
        # Initialize the FLORIS interface fi
        self.grid_resolution = grid_resolution
        self.fi = wfct.floris_utilities.FlorisInterface(json_file)
        #print("getting turbine yaw angles")
        self.yaw = self.fi.get_yaw_angles()
        #print("setting up horizontal plane")
        self.hor_plane = wfct.cut_plane.HorPlane(
            self.fi.get_hub_height_flow_data(),
            self.fi.floris.farm.turbines[0].hub_height
        )
        #print("horizontal plane defined")
        self.set_resolution(grid_resolution)
        #print("calculating u-field")
        self.update_ufield()
        #print("u-field calculated")
        with open(json_file) as WFJSON:
            ## a JSON windfarm object read from a file
            WF = json.load(WFJSON)  # a JSON windfarm object read from a file
        self.wind_speed = WF['farm']['properties']['wind_speed']
        self.wind_direction = WF['farm']['properties']['wind_direction']
        self.turbines = self.fi.floris.farm.flow_field.turbine_map.sorted_in_x_as_list()


    def set_resolution(self, grid_resolution):
        self.grid_resolution = grid_resolution
        wfct.cut_plane.change_resolution(self.hor_plane, self.grid_resolution)
        self.update_ufield()

    def calculate_wake(self):
    # Calculate wake
        self.fi.calculate_wake()
        self.update_ufield()

    def set_incoming(self, v, th):
        self.fi.reinitialize_flow_field(v, th)
        self.wind_speed = v
        self.wind_direction = th
        self.update_ufield()

    def set_yaw(self, yaw):
        self.fi.calculate_wake(yaw)
        self.update_ufield()
        # Initialize the horizontal cut

    def update_ufield(self):
        self.hor_plane = wfct.cut_plane.HorPlane(
            self.fi.get_hub_height_flow_data(),
            self.fi.floris.farm.turbines[0].hub_height
        )
        wfct.cut_plane.change_resolution(self.hor_plane, self.grid_resolution)
        self.u_field = np.asarray(self.hor_plane.u_mesh.reshape(self.grid_resolution[1],self.grid_resolution[0])).transpose()
        self.x_mesh = np.asarray(self.hor_plane.x1_mesh.reshape(self.grid_resolution[1],self.grid_resolution[0])).transpose()
        self.y_mesh = np.asarray(self.hor_plane.x2_mesh.reshape(self.grid_resolution[1],self.grid_resolution[0])).transpose()
