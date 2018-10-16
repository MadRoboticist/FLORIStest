## \file visualization_manager_DJ.py
#
#   This file contains functions which handle visualization of Floris data.
#
# Copyright 2017 NREL
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use
# this file except in compliance with the License. You may obtain a copy of the
# License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed
# under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
# CONDITIONS OF ANY KIND, either express or implied. See the License for the
# specific language governing permissions and limitations under the License.
#

from floris.coordinate import Coordinate
from floris.floris import Floris
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Slider, Button
import numpy as np
from scipy.interpolate import griddata
from copy import deepcopy


## \class VisualizationManager
# Contains functions which handle visualizations of Floris data.
#
# This version of the FLORIS VisualizationManager class was altered
#    to instantiate its own FLORIS object from a JSON windfarm object,
#    and also to include some new functions:
#
#    1. animateDnSerror
#        Plot error against a base FLORIS model
#        compared to a range of different directions and speeds.
#
#    2. plotSensitivityMatrix
#        Plotting the Sensitivity matrix of a Z plane
#        at a given height in the floris model.
#
#    3. reducedSM
#        Plotting the convergence of speed and direction estimates
#        with a reduced sensitivity matrix
#
class VisualizationManager():

    # General plotting functions for default FLORIS

    ## Class constructor
    #   @param WF takes a JSON windfarm object
    #   @param grid_resolution takes a 3-element integer array
    #           with [x_resolution, y_resolution_ z_resolution]
    #   @return a VisualizationManager object
    def __init__(self, WF, grid_resolution=(100, 100, 25)):
        self.WF = WF
        self.figure_count = 0
        self.flowfield = Floris(self.WF).farm.flow_field
        #print(self.flowfield.u_field)
        self.grid_res = grid_resolution
        self.grid_resolution = Coordinate(grid_resolution[0], grid_resolution[1], grid_resolution[2])
        self._initialize_flowfield_for_plotting()
        self._boolBTN = False

    def _set_texts(self, plot_title, horizontal_axis_title, vertical_axis_title):
        fontsize = 15
        plt.title(plot_title, fontsize=fontsize)
        plt.xlabel(horizontal_axis_title, fontsize=fontsize)
        plt.ylabel(vertical_axis_title, fontsize=fontsize)

    def _set_colorbar(self):
        cb = plt.colorbar()
        cb.ax.tick_params(labelsize=15)

    def _set_axis(self):
        plt.axis('equal')
        plt.tick_params(which='both', labelsize=15)

    def _new_figure(self):
        plt.figure(self.figure_count)
        self.figure_count += 1

    def _new_filled_contour(self, mesh1, mesh2, data):
        self._new_figure()
        #vmax = np.amax(data)
        plt.contourf(mesh1, mesh2, data, 50,
                            cmap='gnuplot2', vmin=0, vmax=5)

    def _plot_constant_plane(self, mesh1, mesh2, data, title, xlabel, ylabel):
        self._new_filled_contour(mesh1, mesh2, data)
        self._set_texts(title, xlabel, ylabel)
        self._set_colorbar()
        self._set_axis()

    # FLORIS-specific data manipulation and plotting
    def _initialize_flowfield_for_plotting(self):
        self.flowfield.grid_resolution = self.grid_resolution
        self.flowfield.xmin, self.flowfield.xmax, self.flowfield.ymin, self.flowfield.ymax, self.flowfield.zmin, self.flowfield.zmax = self._set_domain_bounds()
        self.flowfield.x, self.flowfield.y, self.flowfield.z = self._discretize_freestream_domain()
        self.flowfield.initial_flowfield = self.flowfield._initial_flowfield()
        self.flowfield.u_field = self.flowfield._initial_flowfield()
        for turbine in self.flowfield.turbine_map.turbines:
            turbine.plotting = True
        self.flowfield.calculate_wake()

    def _discretize_freestream_domain(self):
        """
            Generate a structured grid for the entire flow field domain.
        """
        x = np.linspace(self.flowfield.xmin, self.flowfield.xmax, self.flowfield.grid_resolution.x)
        y = np.linspace(self.flowfield.ymin, self.flowfield.ymax, self.flowfield.grid_resolution.y)
        z = np.linspace(self.flowfield.zmin, self.flowfield.zmax, self.flowfield.grid_resolution.z)
        return np.meshgrid(x, y, z, indexing="ij")

    def _set_domain_bounds(self):
        coords = self.flowfield.turbine_map.coords
        x = [coord.x for coord in coords]
        y = [coord.y for coord in coords]
        eps = 0.1
        xmin = min(x) - 2 * self.flowfield.max_diameter
        xmax = max(x) + 10 * self.flowfield.max_diameter
        ymin = min(y) - 3 * self.flowfield.max_diameter
        ymax = max(y) + 3 * self.flowfield.max_diameter
        zmin = 0 + eps 
        zmax = 2 * self.flowfield.hub_height
        return xmin, xmax, ymin, ymax, zmin, zmax

    def _add_turbine_marker(self, turbine, coords, wind_direction):
        a = Coordinate(coords.x, coords.y - turbine.rotor_radius)
        b = Coordinate(coords.x, coords.y + turbine.rotor_radius)
        a.rotate_z(turbine.yaw_angle - wind_direction, coords.as_tuple())
        b.rotate_z(turbine.yaw_angle - wind_direction, coords.as_tuple())
        plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1)

    def _plot_constant_z(self, xmesh, ymesh, data):
        self._plot_constant_plane(
            xmesh, ymesh, data, "z plane", "x (m)", "y (m)")

    def _plot_constant_y(self, xmesh, zmesh, data):
        self._plot_constant_plane(
            xmesh, zmesh, data, "y plane", "x (m)", "z (m)")

    def _plot_constant_x(self, ymesh, zmesh, data):
        self._plot_constant_plane(
            ymesh, zmesh, data, "x plane", "y (m)", "z (m)")

    def _add_z_plane(self, percent_height=0.5):
        plane = int(self.flowfield.grid_resolution.z * percent_height)
        self._plot_constant_z(
            self.flowfield.x[:, :, plane],
            self.flowfield.y[:, :, plane],
            self.flowfield.u_field[:, :, plane])
        for coord, turbine in self.flowfield.turbine_map.items():
            self._add_turbine_marker(
                turbine, coord, self.flowfield.wind_direction)

    def _add_y_plane(self, percent_height=0.5):
        plane = int(self.flowfield.grid_resolution.y * percent_height)
        self._plot_constant_y(
            self.flowfield.x[:, plane, :],
            self.flowfield.z[:, plane, :],
            self.flowfield.u_field[:, plane, :])

    def _add_x_plane(self, percent_height=0.5):
        plane = int(self.flowfield.grid_resolution.x * percent_height)
        self._plot_constant_x(
            self.flowfield.y[plane, :, :],
            self.flowfield.z[plane, :, :],
            self.flowfield.u_field[plane, :, :])

    ## plot_z_planes
    #
    # Plots Z planes of the u_field generated from a static floris model
    #
    #   @param planes is a list of numbers between 0-1
    #           which represent a percentage of the
    #           flowfield's height
    #
    # Generates heatmaps of windspeeds on that plane of the windfarm
    #
    def plot_z_planes(self, planes):
        for p in planes:
            self._add_z_plane(p)
        self._show()

    ## \func plot_y_planes
    #
    # \brief Plots Y planes of the u_field generated from a static floris model
    #
    #   @param planes is a list of numbers between 0-1
    #           which represent a percentage of the
    #           flowfield's width
    #
    # Generates heatmaps of windspeeds on those planes of the windfarm
    #
    def plot_y_planes(self, planes):
        for p in planes:
            self._add_y_plane(p)
        self._show()

    ## \func plot_x_planes
    # Plots X planes of the u_field generated from a static floris model
    #
    #   @param planes is a list of numbers between 0-1
    #           which represent a percentage of the
    #           flowfield's length
    #
    # Generates heatmaps of windspeeds on those planes of the windfarm
    #
    def plot_x_planes(self, planes):
        for p in planes:
            self._add_x_plane(p)
        self._show()

    ## Shows plots generated by the Visualization manager
    def _show(self):
        plt.show()

    # Dave's stuff
    ## @func _dirErrorConSpd
    # Calculates error over a u_field at a constant speed
    # over a specified range of direction values
    #
    # @param flowfield takes a Visualization_manager.flowfield object
    # @return page, MAX, MIN, Page, max, min
    #
    #       page[][]:
    #           each page[][0] is an X mesh
    #           each page[][1] is a Y mesh
    #           each page[][2] is an error field at a given step in Drange
    #       MAX is the maximum error value computed
    #       MIN is the minimum error value computed
    #       Page[][]:
    #           each Page[][0] is an X mesh
    #           each Page[][1] is a Y mesh
    #           each Page[][2] is a u_field at a given step in Drange
    #       max is the maximum u_field value computed
    #       min is the minimum u_field value computed
    #
    def _dirErrorConSpd(self, flowfield):

        ref = deepcopy(flowfield)   # A deepcopy of the flowfield to be altered
        ref2 = deepcopy(self)       # A deepcopy of the vman object to hold new plots
        U1 = self.flowfield.u_field # The current visual_manager.flowfield object
        i = self.params.Drange[0]               # A value used to iterate through the values in Drange
        idx = 0                     # A value used to keep track of the index
        # a discretized representation of the u_field at percent_height
        plane = int(self.flowfield.grid_resolution.z * self.params.percent_height)
        # each 'page' contains an X/Y mesh and a u_field error plot
        page = [None] * int((self.params.Drange[1] -
                             self.params.Drange[0]) /
                             self.params.Drange[2] + 1)

        # each Page contains an X/Y mesh and a u_field
        # which corresponds to the error plotted in the same 'page' index
        Page = [None] * int((self.params.Drange[1] -
                             self.params.Drange[0]) /
                             self.params.Drange[2] + 1)
        # MAX holds the largest error value
        MAX = 0
        # MIN holds the minimum error value over the given range
        MIN = 0
        # max holds the maximum u_field value over the given range
        max = 0
        # min holds the minimum u_field value over the given range
        min = 0
        while i <= self.params.Drange[1]:

            # set the new wind direction
            ref.wind_direction = self.flowfield.wind_direction+np.deg2rad(i)
            # a list of the coordinates of the turbines in the field
            turbines = [turbine for _, turbine in ref.turbine_map.items()]
            for k, turbine in enumerate(turbines):
                turbine.yaw_angle = np.deg2rad(i)
            ref.calculate_wake()
            # The reference field from which the error is calculated
            U2 = ref.u_field
            # calculate the error
            ref2.flowfield.u_field = U2-U1

            page[idx]=[ref2.flowfield.x[:, :, plane],
                        ref2.flowfield.y[:, :, plane],
                        ref2.flowfield.u_field[:, :, plane]]
            Page[idx]=[ref.x[:, :, plane],
                        ref.y[:, :, plane],
                        ref.u_field[:, :, plane]]
            if np.amax(page[idx]) > MAX:
                MAX = np.amax(page[idx][2])
            elif np.amin(page[idx]) < MIN:
                MIN = np.amin(page[idx][2])
            if np.amax(Page[idx]) > max:
                max = np.amax(Page[idx][2])
            elif np.amin(Page[idx][2]) < min:
                min = np.amin(Page[idx][2])
            idx += 1
            i += self.params.Drange[2]
        idx -= 1


        return page, MAX, MIN, Page, max, min

    def _grid(self, x, y, z):
        X, Y = np.mgrid[min(x):max(x):self.grid_res[0], min(y):max(y):self.grid_res[1]]
        Z = griddata((x, y), z, (X, Y), method='linear')
        # scipy.interpolate.griddata((x, y), z, (xi, yi), method='cubic')

        # X, Y = np.meshgrid(xi, yi)
        return X, Y, Z
    ## @func animateDnSerror
    # Plots error over a u_field
    # throughout a range of speeds and directions
    #
    # @param self.params.Srange takes a 3-element array
    #               of the following format: [min speed, max speed, step]
    # @param self.params.Drange takes a 3-element array
    #               of the following format: [min direction, max direction, step]
    # @param self.params.percent_height takes a value between zero & one
    #                       which represents a percentage of the
    #                       flowfield's height
    #
    # Generates a plot with two slider bars
    #   which allow you to set different speed/direction combinations
    #   and see the error vs. the original flowfield defined by WF
    #   The plot also includes a button to toggle between error and u_field plots
    #   so you can see the actual plot of the u_field which the error is generated from.
    def animateDnSerror(self):

        Espeed = self.flowfield.wind_speed # estimated wind speed
        Edir = np.rad2deg(self.flowfield.wind_direction) # estimated wind direction
        # title string
        Tstr = 'Flowfield Estimation Error\n(estimated wind direction: ' \
               + str(Edir) + ' degrees)\n' \
               + '(estimated wind speed: ' + str(Espeed) + ' mph)'
        ref = deepcopy(self)    # reference flowfield
        U1 = self.flowfield.u_field # test flowfield
        i = self.params.Srange[0]   # i is used to iterate through the speed range
        errMax = 0 # keep track of the maximum error
        spdMax = 0 # keep track of the maximum speed
        spdMin = 0 # keep track of the minimum speed
        idx = 0 # integer index of current iteration
        # array of error fields
        page = [None] * int((self.params.Srange[1] -
                             self.params.Srange[0]) /
                            self.params.Srange[2] + 1)
        # array of u_fields
        Page = [None] * int((self.params.Srange[1] -
                             self.params.Srange[0]) /
                            self.params.Srange[2] + 1)
        # array of VisualizationManager objects
        vman = [None] * int((self.params.Srange[1] -
                             self.params.Srange[0]) /
                            self.params.Srange[2] + 1)

        while i <= self.params.Srange[1]:
            # set the speed in the JSON object
            self.WF['farm']['properties']['wind_speed'] = i
            # make a new vman object with it
            vman[idx] = VisualizationManager(self.WF, self.grid_res)
            # calculate error fields over the direction range at the current speed
            page[idx], Emax, Emin, Page[idx], Smax, Smin = \
                self._dirErrorConSpd(vman[idx].flowfield)
            # keep track of maximums
            if Emax > errMax:
                errMax = Emax
            if abs(Emin) > errMax:
                errMax = abs(Emin)
            if Smax > spdMax:
                spdMax = Smax
            if Smin < spdMin:
                spdMin = Smin
            # let the user know something is happening
            print("iteration " + str(idx) + " complete")
            # increment i by the speed step
            i += abs(self.params.Srange[2])
            # increment the index
            idx += 1
        # decrement the index since we just left the loop
        idx -= 1
        # print some of the measurements/calculations
        print('max speed: ' + str(spdMax))
        print('min speed: ' + str(spdMin))
        print('max error: ' + str(errMax))

        # start setting up the plot

        # range and resolution of the error plot
        v = np.linspace(errMax*-1, errMax, 100)
        # range and resolution of ticks on the error plot's colorbar
        w = np.linspace(round(errMax)*-1, round(errMax), round(errMax)*4+1, endpoint=True)
        # range and resolution of the u_field plot
        V = np.linspace(spdMin, spdMax, 100)
        # range and resolution of the u_field plot's colorbar
        W = np.linspace(round(spdMin), round(spdMax), round(spdMax) + 1, endpoint=True)
        # put the index in the middle for the speed slider
        idx = round(idx/2)
        # set up an index for the direction slider
        idy = int(abs(self.params.Drange[0])/self.params.Drange[2])
        # start with an initial image
        #im_h = plt.contourf(page[idx][idy][0], page[idx][idy][1], page[idx][idy][2], v,
        #                    cmap='seismic')
        # set plot title
        fontsize = 15


        C_A = plt.gca()
        map = C_A.contourf(page[idx][idy][0], page[idx][idy][1], page[idx][idy][2], v,
                            cmap='seismic')
        plt.suptitle(Tstr, fontsize=fontsize)
        # set up the error colorbar
        cax, _ = mpl.colorbar.make_axes(C_A)
        self.cbr = plt.colorbar(map, cax=cax)
        #cb = mpl.colorbar.ColorbarBase(cax, cmap='seismic')
        self.cbr.set_ticks(w, True)
        self.cbr.ax.tick_params(labelsize=15)
        self.cbr.set_clim(vmin=errMax * -1, vmax=errMax)
        self.cbr.set_label('Error in mph')
        # put the turbines in the plot
        for coord, turbine in ref.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - ref.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - ref.flowfield.wind_direction, coord.as_tuple())
            C_A.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
        # set the location and size of the slider bars
        ax_spd = plt.axes([0.23, 0.02, 0.56, 0.02])
        ax_dir = plt.axes([0.23, 0.04, 0.56, 0.02])
        # set up the speed slider
        slider_spd = Slider(ax_spd,
                            'actual wind speed',
                            self.params.Srange[0],
                            self.params.Srange[1],
                            valinit=Espeed)
        slider_spd.valtext.set_text('{}'.format(Espeed) + ' mph')
        # set up the direction slider
        slider_dir = Slider(ax_dir,
                            'actual wind direction',
                            self.params.Drange[0],
                            self.params.Drange[1],
                            valinit=0.0)
        slider_dir.valtext.set_text('{}'.format(Edir) + ' degrees')
        # set up the error/u_field button
        self._boolBTN = False
        ax_btn = plt.axes([0.01, 0.90, 0.125, 0.05])
        btn = Button(ax_btn, 'toggle\nspeed/error')
        # this is what happens when the button is clicked
        def toggle_plot(event):
            self._boolBTN = not self._boolBTN
            update_plot(event)
        btn.on_clicked(toggle_plot)
        # this gets called when a slider is moved or the button is clicked
        def update_plot(val):
            # set the speed index based on the slider value
            idx = int(round((slider_spd.val -
                             self.params.Srange[0]) /
                            self.params.Srange[2]))
            # set the new speed value
            spdvalue = self.params.Srange[0] + idx * self.params.Srange[2]
            # put the speed value in the slider text
            slider_spd.valtext.set_text('{}'.format(spdvalue) + ' mph')
            # set the direction index based on the slider value
            idy = int(round((slider_dir.val -
                             self.params.Drange[0]) /
                            self.params.Drange[2]))
            # set the new direction value
            dirvalue = idy * self.params.Drange[2] + self.params.Drange[0]
            # put the direction value in the slider text
            slider_dir.valtext.set_text('{}'.format(Edir + dirvalue) + ' degrees')
            # clear the current plot
            C_A.clear()
            #self.cbr.remove()
            # plot the turbines
            for coord, turbine in self.flowfield.turbine_map.items():
                a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
                b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
                a.rotate_z(turbine.yaw_angle - ref.flowfield.wind_direction, coord.as_tuple())
                b.rotate_z(turbine.yaw_angle - ref.flowfield.wind_direction, coord.as_tuple())
                C_A.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
            # set up the plot for either error or u_field plot depending on the button
            if not self._boolBTN:
                # plot the error data
                map = C_A.contourf(page[idx][idy][0], page[idx][idy][1], page[idx][idy][2], v,
                         cmap='seismic', vmin=errMax*-1,vmax=errMax)
                # set the colorbar
                self.cbr = plt.colorbar(map, cax=cax)
                self.cbr.set_clim(vmin=errMax*-1,vmax=errMax)
                self.cbr.set_ticks(w)
                self.cbr.set_cmap('seismic')
                self.cbr.set_label('Error in mph')
                self.cbr.draw_all()
            else:
                # plot the u_field data
                map = C_A.contourf(Page[idx][idy][0], Page[idx][idy][1], Page[idx][idy][2], V,
                         cmap='gnuplot2', vmin=spdMin,vmax=spdMax)
                # set the colorbar
                self.cbr = plt.colorbar(map, cax=cax)
                self.cbr.set_clim(vmin=spdMin,vmax=spdMax)
                self.cbr.set_ticks(W)
                self.cbr.set_cmap('gnuplot2')
                self.cbr.set_label('Speed in mph')
                self.cbr.draw_all()
            # put the title back up
            plt.suptitle(Tstr, fontsize=fontsize)
            # draw the plot
            plt.draw()

        # these set up the callbacks for the sliders
        slider_spd.on_changed(update_plot)
        slider_dir.on_changed(update_plot)
        plt.show()

    ## @func plotSensitivityMatrix
    # Plots the sensitivity matrix of a given u_field
    #
    # @param self.params.v0 takes a wind speed in meters per second
    # @param self.params.d0 takes a wind direction in degrees
    # @param self.params.epSpeed is an epsilon over which to calculate df/dv (derivative of speed)
    # @param self.params.epDir is an epsilon over which to calculate df/dd (derivative of direction)
    # @param self.params.percent_height takes a value between zero & one
    #                       which represents a percentage of the
    #                       flowfield's height
    #
    # Generates a plot with four heatmaps derived from the pseudo-invers of the jacobian:
    #   1. The partial derivative of the floris model with respect to direction (df/dd)
    #   2. The partial derivative of the floris model with respect to speed (df/dv)
    #   3. normalized ||df/dd||+||df/dv||
    #   4. normalized ||df/dd||*||df/dv||
    #
    def plotSensitivityMatrix(self):

        # set wind speed estimate for floris model
        self.WF['farm']['properties']['wind_speed'] = self.params.v0
        # calculate f(v0,d0) = ff
        ff_vman = VisualizationManager(self.WF, self.grid_res)
        ff_vman.flowfield.wind_direction = self.params.d0
        ff_vman.flowfield.calculate_wake()
        ff = deepcopy(ff_vman.flowfield.u_field)

        # print v0,d0 as initial speed and direction estimates
        print('v0 = ' + str(self.params.v0))
        print('d0 = ' + str(self.params.d0))

        # calculate f(v0+ev,d0) = ff_v1
        self.WF['farm']['properties']['wind_speed'] = self.params.v0 + self.params.epSpeed  # add speed epsilon
        ff_v1_vman = VisualizationManager(self.WF, self.grid_res)
        ff_v1_vman.flowfield.wind_direction = self.params.d0
        ff_v1_vman.flowfield.calculate_wake()
        ff_v1 = deepcopy(ff_v1_vman.flowfield.u_field)

        # calculate f(v0,d0+ed) = ff_d1
        ff_d1_vman = deepcopy(ff_vman)
        ff_d1_vman.flowfield.wind_direction = self.params.d0 + self.params.epDir  # add direction epsilon
        ff_d1_vman.flowfield.calculate_wake()  # recalculate wake
        ff_d1 = deepcopy(ff_d1_vman.flowfield.u_field)

        plane = int(ff_d1_vman.flowfield.grid_resolution.z * self.params.percent_height)

        # calculate gradient of f(v,d) = grad_f
        df_dv = (ff_v1 - ff) / self.params.epSpeed  # partial of f WRT speed
        df_dd = (ff_d1 - ff) / self.params.epDir  # partial of f WRT direction
        dfdv = deepcopy(df_dv[:, :, plane])
        dfdd = deepcopy(df_dd[:, :, plane])

        grad_f = np.column_stack([dfdv.flatten(), dfdd.flatten()])  # 2xn gradient of f (jacobian)
        grad_f_pinv = np.linalg.pinv(grad_f)
        # print('grad_f dim: '+str(len(grad_f_pinv))+' rows x '
        #       +str(len(grad_f_pinv[0]))+' cols')

        f, axarr = plt.subplots(2, 2, sharex='col')
        f.suptitle('heatmap of rows from pseudoinverse of sensitivity matrix')

        X = deepcopy(ff_d1_vman.flowfield.x[:, :, plane])
        Y = deepcopy(ff_d1_vman.flowfield.y[:, :, plane])
        # print('dim X: '+str(len(X.flatten()))+' rows x '+str(len(X.flatten()[0]))+'cols')
        dMax = np.amax(np.abs(df_dd[:, :, plane]))
        print('dMax: ' + str(dMax))
        v = np.linspace(0, 1, 100)
        V = np.linspace(0, 1, 5)



        # cf = axarr[0].contourf(ff_d1_vman.flowfield.x[:,:,plane],
        #                     ff_d1_vman.flowfield.y[:, :, plane],
        #                     np.abs(df_dd[:,:,plane])/dMax, v, cmap='gnuplot2')
        Zd = abs(grad_f_pinv[1].flatten()) / np.amax(abs(grad_f_pinv[1]))
        x, y, z = self._grid(X.flatten(), Y.flatten(), Zd)
        cf = axarr[0][0].contourf(x, y, z, v, cmap='gnuplot2')
        # cf = axarr[0].scatter(X.flatten(),Y.flatten(),abs(grad_f_pinv[1].flatten())/
        #                      np.amax(abs(grad_f_pinv[1].flatten())),
        #                      c=abs(grad_f_pinv[1].flatten())/
        #                      np.amax(abs(grad_f_pinv[1].flatten())),
        #                      cmap='gnuplot2')

        for coord, turbine in ff_d1_vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - ff_d1_vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - ff_d1_vman.flowfield.wind_direction, coord.as_tuple())
            axarr[0][0].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')
        cbr = f.colorbar(cf, ax=axarr[0][0])
        cb = mpl.colorbar.ColorbarBase(cbr.ax, cmap='gnuplot2')
        cb.set_clim(vmin=0, vmax=1)
        cb.set_ticks(V)
        cb.set_label('||df/dd||')
        cb.draw_all()

        vMax = np.amax(np.abs(df_dv[:, :, plane]))
        print('vMax: ' + str(vMax))
        w = np.linspace(0, 1, 100)
        W = np.linspace(0, 1, 5)
        # cf2 = axarr[1].contourf(ff_d1_vman.flowfield.x[:,:,plane],
        #                       ff_d1_vman.flowfield.y[:,:,plane],
        #                       np.abs(df_dv[:,:,plane])/vMax, w, cmap='gnuplot2')
        Zv = abs(grad_f_pinv[0].flatten()) / np.amax(abs(grad_f_pinv[0]))
        x, y, z = self._grid(X.flatten(), Y.flatten(), Zv)
        cf2 = axarr[0][1].contourf(x, y, z, w, cmap='gnuplot2')

        for coord, turbine in ff_d1_vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - ff_d1_vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - ff_d1_vman.flowfield.wind_direction, coord.as_tuple())
            axarr[0][1].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')
        cbr2 = f.colorbar(cf2, ax=axarr[0][1])
        cb2 = mpl.colorbar.ColorbarBase(cbr2.ax, cmap='gnuplot2')
        cb2.set_ticks(W)
        cb2.set_label('||df/dv||')
        cb2.draw_all()

        uMax = np.amax(Zd + Zv)
        print('uMax: ' + str(uMax))
        u = np.linspace(0, round(uMax), 100)
        U = np.linspace(0, round(uMax), 5)
        # cf3 = axarr[2].contourf(ff_d1_vman.flowfield.x[:,:,plane],
        #                       ff_d1_vman.flowfield.y[:,:,plane],
        #                       np.abs(df_dv[:,:,plane])/vMax+np.abs(df_dd[:,:,plane])/dMax,
        #                        u, cmap='gnuplot2')

        x, y, z = self._grid(X.flatten(), Y.flatten(), Zd + Zv)
        cf3 = axarr[1][0].contourf(x, y, z, u, cmap='gnuplot2')

        for coord, turbine in ff_d1_vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - ff_d1_vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - ff_d1_vman.flowfield.wind_direction, coord.as_tuple())
            axarr[1][0].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')
        cbr3 = f.colorbar(cf3, ax=axarr[1][0])
        cb3 = mpl.colorbar.ColorbarBase(cbr3.ax, cmap='gnuplot2')
        cb3.set_clim(vmin=0, vmax=round(uMax))
        cb3.set_ticks(U)
        cb3.set_label('||df/dv||+||df/dd||')
        cb3.draw_all()

        nMax = np.amax(Zd * Zv)
        print('nMax: ' + str(nMax))
        n = np.linspace(0, round(nMax), 100)
        N = np.linspace(0, round(nMax), 5)

        x, y, z = self._grid(X.flatten(), Y.flatten(), Zd * Zv)
        cf4 = axarr[1][1].contourf(x, y, z, n, cmap='gnuplot2')

        for coord, turbine in ff_d1_vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - ff_d1_vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - ff_d1_vman.flowfield.wind_direction, coord.as_tuple())
            axarr[1][1].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')
        cbr4 = f.colorbar(cf4, ax=axarr[1][1])
        cb4 = mpl.colorbar.ColorbarBase(cbr4.ax, cmap='gnuplot2')
        cb4.set_clim(vmin=0, vmax=round(nMax))
        cb4.set_ticks(N)
        cb4.set_label('||df/dv||*||df/dd||')
        cb4.draw_all()
        plt.show()


    ## @func plotSensitivityMatrix
    # Plots the convergence of speed and direction estimates
    # based on a reduced sensitivity matrix
    #v0, d0, vBar, dBar, epSpeed, epDir, spErrMax, dirErrMax, iterMax, mask_thresh
    # @param self.params.v0 takes a wind speed in meters per second
    # @param self.params.d0 takes a wind direction in degrees
    # @param self.params.vBar takes a wind speed in meters per second
    # @param self.params.dBar takes a wind direction in degrees
    # @param self.params.epSpeed is an epsilon over which to calculate df/dv (derivative of speed)
    # @param self.params.epDir is an epsilon over which to calculate df/dd (derivative of direction)
    # @param self.params.spErrMax is a speed error threshold for stopping iterations
    # @param self.params.dirErrMax is a direction error threshold for stopping iterations
    # @param self.params.iterMax is the maximum number of iterations to complete
    # @param self.params.mask_thresh threshold value for the normalized sensitivity matrix
    #
    # Generates a plot
    #
    def reducedSM(self):
        # calculate f(vbar,dbar) = xBar
        print('vBar = ' + str(self.params.vBar))
        print('dBar = ' + str(self.params.dBar))
        vdBar = [[self.params.vBar], [self.params.dBar]]
        self.WF['farm']['properties']['wind_speed'] = self.params.vBar
        ff_bar_vman = VisualizationManager(self.WF, self.grid_res)
        ff_bar_vman.flowfield.wind_direction = self.params.dBar
        ff_bar_vman.flowfield.calculate_wake()
        plane = int(ff_bar_vman.flowfield.grid_resolution.z * self.params.percent_height)
        X = deepcopy(ff_bar_vman.flowfield.x)
        Y = deepcopy(ff_bar_vman.flowfield.y)
        ff_bar = deepcopy(ff_bar_vman.flowfield.u_field)
        ff = deepcopy(ff_bar)
        xBar = np.vstack(ff.flatten())  # 1xn

        # initialize error
        vdErr = [[self.params.vBar - self.params.v0],
                 [self.params.dBar - self.params.d0]]  # initial error
        speedError = list()  # list for history of speed error
        directionError = list()  # list for history of direction error
        print('Initial speed error = ' + str(vdErr[0]))
        print('Initial direction error = ' + str(vdErr[1]))
        speedError.append(vdErr[0][0])
        directionError.append(vdErr[1][0])

        # set up first iteration
        iters = 0
        iterations = list()  # iteration list for graphing
        iterations.append(iters)
        V_k = list()  # list for history of vk
        D_k = list()  # list for history of dk
        V_k.append(self.params.v0)
        D_k.append(self.params.d0)
        vd_k = [[self.params.v0], [self.params.d0]]  # initial estimate
        err_field = list()
        masks = list()
        errMax = 0

        while iters < self.params.iterMax and \
                ( abs(speedError[iters]) > self.params.spErrMax or \
                  abs(directionError[iters]) > self.params.dirErrMax ):
            # set new wind speed estimate for floris model
            self.WF['farm']['properties']['wind_speed'] = vd_k[0][0]
            # calculate f(vk,dk) = temp_x_k
            temp0_vman = VisualizationManager(self.WF, self.grid_res)
            temp0_vman.flowfield.wind_direction = vd_k[1][0]  # new wind direction estimate
            temp0_vman.flowfield.calculate_wake()  # recalculate field for new wind direction
            temp_x_k = deepcopy(temp0_vman.flowfield.u_field)
            err_field.append(abs(ff_bar[:,:,plane]-temp_x_k[:,:,plane]))
            if np.amax(err_field[iters]) > errMax:
                errMax = np.amax(err_field[iters])
            # print vk,dk as new speed and direction estimates
            print('vk = ' + str(vd_k[0][0]))
            print('dk = ' + str(vd_k[1][0]))

            # calculate f(vk+ev,dk) = temp_x_ev
            self.WF['farm']['properties']['wind_speed'] = vd_k[0][0] + self.params.epSpeed  # add speed epsilon to current estimate
            temp1_vman = VisualizationManager(self.WF, self.grid_res)
            temp1_vman.flowfield.wind_direction = vd_k[1][0]  # set wind direction to current estimate
            temp1_vman.flowfield.calculate_wake()
            temp_x_ev = deepcopy(temp1_vman.flowfield.u_field)

            # calculate f(vk,dk+ed) = temp_x_ed
            temp2_vman = deepcopy(temp0_vman)  # we can use same flowfield from temp0
            temp2_vman.flowfield.wind_direction = vd_k[1][0] + self.params.epDir  # set wind direction to current estimate + epsilon
            temp2_vman.flowfield.calculate_wake()  # recalculate wake for new wind speed
            temp_x_ed = deepcopy(temp2_vman.flowfield.u_field)

            # calculate gradient
            temp_df_dv = (temp_x_ev - temp_x_k) / self.params.epSpeed  # partial of f WRT speed
            temp_df_dd = (temp_x_ed - temp_x_k) / self.params.epDir  # partial of f WRT direction
            temp_grad_fk = np.column_stack([temp_df_dv.flatten(), temp_df_dd.flatten()])  # 2xn gradient of f (jacobian)
            temp_sens_mat_pinv = np.linalg.pinv(temp_grad_fk)

            # calculate mask
            temp_dMax = np.amax(np.abs(temp_sens_mat_pinv[1]))
            temp_vMax = np.amax(np.abs(temp_sens_mat_pinv[0]))
            temp_Zd = np.abs(temp_sens_mat_pinv[1] / temp_dMax)
            temp_Zv = np.abs(temp_sens_mat_pinv[0] / temp_vMax)
            temp_Z_mask = temp_Zd * temp_Zv / np.amax(temp_Zd * temp_Zv)
            # print(np.amax(temp_Z_mask))
            temp_Z_mask = np.where(temp_Z_mask > self.params.mask_thresh, 1, 0)
            # print(temp_Z_mask)
            temp_sens_mat_pinv[0] = temp_sens_mat_pinv[0] * temp_Z_mask
            temp_sens_mat_pinv[1] = temp_sens_mat_pinv[1] * temp_Z_mask
            masks.append(temp_Z_mask)

            # calculate pseudoinverse[gradient{f(vk,dk)}]*{xBar-f(vk,dk)} = adj_vd (adjustment to current v&d estimates)
            adj_vd = np.matmul(temp_sens_mat_pinv, xBar - np.vstack(temp_x_k.flatten()))
            # print(adj_vd)
            # calculate v_k+1 and d_k+1
            vd_kp1 = vd_k + adj_vd
            # update vd_k for next iteration
            vd_k = deepcopy(vd_kp1)
            # actually starting loop with 2nd iteration
            iters = iters + 1
            iterations.append(iters)
            print('\n\niteration ' + str(iters) + ' complete.')

            # calculate error = [[vbar],[dbar]]-[[vk],[dk]]
            vdErr = vdBar - vd_k
            print('Speed error: ' + str(vdErr[0]))
            print('Direction error: ' + str(vdErr[1]))
            speedError.append(vdErr[0][0])
            directionError.append(vdErr[1][0])
            V_k.append(vd_k[0][0])
            D_k.append(vd_k[1][0])

            # delete temporary objects
            del temp0_vman
            del temp_x_k
            del temp1_vman
            del temp_x_ev
            del temp2_vman
            del temp_x_ed
            del temp_df_dv
            del temp_df_dd
            del temp_grad_fk
            del adj_vd
            del vd_kp1
        print('Max error = '+str(errMax))
        print('Total iterations: ' + str(iters))

        f, axarr = plt.subplots(2, 3, sharex='col')
        sld_ax = plt.axes([0.23, 0.02, 0.56, 0.02])
        sld = Slider(sld_ax,
                     'iterations',
                     0,iters-1,valinit=0)
        sld.valtext.set_text('iteration 0')
        f.suptitle('Speed estimate: ' + str(self.params.v0) + ' m/s, Actual: ' + str(self.params.vBar) +
                   ' m/s\nDirection estimate: ' + str(np.rad2deg(self.params.d0)) + '\N{DEGREE SIGN}' +
                   ', Actual: ' + str(np.rad2deg(self.params.dBar)) + '\N{DEGREE SIGN}' +
                   '\nThreshold: ' + str(self.params.mask_thresh) + ' Iterations: ' + str(iters) +
                   '\nFinal error: e\N{GREEK SMALL LETTER THETA} = ' + str(vdErr[1][0]) +
                   ' ev = ' + str(vdErr[0][0]))
        axarr[0][0].plot(iterations, V_k)
        axarr[0][0].set_title('speed estimate')
        axarr[0][1].plot(iterations, speedError)
        axarr[0][1].set_title('speed error')
        v = np.linspace(0, errMax, 100)
        V = np.linspace(0, errMax, 5)
        cont = axarr[0][2].contourf(X[:,:,plane], Y[:,:,plane], err_field[0],v,cmap='gnuplot2')
        cb = plt.colorbar(cont, ax=axarr[0][2])
        cb.set_clim(vmin=0, vmax=errMax)
        cb.set_ticks(V, True)
        cb.set_label('u_field error')
        cb.draw_all()
        axarr[0][2].set_title('u_field error')
        axarr[1][0].plot(iterations, np.rad2deg(D_k))
        axarr[1][0].set_title('direction estimate')
        axarr[1][0].set_xlabel('iterations')
        axarr[1][1].plot(iterations, np.rad2deg(directionError))
        axarr[1][1].set_title('direction error')
        axarr[1][1].set_xlabel('iterations')
        x, y = deepcopy(X), deepcopy(Y)
        scat = axarr[1][2].scatter(x.flatten(),y.flatten(),masks[0],color='black')
        cb2 = plt.colorbar(mappable=cont, ax=axarr[1][2])
        cb2.set_clim(vmin=0, vmax=errMax)
        cb2.set_ticks(V)
        cb2.draw_all()
        cb2.remove()
        for coord, turbine in self.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.flowfield.wind_direction, coord.as_tuple())
            axarr[0][2].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='white')
            axarr[1][2].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
        axarr[1][2].set_title('mask')

        # f.tight_layout()
        plt.subplots_adjust(left=0.05,
                            bottom=0.15,
                            right=0.95,
                            top=0.83,
                            wspace=0.27,
                            hspace=0.19)
        def update_plot(val):

            idx = int(round(sld.val))
            sld.valtext.set_text('iteration '+'{}'.format(idx))
            axarr[0][2].clear()
            axarr[0][2].contourf(X[:, :, plane], Y[:, :, plane], err_field[idx], v, cmap='gnuplot2')
            axarr[1][2].clear()
            axarr[1][2].scatter(x,y,masks[idx],color='black')
            for coord, turbine in self.flowfield.turbine_map.items():
                a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
                b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
                a.rotate_z(turbine.yaw_angle - self.flowfield.wind_direction, coord.as_tuple())
                b.rotate_z(turbine.yaw_angle - self.flowfield.wind_direction, coord.as_tuple())
                axarr[0][2].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='white')
                axarr[1][2].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
            plt.draw()
        sld.on_changed(update_plot)
        plt.show()

    ## params is a class used to pass parameters to various functions
    #
    #   animateDnSerror: Srange, Drange, (percent_height)
    #
    #   plotSensitivityMatrix: v0, d0,j epSpeed, epDir, (percent_height)
    #
    #   reducedSM: v0, d0, vBar, dBar, epSpeed, epDir, spErrMax, dirErrMax, iterMax, mask_thresh
    #
    class params:
        ## v0 serves as the initial speed estimate
        v0 = 8.0  # initial speed estimate
        ## d0 serves as the initial direction estimate
        d0 = np.deg2rad(0.0)  # initial direction estimate
        ## speed epsilon for calculating the sensitivity matrix
        epSpeed = 0.001  # speed epsilon (ev)
        ## direction epsilon for calculating the sensitivity matrix
        epDir = 0.0001  # direction epsilon (ed)
        ## speed error threshold for stopping iterations
        spErrMax = 0.1  # speed error threshold
        ## direction error threshold for stopping iterations
        dirErrMax = 0.01  # direction error threshold
        ## iteration threshold for stopping iterations
        iterMax = 40  # iteration threshold
        ## vBar serves as the 'actual' wind speed
        vBar = 15.0 # actual wind speed
        ## dBar serves as the 'actual' wind direction
        dBar = np.deg2rad(10.0) # actual wind direction
        ## eliminate values in the normalized sensitivity matrix lower than 0.35
        mask_thresh = 0.35 # masking threshold for normalized sensitivity matrix
        ## a value between zero & one
        #  which represents a percentage of the
        #  flowfield's height
        percent_height = 0.5 # percentage of flowfield height
        ## average wind speed in meters per second
        speed = 8.0  # average wind speed
        ## a range of speeds: [min speed, max speed, step]
        Srange = [speed - 2.0, speed + 2.0, 0.25]  # {min speed, max speed, step}
        ## a range of directions: [min direction err, max direction err, step]
        Drange = [-2.0, 2.0, 0.25]  # [min dir error, max dir error, step]
        ## a JSON windfarm object
