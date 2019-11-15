## \file VIZMAN_new.py
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
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Slider, Button
import numpy as np
from scipy.interpolate import griddata
import scipy.io as scio
from copy import deepcopy
from matplotlib import animation as anim
from old.readVTK import VTKreader

## \class VisualizationManager
#
#    1. animateDnSerror
#        Plot error against a base FLORIS model
#        compared to a range of different directions and speeds.
#
#    2. calcSensitivityMatrix
#        calculate the Sensitivity matrix of a Z plane
#        of the floris model at hub height.
#
#    3. reducedSM
#        Plotting the convergence of speed and direction estimates
#        with a reduced sensitivity matrix
#
class VisualizationManager():

    ## Class constructor
    #   @param WF takes a JSON windfarm object
    #   @param grid_resolution takes a 3-element integer array
    #           with [x_resolution, y_resolution_ z_resolution]
    #   @return VisualizationManager object
    def __init__(self, WF, grid_resolution=(100, 100, 25)):
        ## @var WF
        #   A JSON object used to define the wind farm
        self.WF = WF
        ## @var flowfield
        #   holds an instance of the flowfield instantiated from WF
        self.u_field = WF.u_field
        ## @var grid_resolution
        # holds a copy of the grid_resolution passed to the constructor
        self.grid_resolution = grid_resolution
        self._boolBTN = False

    # Dave's stuff

    # animateDnSerror
    ## @brief Plots error over a u_field
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

        Espeed = self.WF.wind_speed # estimated wind speed
        Edir = self.WF.wind_direction # estimated wind direction
        self.WF.set_incoming(Espeed, Edir)
        Uref = deepcopy(self.WF.u_field)
        Smax = np.amax(Uref)
        Smin = np.amin(Uref)
        # title string
        Tstr = 'Flowfield Estimation Error\n(estimated wind direction: ' \
               + str(Edir) + ' degrees)\n' \
               + '(estimated wind speed: ' + str(Espeed) + ' mph)'

        errMax = 0 # keep track of the maximum error
        spdMax = 0 # keep track of the maximum speed
        spdMin = 0 # keep track of the minimum speed
        idx = 0 # integer index of current speed iteration
        idy = 0 # integer index of current direction iteration
        # array of error fields
        idx_max = int((self.params.Srange[1] -
                             self.params.Srange[0]) /
                            self.params.Srange[2] + 1)
        idy_max = int((self.params.Drange[1] -
                       self.params.Drange[0]) /
                      self.params.Drange[2] + 1)

        # array of arrays of error fields
        page = [None] * idx_max
        # array of arrays of u_fields
        Page = [None] * idx_max
        i = self.params.Srange[0]  # i is used to iterate through the speed range
        while i <= self.params.Srange[1]:
            idy = 0
            j = self.params.Drange[0]
            temp_page = [None] * idy_max
            Temp_Page = [None] * idy_max

            while j <= self.params.Drange[1]:
                # set the speed in the JSON object
                self.WF.set_incoming(i,j)
                U = deepcopy(self.WF.u_field)
                Uerr = U-Uref

                temp_page[idy]=deepcopy(Uerr)
                Temp_Page[idy]=deepcopy(U)

                # calculate error fields over the direction range at the current speed
                Emax = np.amax(Uerr)
                Emin = np.amin(Uerr)
                # keep track of maximums
                if Emax > errMax:
                    errMax = Emax
                if abs(Emin) > errMax:
                    errMax = abs(Emin)
                if Smax > spdMax:
                    spdMax = Smax
                if Smin < spdMin:
                    spdMin = Smin



                # increment j by the direction step
                j += abs(self.params.Drange[2])
                print("direction iteration "+str(idy+1)+" of "+str(idy_max)+" complete.")
                idy += 1
            # increment i by the speed step
            i += abs(self.params.Srange[2])
                # increment the index
            print("speed iteration " + str(idx+1) +" of "+str(idx_max)+" complete")
            page[idx]=deepcopy(temp_page)
            Page[idx]=deepcopy(Temp_Page)
            idx += 1
        # print some of the measurements/calculations
        print('max speed: ' + str(spdMax))
        print('min speed: ' + str(spdMin))
        print('max error: ' + str(errMax))
        idx -= 1
        idy -= 1
        # start setting up the plot

        # range and resolution of the error plot
        v = np.linspace(errMax*-1, errMax, 100)
        # range and resolution of ticks on the error plot's colorbar
        w = np.linspace(round(errMax)*-1, round(errMax), 3, endpoint=True)
        # range and resolution of the u_field plot
        V = np.linspace(spdMin, spdMax, 100)
        # range and resolution of the u_field plot's colorbar
        W = np.linspace(round(spdMin), round(spdMax), 5, endpoint=True)
        # put the index in the middle for the speed slider
        idx = round(idx/2)
        # set up an index for the direction slider
        idy = round(idy/2)
        # set font size
        fontsize = 14


        C_A = plt.gca()
        map = C_A.contourf(self.WF.x_mesh, self.WF.y_mesh, page[idx][idy], v,
                            cmap='seismic')
        plt.suptitle(Tstr, fontsize=fontsize)
        C_A.set_aspect('equal')
        C_A.text(-220, 290, '8.0 m/s, 0.0\N{DEGREE SIGN}', fontsize=fontsize)
        C_A.arrow(-220, 225, 15 * 8 * np.cos(np.deg2rad(0)),
                  -15 * 8 * np.sin(np.deg2rad(0)), head_width=50, head_length=50)
        C_A.set_xlabel('meters', fontsize=fontsize)
        C_A.set_ylabel('meters', fontsize=fontsize)
        #divider = make_axes_locatable(C_A)
        #cax = divider.append_axes("right", size="5%", pad=0.05)

        cax, _ = mpl.colorbar.make_axes(C_A)
        ## @var cbr
        # holds a copy of the colorbar
        # plt.colorbar(im, fraction=0.046, pad=0.04)
        self._cbr = plt.colorbar(map, cax=cax)
        #cb = mpl.colorbar.ColorbarBase(cax, cmap='seismic')
        self._cbr.set_ticks(w, True)
        self._cbr.ax.tick_params(labelsize=14)
        self._cbr.set_clim(vmin=errMax * -1, vmax=errMax)
        self._cbr.set_label('Error in m/s', fontsize=fontsize)
        # put the turbines in the plot
        for coord, turbine in self.WF.turbines:
            x_0 = coord.x1 + np.sin(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            x_1 = coord.x1 - np.sin(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            y_0 = coord.x2 - np.cos(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            y_1 = coord.x2 + np.cos(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            C_A.plot([x_0, x_1], [y_0, y_1], color='k', linewidth=1)

        for tick in C_A.xaxis.get_major_ticks():
            tick.label.set_fontsize(14)
        for tick in C_A.yaxis.get_major_ticks():
            tick.label.set_fontsize(14)
            # set the location and size of the slider bars
        ax_spd = plt.axes([0.23, 0.02, 0.56, 0.02])
        ax_dir = plt.axes([0.23, 0.04, 0.56, 0.02])
        # set up the speed slider
        slider_spd = Slider(ax_spd,
                            'actual wind speed',
                            self.params.Srange[0],
                            self.params.Srange[1],
                            valinit=Espeed)
        slider_spd.valtext.set_text('{}'.format(Espeed) + ' m/s')
        # set up the direction slider
        slider_dir = Slider(ax_dir,
                            'actual wind direction',
                            self.params.Drange[0],
                            self.params.Drange[1],
                            valinit=Edir)
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
            slider_spd.valtext.set_text('{}'.format(spdvalue) + ' m/s')
            # set the direction index based on the slider value
            idy = int(round((slider_dir.val -
                             self.params.Drange[0]) /
                            self.params.Drange[2]))
            # set the new direction value
            dirvalue = idy * self.params.Drange[2] + self.params.Drange[0]
            # put the direction value in the slider text
            slider_dir.valtext.set_text('{}'.format(dirvalue) + ' degrees')
            # clear the current plot
            C_A.clear()
            #self._cbr.remove()
            # plot the turbines
            for coord, turbine in self.WF.turbines:
                x_0 = coord.x1 + np.sin(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
                x_1 = coord.x1 - np.sin(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
                y_0 = coord.x2 - np.cos(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
                y_1 = coord.x2 + np.cos(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
                C_A.plot([x_0, x_1], [y_0, y_1], color='k', linewidth=1)

            # set up the plot for either error or u_field plot depending on the button
            if not self._boolBTN:
                # plot the error data
                map = C_A.contourf(self.WF.x_mesh, self.WF.y_mesh, page[idx][idy], v,
                         cmap='seismic', vmin=errMax*-1,vmax=errMax)

                # set the colorbar
                self._cbr = plt.colorbar(map, cax=cax)
                self._cbr.set_clim(vmin=errMax*-1,vmax=errMax)
                self._cbr.set_ticks(w)
                self._cbr.set_cmap('seismic')
                self._cbr.set_label('Error in m/s', fontsize=fontsize)
                self._cbr.draw_all()
            else:
                # plot the u_field data
                map = C_A.contourf(self.WF.x_mesh, self.WF.y_mesh, Page[idx][idy], V,
                         cmap='gnuplot2', vmin=spdMin,vmax=spdMax)
                # set the colorbar
                self._cbr = plt.colorbar(map, cax=cax)
                self._cbr.set_clim(vmin=spdMin,vmax=spdMax)
                self._cbr.set_ticks(W)
                self._cbr.set_cmap('gnuplot2')
                self._cbr.set_label('Speed in m/s')
                self._cbr.draw_all()
            for tick in C_A.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            for tick in C_A.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            # put the title back up
            C_A.set_xlabel('meters', fontsize=fontsize)
            C_A.set_ylabel('meters', fontsize=fontsize)
            C_A.text(-220,290,str(spdvalue)+' m/s, '+str(dirvalue)+'\N{DEGREE SIGN}',fontsize=fontsize)
            C_A.arrow(-220,225,15*spdvalue*np.cos(np.deg2rad(dirvalue)),-15*spdvalue*np.sin(np.deg2rad(dirvalue)), head_width=50, head_length=50)
            plt.suptitle(Tstr, fontsize=fontsize)
            # draw the plot
            plt.draw()

        # these set up the callbacks for the sliders
        slider_spd.on_changed(update_plot)
        slider_dir.on_changed(update_plot)
        plt.show()

    # calcSensitivityMatrix
    ## @brief calculates the sensitivity matrix of a given u_field
    #
    # @param self.params.v0 takes a wind speed in meters per second
    # @param self.params.d0 takes a wind direction in degrees
    # @param self.params.epSpeed is an epsilon over which to calculate df/dv (derivative of speed)
    # @param self.params.epDir is an epsilon over which to calculate df/dd (derivative of direction)
    # @param self.params.percent_height takes a value between zero & one
    #                       which represents a percentage of the
    #                       flowfield's height
    # @return returns the sensitivity matrix of the u_field
    def calcSensitivityMatrix(self, v0, d0):
        # set wind speed estimate for floris model
        self.WF.set_incoming(v0, np.rad2deg(d0))
        # calculate f(v0,d0) = ff
        ff = deepcopy(self.WF.u_field)

        # print v0,d0 as initial speed and direction estimates
        # print('v0 = ' + str(self.params.v0))
        # print('d0 = ' + str(self.params.d0))

        # calculate f(v0+ev,d0) = ff_v1
        self.WF.set_incoming(v0 + self.params.epSpeed, np.rad2deg(d0))   # add speed epsilon
        ff_v1 = deepcopy(self.WF.u_field)

        # calculate f(v0,d0+ed) = ff_d1
        self.WF.set_incoming(v0,np.rad2deg(d0 + self.params.epDir))  # add direction epsilon
        ff_d1 = deepcopy(self.WF.u_field)

        # calculate gradient of f(v,d) = grad_f
        dfdv = (ff_v1 - ff) / self.params.epSpeed  # partial of f WRT speed
        dfdd = (ff_d1 - ff) / self.params.epDir  # partial of f WRT direction

        sens_mat = [dfdv, dfdd]  # 2xn gradient of f (jacobian)

        return sens_mat

    def calcSensitivityMatrix_vy(self, v0, y0):
        # set wind speed estimate for floris model
        self.WF.set_vy(v0, np.rad2deg(y0))
        # calculate f(v0,d0) = ff
        ff = deepcopy(self.WF.u_field)

        # print v0,d0 as initial speed and direction estimates
        # print('v0 = ' + str(self.params.v0))
        # print('d0 = ' + str(self.params.d0))

        # calculate f(v0+ev,d0) = ff_v1
        self.WF.set_vy(v0 + self.params.epSpeed, np.rad2deg(y0))   # add speed epsilon
        ff_v1 = deepcopy(self.WF.u_field)

        # calculate f(v0,d0+ed) = ff_d1
        self.WF.set_vy(v0,np.rad2deg(y0 + self.params.epYaw))  # add direction epsilon
        ff_d1 = deepcopy(self.WF.u_field)

        # calculate gradient of f(v,d) = grad_f
        dfdv = (ff_v1 - ff) / self.params.epSpeed  # partial of f WRT speed
        dfdy = (ff_d1 - ff) / self.params.epYaw  # partial of f WRT direction

        sens_mat = [dfdv, dfdy]  # 2xn gradient of f (jacobian)

        return sens_mat

    # reducedSM
    ## @brief Plots the convergence of speed and direction estimates
    # based on a reduced sensitivity matrix
    #
    # @param ANIMATE=False If set to true, the result plays automatically
    # @param FILE=None If a file name is given, an mp4 is created of the animation
    # @param xBar=None If an xBar is given, it will be used as the measurement set
    # @param self.params.v0 takes a wind speed in meters per second
    # @param self.params.d0 takes a wind direction in degrees
    # @param self.params.vBar takes a wind speed in meters per second
    # @param self.params.dBar takes a wind direction in degrees
    # @param self.params.epSpeed is an epsilon over which to calculate df/dv (derivative of speed)
    # @param self.params.epDir is an epsilon over which to calculate df/dd (derivative of direction)
    # @param self.params.spErrMax is a speed error threshold for stopping iterations
    # @param self.params.dirErrMax is a direction error threshold for stopping iterations
    # @param self.params.iterMax is the maximum number of iterations to complete
    # @param self.params.mask_thresh threshold value for reducing the normalized sensitivity matrix
    #
    #
    # Generates a figure with four plots:
    #   1.  top left:       speed estimate error vs. iterations
    #   2.  bottom left:    direction estimate error vs. iterations
    #   3.  top right:      u_field estimation error. slider on bottom
    #                       allows you to see the map at different iterations
    #   4.  bottom right:   the mask being applied to the sensitivity matrix
    #                       for these calculations. the slider on the bottom
    #                       allows you to see the mask at different iterations
    #
    def reducedSM(self, MAT=False, ANIMATE=False, FILE=None, xBar=None, SHOW=True):
        # calculate f(vbar,dbar) = xBar
        print('vBar = ' + str(self.params.vBar))
        print('dBar = ' + str(self.params.dBar))
        vdBar = [[self.params.vBar], [self.params.dBar]]
        X = self.WF.x_mesh
        Y = self.WF.y_mesh

        TV = False # assume xBar is not time-varying
        if(xBar is None):
            print("no xBar provided, using vBar and dBar with FLORIS")
            self.WF.set_incoming(self.params.vBar,np.rad2deg(self.params.dBar))
            ff_bar = deepcopy(self.WF.u_field)
            ff = deepcopy(ff_bar)
            self.xBar = np.vstack(ff.flatten())  # 1xn
        elif xBar.ndim==2:
            print("Static xBar provided")
            ff_bar = deepcopy(xBar)
            ff = deepcopy(ff_bar)
            self.xBar = np.vstack(ff.flatten())
        else:
            print("Time varying xBar provided")
            TV = True
            self.xBar = deepcopy(xBar)

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
        ref_err_field_norm = list()
        err_field_norm = list()
        masks = list()
        mask_per = list()
        grid_sz=self.grid_resolution[0]*self.grid_resolution[1]
        errMax = 0
        self.WF.set_incoming(vdBar[0][0], np.rad2deg(vdBar[1][0]))
        ref_x_k = deepcopy(self.WF.u_field)
        while iters < self.params.iterMax and \
                ( abs(speedError[iters]) > self.params.spErrMax or \
                  abs(directionError[iters]) > self.params.dirErrMax ):
            # set new wind speed estimate for floris model
            self.WF.set_incoming(vd_k[0][0],np.rad2deg(vd_k[1][0]))
            # calculate f(vk,dk) = temp_x_k
            temp_x_k = deepcopy(self.WF.u_field)
            #print("temp_xk.shape: ",temp_x_k.shape)
            if TV: #if it is a time-varying SOWFA file, get the current ff_bar
                ff_bar = deepcopy(self.xBar[iters])
            err_field.append(abs(ff_bar-temp_x_k))
            err_field_norm.append(np.linalg.norm(abs(ff_bar-temp_x_k)))
            ref_err_field_norm.append(np.linalg.norm(abs(ff_bar-ref_x_k)))
            if np.amax(err_field[iters]) > errMax:
                errMax = np.amax(err_field[iters])
            # print vk,dk as new speed and direction estimates
            print('vk = ' + str(vd_k[0][0]))
            print('dk = ' + str(vd_k[1][0]))

            # calculate f(vk+ev,dk) = temp_x_ev
            self.WF.set_incoming(vd_k[0][0] + self.params.epSpeed, np.rad2deg(vd_k[1][0]))  # add speed epsilon to current estimate
            temp_x_ev = deepcopy(self.WF.u_field)

            # calculate f(vk,dk+ed) = temp_x_ed
            self.WF.set_incoming(vd_k[0][0],np.rad2deg(vd_k[1][0] + self.params.epDir))  # set wind direction to current estimate + epsilon
            temp_x_ed = deepcopy(self.WF.u_field)

            # calculate gradient
            temp_df_dv = (temp_x_ev - temp_x_k) / self.params.epSpeed  # partial of f WRT speed
            temp_df_dd = (temp_x_ed - temp_x_k) / self.params.epDir  # partial of f WRT direction
            temp_Xk = temp_x_k.flatten()

            # calculate mask
            temp_dMax = np.amax(np.abs(temp_df_dd))
            temp_vMax = np.amax(np.abs(temp_df_dv))
            temp_Zd = np.abs(temp_df_dd / temp_dMax)
            temp_Zv = np.abs(temp_df_dv / temp_vMax)
            temp_Z_mask = temp_Zd * temp_Zv / np.amax(temp_Zd * temp_Zv)
            temp_Z_mask = np.where(temp_Z_mask >= self.params.mask_thresh, 1, 0).flatten()
            sens_mat0 = list()
            sens_mat1 = list()
            Xbar = list()
            Xk = list()
            # removal of lines and columns
            for i in range(len(temp_Z_mask)):
                if temp_Z_mask[i] == 1:
                    if TV:
                        temp_xBar = np.vstack(self.xBar[iters].flatten())
                        Xbar.append(temp_xBar[i])
                        del temp_xBar
                    else:
                        Xbar.append(self.xBar[i])

                    Xk.append(temp_Xk[i])
                    sens_mat0.append(temp_df_dv.flatten()[i])
                    sens_mat1.append(temp_df_dd.flatten()[i])
            Xbar = np.vstack(np.array(Xbar))
            Xk = np.vstack(np.array(Xk))
            sens_mat0 = np.vstack(np.array(sens_mat0))
            sens_mat1 = np.vstack(np.array(sens_mat1))
            sens_mat_pinv = np.linalg.pinv(np.column_stack([sens_mat0, sens_mat1]))
            masks.append(temp_Z_mask)
            mask_per.append(100*np.sum(temp_Z_mask)/grid_sz)
            # calculate pseudoinverse[gradient{f(vk,dk)}]*{xBar-f(vk,dk)} = adj_vd (adjustment to current v&d estimates)
            adj_vd = np.matmul(sens_mat_pinv, Xbar - Xk)
            # calculate v_k+1 and d_k+1
            vd_kp1 = vd_k + adj_vd


            # update vd_k for next iteration
            vd_k = deepcopy(vd_kp1)
            iters = iters + 1
            iterations.append(iters)
            print('\n\niteration ' + str(iters) + ' complete.')
            # calculate error = [[vbar],[dbar]]-[[vk],[dk]]
            vdErr = vdBar - vd_k
            print('Speed error: ' + str(vdErr[0][0]))
            print('Direction error: ' + str(np.rad2deg(vdErr[1][0])) + '\N{DEGREE SIGN}')
            speedError.append(vdErr[0][0])
            directionError.append(vdErr[1][0])
            V_k.append(vd_k[0][0])
            D_k.append(vd_k[1][0])


            # delete temporary objects
            del temp_x_k
            del temp_x_ev
            del temp_x_ed
            del temp_df_dv
            del temp_df_dd
            del adj_vd
            del vd_kp1

        self.WF.set_incoming(vd_k[0][0],np.rad2deg(vd_k[1][0]))
        temp_x_k = deepcopy(self.WF.u_field)
        err_field.append(abs(ff_bar - temp_x_k))
        ref_err_field_norm.append(np.linalg.norm(ff_bar-ref_x_k))
        err_field_norm.append(np.linalg.norm(ff_bar-temp_x_k))
        if np.amax(err_field[iters]) > errMax:
            errMax = np.amax(err_field[iters])
        masks.append(temp_Z_mask)
        print('Max error = '+str(errMax))
        print('Total iterations: ' + str(iters))

        #### PLOTTING #####
        fontsize=14
        f = plt.figure(figsize=(10, 7.5))
        gs = f.add_gridspec(2, 3)
        sld_ax = plt.axes([0.23, 0.02, 0.56, 0.02])
        sld = Slider(sld_ax,
                     'iterations',
                     0,iters,valinit=0)
        sld.valtext.set_text('iteration 0')
        f.suptitle('Speed estimate: ' + str(self.params.v0) + ' m/s, Actual: ' + str(self.params.vBar) +
                   ' m/s\nDirection estimate: ' + str(np.rad2deg(self.params.d0)) + '\N{DEGREE SIGN}' +
                   ', Actual: ' + str(np.rad2deg(self.params.dBar)) + '\N{DEGREE SIGN}' +
                   '\nThreshold: ' + str(self.params.mask_thresh) + ' Iterations: ' + str(iters) +
                   '\nFinal error: e\N{GREEK SMALL LETTER THETA} = ' + str(np.rad2deg(vdErr[1][0])) +
                   '\N{DEGREE SIGN} ev = ' + str(vdErr[0][0]))
        spdErrAx = f.add_subplot(gs[0,0], title='speed error')
        spdErrAx.plot(iterations, speedError)
        dirErrAx = f.add_subplot(gs[1,0], title='direction error')
        dirErrAx.plot(iterations, np.rad2deg(directionError))
        dirErrAx.set_xlabel('iterations')
        errMapAx = f.add_subplot(gs[0,1:], title='u_field error')
        v = np.linspace(0, errMax, 100)
        V = np.linspace(0, errMax, 5)
        cont = errMapAx.contourf(X, Y, err_field[0],v,cmap='gnuplot2')

        cb = plt.colorbar(cont, ax=errMapAx)
        cb.set_clim(vmin=0, vmax=errMax)
        cb.set_ticks(V, True)
        cb.set_label('u_field error')
        cb.draw_all()

        maskAx = f.add_subplot(gs[1,1:])
        for tick in maskAx.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        for tick in maskAx.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        x, y = deepcopy(X), deepcopy(Y)
        scat = maskAx.scatter(x.flatten(),y.flatten(),masks[0],color='black')
        maskAx.text(-220, 290, str(self.params.v0) + ' m/s, ' + \
                      str(np.rad2deg(self.params.d0)) + '\N{DEGREE SIGN}', fontsize=fontsize)
        maskAx.arrow(-220, 225, 15 * self.params.v0 * np.cos(self.params.d0),
                       -15 * self.params.v0 * np.sin(self.params.d0), head_width=50, head_length=50)
        maskAx.set_xlabel('meters', fontsize=fontsize)
        ratio_default = (maskAx.get_xlim()[1] - maskAx.get_xlim()[0]) / \
                        (maskAx.get_ylim()[1] - maskAx.get_ylim()[0])
        print(ratio_default)
        errMapAx.set_aspect('equal')
        maskAx.set_aspect('equal')
        cb2 = plt.colorbar(mappable=cont, ax=maskAx)
        cb2.set_clim(vmin=0, vmax=errMax)
        cb2.set_ticks(V)
        cb2.draw_all()
        cb2.remove()
        i=0
        for coord, turbine in self.WF.turbines:
            x_0 = coord.x1 + np.sin(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            x_1 = coord.x1 - np.sin(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            y_0 = coord.x2 - np.cos(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            y_1 = coord.x2 + np.cos(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            errMapAx.plot([x_0, x_1], [y_0, y_1], color='white', linewidth=1)
            maskAx.plot([x_0, x_1], [y_0, y_1], color='black', linewidth=1)
            for tick in maskAx.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            for tick in maskAx.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            i += 1
        if(MAT):
            scio.savemat('mat/' + str(i)+'turbs_'+str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '-' + str(np.rad2deg(self.params.d0)) +\
                         '_'+str(self.params.mask_thresh)+'_spdERR.mat',
                         mdict={'spdErrTh'+str(int(self.params.mask_thresh*10)): speedError})
            scio.savemat('mat/' + str(i) + 'turbs_' + str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '-' + str(np.rad2deg(self.params.d0)) +\
                         '_'+str(self.params.mask_thresh)+'_dirERR.mat',
                         mdict={'dirErrTh'+str(int(self.params.mask_thresh*10)): directionError})
            scio.savemat('mat/'+str(i)+'turbs_'+str(self.params.vBar)+'-'+str(self.params.v0)+'_'+
                         str(np.rad2deg(self.params.dBar))+'_refERRnorm.mat',
                         mdict={'refERRnorm'+str(int(self.params.mask_thresh*10)): ref_err_field_norm})
            scio.savemat('mat/' + str(i) + 'turbs_' + str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '_ERRnorm.mat',
                         mdict={'ERRnorm' + str(int(self.params.mask_thresh * 10)): err_field_norm})
        print('max coverage: '+str(np.amax(mask_per))+'%')
        print('min coverage: '+str(np.amin(mask_per))+'%')
        print('average coverage: ' + str(np.average(mask_per)) + '%')
        plt.subplots_adjust(left=0.05,
                            bottom=0.15,
                            right=0.95,
                            top=0.83,
                            wspace=0.27,
                            hspace=0.19)
        def update_plot(val):

            idx = int(round(sld.val))
            sld.valtext.set_text('iteration '+'{}'.format(idx))
            spdErrAx.clear()
            dirErrAx.clear()
            errMapAx.clear()
            maskAx.clear()
            errMapAx.contourf(X, Y, err_field[idx], v, cmap='gnuplot2')
            maskAx.scatter(x,y,masks[idx],color='black')
            maskAx.text(-220, 290, str(round(V_k[idx],2)) + ' m/s, ' + \
                        str(round(np.rad2deg(D_k[idx]),2)) + '\N{DEGREE SIGN}', fontsize=fontsize)
            maskAx.arrow(-220, 225, 15 * V_k[idx] * np.cos(D_k[idx]),
                         -15 * V_k[idx] * np.sin(D_k[idx]), head_width=50, head_length=50)
            maskAx.set_xlabel('meters', fontsize=fontsize)
            for tick in maskAx.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            for tick in maskAx.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            dirErrAx.plot(iterations, np.rad2deg(directionError))
            dirErrAx.axvline(idx, color='red')
            dirErrAx.set_title('direction error')
            spdErrAx.plot(iterations, speedError)
            spdErrAx.axvline(idx, color='red')
            spdErrAx.set_title('speed error')

            for coord, turbine in self.WF.turbines:
                x_0 = coord.x1 + np.sin(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
                x_1 = coord.x1 - np.sin(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
                y_0 = coord.x2 - np.cos(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
                y_1 = coord.x2 + np.cos(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
                errMapAx.plot([x_0, x_1], [y_0, y_1], color='white', linewidth=1)
                maskAx.plot([x_0, x_1], [y_0, y_1], color='black', linewidth=1)

            plt.draw()
            print(sum(masks[idx]))
        sld.on_changed(update_plot)
        def animate(frame, *fargs):
            # if it's not at the end, increment the slider value
            if sld.val < sld.valmax-1:
                temp = sld.val
                sld.set_val(temp + 1)
            else:
            # if it's at the end, set it to the beginning
                sld.set_val(sld.valmin)

        # set the animate function to the FuncAnimation function for animation
        if ANIMATE:
            an = anim.FuncAnimation(f, animate, interval=100,frames=sld.valmax*5)
            # render to video. to make it play faster, increase fps
            if FILE:
                an.save(FILE+'.mp4',fps=7,dpi=300)
        if SHOW:
            plt.show()

    # reducedSM_vdy
    ## @brief Plots the convergence of speed, direction and yaw estimates
    # based on a reduced sensitivity matrix
    #
    # @param ANIMATE=False If set to true, the result plays automatically
    # @param FILE=None If a file name is given, an mp4 is created of the animation
    # @param xBar=None If an xBar is given, it will be used as the measurement set
    # @param self.params.v0 takes a wind speed in meters per second
    # @param self.params.d0 takes a wind direction in radians
    # @param self.params.y0 takes a turbine yaw in radians
    # @param self.params.vBar takes a wind speed in meters per second
    # @param self.params.dBar takes a wind direction in radians
    # @param self.params.yBar takes a turbine yaw in radians
    # @param self.params.epSpeed is an epsilon over which to calculate df/dv (derivative of speed)
    # @param self.params.epDir is an epsilon over which to calculate df/dd (derivative of direction)
    # @param self.params.epYaw is an epsilon over which to calculate df/dy (derivative of yaw)
    # @param self.params.spErrMax is a speed error threshold for stopping iterations
    # @param self.params.dirErrMax is a direction error threshold for stopping iterations
    # @param self.params.yawErrMax is a yaw error threshold for stoppin iterations
    # @param self.params.iterMax is the maximum number of iterations to complete
    # @param self.params.mask_thresh threshold value for reducing the normalized sensitivity matrix
    #
    #
    # Generates a figure with five plots:
    #   1.  top left:       speed estimate error vs. iterations
    #   2.  center left:    direction estimate error vs. iterations
    #   3.  bottom left:    yaw estimate error vs. iterations
    #   4.  top right:      u_field estimation error. slider on bottom
    #                       allows you to see the map at different iterations
    #   5.  bottom right:   the mask being applied to the sensitivity matrix
    #                       for these calculations. the slider on the bottom
    #                       allows you to see the mask at different iterations
    #
    def reducedSM_vdy(self, MAT=False, ANIMATE=False, FILE=None, xBar=None, SHOW=True):
        # calculate f(vbar,dbar) = xBar
        print('vBar = ' + str(self.params.vBar)+' m/s')
        print('dBar = ' + str(np.deg2rad(self.params.dBar))+'\N{DEGREE SIGN}')
        print('yBar = ' + str(np.rad2deg(self.params.yBar))+'\N{DEGREE SIGN}')
        vdyBar = [[self.params.vBar], [self.params.dBar],[self.params.yBar]]
        X = self.WF.x_mesh
        Y = self.WF.y_mesh

        TV = False  # assume xBar is not time-varying
        if (xBar is None):
            print("no xBar provided, using vBar and dBar with FLORIS")
            self.WF.set_vdy(self.params.vBar, np.rad2deg(self.params.dBar),np.rad2deg(self.params.yBar))
            ff_bar = deepcopy(self.WF.u_field)
            ff = deepcopy(ff_bar)
            self.xBar = np.vstack(ff.flatten())  # 1xn
        elif xBar.ndim == 2:
            print("Static xBar provided")
            ff_bar = deepcopy(xBar)
            ff = deepcopy(ff_bar)
            self.xBar = np.vstack(ff.flatten())
        else:
            print("Time varying xBar provided")
            TV = True
            self.xBar = deepcopy(xBar)

        # initialize error
        vdyErr = [[self.params.vBar - self.params.v0],
                 [self.params.dBar - self.params.d0],
                 [self.params.yBar - self.params.y0]]  # initial error
        speedError = list()  # list for history of speed error
        directionError = list()  # list for history of direction error
        yawError = list()
        print('Initial speed error = ' + str(vdyErr[0][0])+' m/s')
        print('Initial direction error = ' + str(np.rad2deg(vdyErr[1][0]))+'\N{DEGREE SIGN}')
        print('Initial yaw error = ' + str(np.rad2deg(vdyErr[2][0]))+'\N{DEGREE SIGN}')
        speedError.append(vdyErr[0][0])
        directionError.append(vdyErr[1][0])
        yawError.append(vdyErr[2][0])
        # set up first iteration
        iters = 0
        iterations = list()  # iteration list for graphing
        iterations.append(iters)
        V_k = list()  # list for history of vk
        D_k = list()  # list for history of dk
        Y_k = list()  # list for history of yk
        V_k.append(self.params.v0)
        D_k.append(self.params.d0)
        Y_k.append(self.params.y0)
        vdy_k = [[self.params.v0], [self.params.d0], [self.params.y0]]  # initial estimate
        err_field = list()
        ref_err_field_norm = list()
        err_field_norm = list()
        masks = list()
        mask_per = list()
        grid_sz = self.grid_resolution[0] * self.grid_resolution[1]
        errMax = 0
        self.WF.set_vdy(vdyBar[0][0], np.rad2deg(vdyBar[1][0]), np.rad2deg(vdyBar[2][0]))
        ref_x_k = deepcopy(self.WF.u_field)
        while iters < self.params.iterMax and \
                (abs(speedError[iters]) > self.params.spErrMax or \
                 abs(directionError[iters]) > self.params.dirErrMax or \
                        abs(yawError[iters]) > self.params.yawErrMax):
            # set new wind speed estimate for floris model
            self.WF.set_vdy(vdy_k[0][0], np.rad2deg(vdy_k[1][0]), np.rad2deg(vdy_k[2][0]))
            # calculate f(vk,dk) = temp_x_k
            temp_x_k = deepcopy(self.WF.u_field)
            # print("temp_xk.shape: ",temp_x_k.shape)
            if TV:  # if it is a time-varying SOWFA file, get the current ff_bar
                ff_bar = deepcopy(self.xBar[iters])
            err_field.append(abs(ff_bar - temp_x_k))
            err_field_norm.append(np.linalg.norm(abs(ff_bar - temp_x_k)))
            ref_err_field_norm.append(np.linalg.norm(abs(ff_bar - ref_x_k)))
            if np.amax(err_field[iters]) > errMax:
                errMax = np.amax(err_field[iters])
            # print vk,dk as new speed and direction estimates
            print('vk = ' + str(vdy_k[0][0])+' m/s')
            print('dk = ' + str(np.rad2deg(vdy_k[1][0]))+'\N{DEGREE SIGN}')
            print('yk = ' + str(np.rad2deg(vdy_k[2][0]))+'\N{DEGREE SIGN}')
            # calculate f(vk+ev,dk,yk) = temp_x_ev
            self.WF.set_vdy(vdy_k[0][0] + self.params.epSpeed,
                                 np.rad2deg(vdy_k[1][0]),
                                 np.rad2deg(vdy_k[2][0]))  # add speed epsilon to current estimate
            temp_x_ev = deepcopy(self.WF.u_field)

            # calculate f(vk,dk+ed,yk) = temp_x_ed
            self.WF.set_vdy(vdy_k[0][0],
                           np.rad2deg(vdy_k[1][0] + self.params.epDir),
                           np.rad2deg(vdy_k[2][0]))  # set wind direction to current estimate + epsilon
            temp_x_ed = deepcopy(self.WF.u_field)

            # calculate f(vk,dk,yk+ey) = temp_x_ed
            self.WF.set_vdy(vdy_k[0][0],
                           np.rad2deg(vdy_k[1][0]),
                           np.rad2deg(vdy_k[2][0] + self.params.epYaw))  # set wind direction to current estimate + epsilon
            temp_x_ey = deepcopy(self.WF.u_field)

            # calculate gradient
            temp_df_dv = (temp_x_ev - temp_x_k) / self.params.epSpeed  # partial of f WRT speed
            temp_df_dd = (temp_x_ed - temp_x_k) / self.params.epDir  # partial of f WRT direction
            temp_df_dy = (temp_x_ey - temp_x_k) / self.params.epYaw  # partial of f WRT yaw
            temp_Xk = temp_x_k.flatten()

            # calculate mask
            temp_dMax = np.amax(np.abs(temp_df_dd))
            temp_vMax = np.amax(np.abs(temp_df_dv))
            temp_yMax = np.amax(np.abs(temp_df_dy))
            temp_Zd = np.abs(temp_df_dd / temp_dMax)
            temp_Zv = np.abs(temp_df_dv / temp_vMax)
            temp_Zy = np.abs(temp_df_dy / temp_yMax)
            temp_Z_mask = temp_Zd * temp_Zv * temp_Zy/ np.amax(temp_Zd * temp_Zv * temp_Zy)
            temp_Z_mask = np.where(temp_Z_mask >= self.params.mask_thresh, 1, 0).flatten()
            sens_mat0 = list()
            sens_mat1 = list()
            sens_mat2 = list()
            Xbar = list()
            Xk = list()
            # removal of lines and columns
            for i in range(len(temp_Z_mask)):
                if temp_Z_mask[i] == 1:
                    if TV:
                        temp_xBar = np.vstack(self.xBar[iters].flatten())
                        Xbar.append(temp_xBar[i])
                        del temp_xBar
                    else:
                        Xbar.append(self.xBar[i])

                    Xk.append(temp_Xk[i])
                    sens_mat0.append(temp_df_dv.flatten()[i])
                    sens_mat1.append(temp_df_dd.flatten()[i])
                    sens_mat2.append(temp_df_dy.flatten()[i])
            Xbar = np.vstack(np.array(Xbar))
            Xk = np.vstack(np.array(Xk))
            sens_mat0 = np.vstack(np.array(sens_mat0))
            sens_mat1 = np.vstack(np.array(sens_mat1))
            sens_mat2 = np.vstack(np.array(sens_mat2))
            sens_mat_pinv = np.linalg.pinv(np.column_stack([sens_mat0, sens_mat1, sens_mat2]))
            masks.append(temp_Z_mask)
            mask_per.append(100 * np.sum(temp_Z_mask) / grid_sz)
            # calculate pseudoinverse[gradient{f(vk,dk)}]*{xBar-f(vk,dk)} = adj_vd (adjustment to current v&d estimates)
            adj_vdy = np.matmul(sens_mat_pinv, Xbar - Xk)
            # calculate v_k+1 and d_k+1
            vdy_kp1 = vdy_k + adj_vdy

            # update vd_k for next iteration
            vdy_k = deepcopy(vdy_kp1)
            iters = iters + 1
            iterations.append(iters)
            print('\n\niteration ' + str(iters) + ' complete.')
            # calculate error = [[vbar],[dbar]]-[[vk],[dk]]
            vdyErr = vdyBar - vdy_k
            print('Speed error: ' + str(vdyErr[0][0])+' m/s')
            print('Direction error: ' + str(np.rad2deg(vdyErr[1][0])) + '\N{DEGREE SIGN}')
            print('Yaw error: ' + str(np.rad2deg(vdyErr[2][0])) + '\N{DEGREE SIGN}')
            speedError.append(vdyErr[0][0])
            directionError.append(vdyErr[1][0])
            yawError.append(vdyErr[2][0])
            V_k.append(vdy_k[0][0])
            D_k.append(vdy_k[1][0])
            Y_k.append(vdy_k[2][0])

            # delete temporary objects
            del temp_x_k
            del temp_x_ev
            del temp_x_ed
            del temp_x_ey
            del temp_df_dv
            del temp_df_dd
            del temp_df_dy
            del adj_vdy
            del vdy_kp1

        self.WF.set_vdy(vdy_k[0][0], np.rad2deg(vdy_k[1][0]), np.rad2deg(vdy_k[2][0]))
        temp_x_k = deepcopy(self.WF.u_field)
        err_field.append(abs(ff_bar - temp_x_k))
        ref_err_field_norm.append(np.linalg.norm(ff_bar - ref_x_k))
        err_field_norm.append(np.linalg.norm(ff_bar - temp_x_k))
        if np.amax(err_field[iters]) > errMax:
            errMax = np.amax(err_field[iters])
        masks.append(temp_Z_mask)
        print('Max error = ' + str(errMax))
        print('Total iterations: ' + str(iters))

        #### PLOTTING #####
        fontsize = 14
        f = plt.figure(figsize=(10, 7.5))
        gs = f.add_gridspec(6, 3)
        sld_ax = plt.axes([0.23, 0.02, 0.56, 0.02])
        sld = Slider(sld_ax,
                     'iterations',
                     0, iters, valinit=0)
        sld.valtext.set_text('iteration 0')
        f.suptitle('Speed estimate: ' + str(self.params.v0) + ' m/s, Actual: ' + str(self.params.vBar) +
                   ' m/s\nDirection estimate: ' + str(np.rad2deg(self.params.d0)) + '\N{DEGREE SIGN}' +
                   ', Actual: ' + str(np.rad2deg(self.params.dBar)) + '\N{DEGREE SIGN}' +
                   '\nYaw estimate: ' + str(self.params.y0) + '\N{DEGREE SIGN}' +
                   ', Actual: ' + str(np.rad2deg(self.params.yBar)) + '\N{DEGREE SIGN}' +
                   '\nThreshold: ' + str(self.params.mask_thresh) + ' Iterations: ' + str(iters) +
                   '\nFinal error: e\N{GREEK SMALL LETTER THETA} = ' + str(np.rad2deg(vdyErr[1][0])) +
                   '\N{DEGREE SIGN}, e\N{GREEK SMALL LETTER PHI} = ' + str(np.rad2deg(vdyErr[2][0])) +
                   '\N{DEGREE SIGN}, ev = ' + str(vdyErr[0][0]) + 'm/s, ')
        spdErrAx = f.add_subplot(gs[0:2, 0], title='speed error')
        spdErrAx.plot(iterations, speedError)
        dirErrAx = f.add_subplot(gs[2:4, 0], title='direction error')
        dirErrAx.plot(iterations, np.rad2deg(directionError))
        yawErrAx = f.add_subplot(gs[4:, 0], title='yaw error')
        yawErrAx.plot(iterations, np.rad2deg(yawError))
        yawErrAx.set_xlabel('iterations')
        errMapAx = f.add_subplot(gs[0:3, 1:], title='u_field error')
        v = np.linspace(0, errMax, 100)
        V = np.linspace(0, errMax, 5)
        cont = errMapAx.contourf(X, Y, err_field[0], v, cmap='gnuplot2')

        cb = plt.colorbar(cont, ax=errMapAx)
        cb.set_clim(vmin=0, vmax=errMax)
        cb.set_ticks(V, True)
        cb.set_label('u_field error')
        cb.draw_all()

        maskAx = f.add_subplot(gs[3:, 1:])
        for tick in maskAx.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        for tick in maskAx.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        x, y = deepcopy(X), deepcopy(Y)
        scat = maskAx.scatter(x.flatten(), y.flatten(), masks[0], color='black')
        maskAx.text(-220, 290, str(self.params.v0) + ' m/s, ' + \
                    str(np.rad2deg(self.params.d0)) + '\N{DEGREE SIGN}, ' + \
                    str(np.rad2deg(self.params.y0)) + '\N{DEGREE SIGN}', fontsize=fontsize)
        maskAx.arrow(-220, 225, 15 * self.params.v0 * np.cos(self.params.d0),
                     -15 * self.params.v0 * np.sin(self.params.d0), head_width=50, head_length=50)
        maskAx.set_xlabel('meters', fontsize=fontsize)
        ratio_default = (maskAx.get_xlim()[1] - maskAx.get_xlim()[0]) / \
                        (maskAx.get_ylim()[1] - maskAx.get_ylim()[0])
        print(ratio_default)
        errMapAx.set_aspect('equal')
        maskAx.set_aspect('equal')
        cb2 = plt.colorbar(mappable=cont, ax=maskAx)
        cb2.set_clim(vmin=0, vmax=errMax)
        cb2.set_ticks(V)
        cb2.draw_all()
        cb2.remove()
        i = 0
        for coord, turbine in self.WF.turbines:
            x_0 = coord.x1 + np.sin(self.params.y0) * turbine.rotor_radius
            x_1 = coord.x1 - np.sin(self.params.y0) * turbine.rotor_radius
            y_0 = coord.x2 - np.cos(self.params.y0) * turbine.rotor_radius
            y_1 = coord.x2 + np.cos(self.params.y0) * turbine.rotor_radius
            X_0 = coord.x1 + np.sin(self.params.yBar) * turbine.rotor_radius
            X_1 = coord.x1 - np.sin(self.params.yBar) * turbine.rotor_radius
            Y_0 = coord.x2 - np.cos(self.params.yBar) * turbine.rotor_radius
            Y_1 = coord.x2 + np.cos(self.params.yBar) * turbine.rotor_radius

            errMapAx.plot([x_0, x_1], [y_0, y_1], color='lime', linewidth=1)
            maskAx.plot([x_0, x_1], [y_0, y_1], color='red', linewidth=1)
            errMapAx.plot([X_0, X_1], [Y_0, Y_1], color='white', linewidth=1)
            maskAx.plot([X_0, X_1], [Y_0, Y_1], color='black', linewidth=1)

            for tick in maskAx.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            for tick in maskAx.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            i += 1
        if (MAT):
            scio.savemat('mat/' + str(i) + 'turbs_' + str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '-' + str(np.rad2deg(self.params.d0)) + \
                         str(np.rad2deg(self.params.yBar)) + '-' + str(np.rad2deg(self.params.y0)) + \
                         '_' + str(self.params.mask_thresh) + '_spdERR.mat',
                         mdict={'spdErrTh' + str(int(self.params.mask_thresh * 10)): speedError})
            scio.savemat('mat/' + str(i) + 'turbs_' + str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '-' + str(np.rad2deg(self.params.d0)) + \
                         str(np.rad2deg(self.params.yBar)) + '-' + str(np.rad2deg(self.params.y0)) + \
                         '_' + str(self.params.mask_thresh) + '_dirERR.mat',
                         mdict={'dirErrTh' + str(int(self.params.mask_thresh * 10)): directionError})
            scio.savemat('mat/' + str(i) + 'turbs_' + str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '-' + str(np.rad2deg(self.params.d0)) + \
                         str(np.rad2deg(self.params.yBar)) + '-' + str(np.rad2deg(self.params.y0)) + \
                         '_' + str(self.params.mask_thresh) + '_yawERR.mat',
                         mdict={'yawErrTh' + str(int(self.params.mask_thresh * 10)): yawError})

            scio.savemat('mat/' + str(i) + 'turbs_' + str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '_refERRnorm.mat',
                         mdict={'refERRnorm' + str(int(self.params.mask_thresh * 10)): ref_err_field_norm})
            scio.savemat('mat/' + str(i) + 'turbs_' + str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '_ERRnorm.mat',
                         mdict={'ERRnorm' + str(int(self.params.mask_thresh * 10)): err_field_norm})
        print('max coverage: ' + str(np.amax(mask_per)) + '%')
        print('min coverage: ' + str(np.amin(mask_per)) + '%')
        print('average coverage: ' + str(np.average(mask_per)) + '%')
        plt.subplots_adjust(left=0.05,
                            bottom=0.15,
                            right=0.95,
                            top=0.83,
                            wspace=0.27,
                            hspace=0.19)

        def update_plot(val):

            idx = int(round(sld.val))
            sld.valtext.set_text('iteration ' + '{}'.format(idx))
            spdErrAx.clear()
            dirErrAx.clear()
            yawErrAx.clear()
            errMapAx.clear()
            maskAx.clear()
            errMapAx.contourf(X, Y, err_field[idx], v, cmap='gnuplot2')
            maskAx.scatter(x, y, masks[idx], color='black')
            maskAx.text(-220, 290, str(round(V_k[idx], 2)) + ' m/s, ' + \
                        str(round(np.rad2deg(D_k[idx]), 2)) + '\N{DEGREE SIGN}, ' + \
                        str(round(np.rad2deg(Y_k[idx]), 2)) + '\N{DEGREE SIGN}', fontsize=fontsize)
            maskAx.arrow(-220, 225, 15 * V_k[idx] * np.cos(D_k[idx]),
                         -15 * V_k[idx] * np.sin(D_k[idx]), head_width=50, head_length=50)
            maskAx.set_xlabel('meters', fontsize=fontsize)
            for tick in maskAx.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            for tick in maskAx.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            dirErrAx.plot(iterations, np.rad2deg(directionError))
            dirErrAx.axvline(idx, color='red')
            dirErrAx.set_title('direction error')
            spdErrAx.plot(iterations, speedError)
            spdErrAx.axvline(idx, color='red')
            spdErrAx.set_title('speed error')
            yawErrAx.plot(iterations, yawError)
            yawErrAx.axvline(idx, color='red')
            yawErrAx.set_title('yaw error')

            for coord, turbine in self.WF.turbines:
                x_0 = coord.x1 + np.sin(Y_k[idx]) * turbine.rotor_radius
                x_1 = coord.x1 - np.sin(Y_k[idx]) * turbine.rotor_radius
                y_0 = coord.x2 - np.cos(Y_k[idx]) * turbine.rotor_radius
                y_1 = coord.x2 + np.cos(Y_k[idx]) * turbine.rotor_radius
                X_0 = coord.x1 + np.sin(self.params.yBar) * turbine.rotor_radius
                X_1 = coord.x1 - np.sin(self.params.yBar) * turbine.rotor_radius
                Y_0 = coord.x2 - np.cos(self.params.yBar) * turbine.rotor_radius
                Y_1 = coord.x2 + np.cos(self.params.yBar) * turbine.rotor_radius

                errMapAx.plot([x_0, x_1], [y_0, y_1], color='lime', linewidth=1)
                maskAx.plot([x_0, x_1], [y_0, y_1], color='red', linewidth=1)
                errMapAx.plot([X_0, X_1], [Y_0, Y_1], color='white', linewidth=1)
                maskAx.plot([X_0, X_1], [Y_0, Y_1], color='black', linewidth=1)
                plt.draw()
            print(sum(masks[idx]))

        sld.on_changed(update_plot)

        def animate(frame, *fargs):
            # if it's not at the end, increment the slider value
            if sld.val < sld.valmax - 1:
                temp = sld.val
                sld.set_val(temp + 1)
            else:
                # if it's at the end, set it to the beginning
                sld.set_val(sld.valmin)

        # set the animate function to the FuncAnimation function for animation
        if ANIMATE:
            an = anim.FuncAnimation(f, animate, interval=100, frames=sld.valmax * 5)
            # render to video. to make it play faster, increase fps
            if FILE:
                an.save(FILE + '.mp4', fps=7, dpi=300)
        if SHOW:
            plt.show()

    # reducedSM_vy
    ## @brief Plots the convergence of speed and yaw estimates
    # based on a reduced sensitivity matrix
    #
    # @param ANIMATE=False If set to true, the result plays automatically
    # @param FILE=None If a file name is given, an mp4 is created of the animation
    # @param xBar=None If an xBar is given, it will be used as the measurement set
    # @param self.params.v0 takes a wind speed in meters per second
    # @param self.params.y0 takes a turbine yaw in radians
    # @param self.params.vBar takes a wind speed in meters per second
    # @param self.params.yBar takes a turbine yaw in radians
    # @param self.params.epSpeed is an epsilon over which to calculate df/dv (derivative of speed)
    # @param self.params.epYaw is an epsilon over which to calculate df/dy (derivative of yaw)
    # @param self.params.spErrMax is a speed error threshold for stopping iterations
    # @param self.params.yawErrMax is a direction error threshold for stopping iterations
    # @param self.params.iterMax is the maximum number of iterations to complete
    # @param self.params.mask_thresh threshold value for reducing the normalized sensitivity matrix
    #
    #
    # Generates a figure with four plots:
    #   1.  top left:       speed estimate error vs. iterations
    #   2.  bottom left:    yaw estimate error vs. iterations
    #   3.  top right:      u_field estimation error. slider on bottom
    #                       allows you to see the map at different iterations
    #   4.  bottom right:   the mask being applied to the sensitivity matrix
    #                       for these calculations. the slider on the bottom
    #                       allows you to see the mask at different iterations
    #
    def reducedSM_vy(self, MAT=False, ANIMATE=False, FILE=None, xBar=None, SHOW=True):
        # calculate f(vbar,dbar) = xBar
        print('vBar = ' + str(self.params.vBar) + ' m/s')
        print('yBar = ' + str(np.rad2deg(self.params.yBar)) + '\N{DEGREE SIGN}')
        vyBar = [[self.params.vBar], [self.params.yBar]]
        X = self.WF.x_mesh
        Y = self.WF.y_mesh

        TV = False  # assume xBar is not time-varying
        if (xBar is None):
            print("no xBar provided, using vBar and dBar with FLORIS")
            self.WF.set_vy(self.params.vBar, np.rad2deg(self.params.yBar))
            ff_bar = deepcopy(self.WF.u_field)
            ff = deepcopy(ff_bar)
            self.xBar = np.vstack(ff.flatten())  # 1xn
        elif xBar.ndim == 2:
            print("Static xBar provided")
            ff_bar = deepcopy(xBar)
            ff = deepcopy(ff_bar)
            self.xBar = np.vstack(ff.flatten())
        else:
            print("Time varying xBar provided")
            TV = True
            self.xBar = deepcopy(xBar)

        # initialize error
        vyErr = [[self.params.vBar - self.params.v0],
                 [self.params.yBar - self.params.y0]]  # initial error
        speedError = list()  # list for history of speed error
        yawError = list()
        print('Initial speed error = ' + str(vyErr[0][0]) + ' m/s')
        print('Initial yaw error = ' + str(np.rad2deg(vyErr[1][0])) + '\N{DEGREE SIGN}')
        speedError.append(vyErr[0][0])
        yawError.append(vyErr[1][0])
        # set up first iteration
        iters = 0
        iterations = list()  # iteration list for graphing
        iterations.append(iters)
        V_k = list()  # list for history of vk
        Y_k = list()  # list for history of yk
        V_k.append(self.params.v0)
        Y_k.append(self.params.y0)
        vy_k = [[self.params.v0], [self.params.y0]]  # initial estimate
        err_field = list()
        ref_err_field_norm = list()
        err_field_norm = list()
        masks = list()
        mask_per = list()
        grid_sz = self.grid_resolution[0] * self.grid_resolution[1]
        errMax = 0
        self.WF.set_vy(vyBar[0][0], np.rad2deg(vyBar[1][0]))
        ref_x_k = deepcopy(self.WF.u_field)
        while iters < self.params.iterMax and \
                (abs(speedError[iters]) > self.params.spErrMax or \
                 abs(yawError[iters]) > self.params.yawErrMax):
            # set new wind speed estimate for floris model
            self.WF.set_vy(vy_k[0][0], np.rad2deg(vy_k[1][0]))
            # calculate f(vk,dk) = temp_x_k
            temp_x_k = deepcopy(self.WF.u_field)
            # print("temp_xk.shape: ",temp_x_k.shape)
            if TV:  # if it is a time-varying SOWFA file, get the current ff_bar
                ff_bar = deepcopy(self.xBar[iters])
            err_field.append(abs(ff_bar - temp_x_k))
            err_field_norm.append(np.linalg.norm(abs(ff_bar - temp_x_k)))
            ref_err_field_norm.append(np.linalg.norm(abs(ff_bar - ref_x_k)))
            if np.amax(err_field[iters]) > errMax:
                errMax = np.amax(err_field[iters])
            # print vk,dk as new speed and direction estimates
            print('vk = ' + str(vy_k[0][0]) + ' m/s')
            print('yk = ' + str(np.rad2deg(vy_k[1][0])) + '\N{DEGREE SIGN}')
            # calculate f(vk+ev,yk) = temp_x_ev
            self.WF.set_vy(vy_k[0][0] + self.params.epSpeed,
                        np.rad2deg(vy_k[1][0]))  # add speed epsilon to current estimate
            temp_x_ev = deepcopy(self.WF.u_field)

            # calculate f(vk,yk+ey) = temp_x_ed
            self.WF.set_vy(vy_k[0][0],
                        np.rad2deg(vy_k[1][0] + self.params.epYaw))  # set wind direction to current estimate + epsilon
            temp_x_ey = deepcopy(self.WF.u_field)

            # calculate gradient
            temp_df_dv = (temp_x_ev - temp_x_k) / self.params.epSpeed  # partial of f WRT speed
            temp_df_dy = (temp_x_ey - temp_x_k) / self.params.epYaw  # partial of f WRT yaw
            temp_Xk = temp_x_k.flatten()

            # calculate mask
            temp_vMax = np.amax(np.abs(temp_df_dv))
            temp_yMax = np.amax(np.abs(temp_df_dy))
            temp_Zv = np.abs(temp_df_dv / temp_vMax)
            temp_Zy = np.abs(temp_df_dy / temp_yMax)
            temp_Z_mask = temp_Zv * temp_Zy / np.amax(temp_Zv * temp_Zy)
            temp_Z_mask = np.where(temp_Z_mask >= self.params.mask_thresh, 1, 0).flatten()
            sens_mat0 = list()
            sens_mat1 = list()
            Xbar = list()
            Xk = list()
            # removal of lines and columns
            for i in range(len(temp_Z_mask)):
                if temp_Z_mask[i] == 1:
                    if TV:
                        temp_xBar = np.vstack(self.xBar[iters].flatten())
                        Xbar.append(temp_xBar[i])
                        del temp_xBar
                    else:
                        Xbar.append(self.xBar[i])

                    Xk.append(temp_Xk[i])
                    sens_mat0.append(temp_df_dv.flatten()[i])
                    sens_mat1.append(temp_df_dy.flatten()[i])
            Xbar = np.vstack(np.array(Xbar))
            Xk = np.vstack(np.array(Xk))
            sens_mat0 = np.vstack(np.array(sens_mat0))
            sens_mat1 = np.vstack(np.array(sens_mat1))
            sens_mat_pinv = np.linalg.pinv(np.column_stack([sens_mat0, sens_mat1]))
            masks.append(temp_Z_mask)
            mask_per.append(100 * np.sum(temp_Z_mask) / grid_sz)
            # calculate pseudoinverse[gradient{f(vk,dk)}]*{xBar-f(vk,dk)} = adj_vd (adjustment to current v&d estimates)
            adj_vy = np.matmul(sens_mat_pinv, Xbar - Xk)
            # calculate v_k+1 and d_k+1
            vy_kp1 = vy_k + adj_vy

            # update vd_k for next iteration
            vy_k = deepcopy(vy_kp1)
            iters = iters + 1
            iterations.append(iters)
            print('\n\niteration ' + str(iters) + ' complete.')
            # calculate error = [[vbar],[dbar]]-[[vk],[dk]]
            vyErr = vyBar - vy_k
            print('Speed error: ' + str(vyErr[0][0]) + ' m/s')
            print('Yaw error: ' + str(np.rad2deg(vyErr[1][0])) + '\N{DEGREE SIGN}')
            speedError.append(vyErr[0][0])
            yawError.append(vyErr[1][0])
            V_k.append(vy_k[0][0])
            Y_k.append(vy_k[1][0])

            # delete temporary objects
            del temp_x_k
            del temp_x_ev
            del temp_x_ey
            del temp_df_dv
            del temp_df_dy
            del adj_vy
            del vy_kp1

        self.WF.set_vy(vy_k[0][0], np.rad2deg(vy_k[1][0]))
        temp_x_k = deepcopy(self.WF.u_field)
        err_field.append(abs(ff_bar - temp_x_k))
        ref_err_field_norm.append(np.linalg.norm(ff_bar - ref_x_k))
        err_field_norm.append(np.linalg.norm(ff_bar - temp_x_k))
        if np.amax(err_field[iters]) > errMax:
            errMax = np.amax(err_field[iters])
        masks.append(temp_Z_mask)
        print('Max error = ' + str(errMax))
        print('Total iterations: ' + str(iters))

        #### PLOTTING #####
        fontsize = 14
        f = plt.figure(figsize=(10, 7.5))
        gs = f.add_gridspec(6, 3)
        sld_ax = plt.axes([0.23, 0.02, 0.56, 0.02])
        sld = Slider(sld_ax,
                     'iterations',
                     0, iters, valinit=0)
        sld.valtext.set_text('iteration 0')
        f.suptitle('Speed estimate: ' + str(self.params.v0) + ' m/s, Actual: ' + str(self.params.vBar) +
                   '\nYaw estimate: ' + str(self.params.y0) + '\N{DEGREE SIGN}' +
                   ', Actual: ' + str(np.rad2deg(self.params.yBar)) + '\N{DEGREE SIGN}' +
                   '\nThreshold: ' + str(self.params.mask_thresh) + ' Iterations: ' + str(iters) +
                   '\nFinal error: e\N{GREEK SMALL LETTER PHI} = ' + str(np.rad2deg(vyErr[1][0])) +
                   '\N{DEGREE SIGN}, ev = ' + str(vyErr[0][0]) + 'm/s, ')
        spdErrAx = f.add_subplot(gs[0:3, 0], title='speed error')
        spdErrAx.plot(iterations, speedError)
        yawErrAx = f.add_subplot(gs[3:, 0], title='yaw error')
        yawErrAx.plot(iterations, np.rad2deg(yawError))
        yawErrAx.set_xlabel('iterations')
        errMapAx = f.add_subplot(gs[0:3, 1:], title='u_field error')
        v = np.linspace(0, errMax, 100)
        V = np.linspace(0, errMax, 5)
        cont = errMapAx.contourf(X, Y, err_field[0], v, cmap='gnuplot2')

        cb = plt.colorbar(cont, ax=errMapAx)
        cb.set_clim(vmin=0, vmax=errMax)
        cb.set_ticks(V, True)
        cb.set_label('u_field error')
        cb.draw_all()

        maskAx = f.add_subplot(gs[3:, 1:])
        for tick in maskAx.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        for tick in maskAx.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        x, y = deepcopy(X), deepcopy(Y)
        scat = maskAx.scatter(x.flatten(), y.flatten(), masks[0], color='black')
        maskAx.text(-220, 290, str(self.params.v0) + ' m/s, ' + \
                    str(np.rad2deg(self.params.y0)) + '\N{DEGREE SIGN}', fontsize=fontsize)
        maskAx.arrow(-220, 225, 15 * self.params.v0 * np.cos(self.params.y0),
                     -15 * self.params.v0 * np.sin(self.params.y0), head_width=50, head_length=50)
        maskAx.set_xlabel('meters', fontsize=fontsize)
        ratio_default = (maskAx.get_xlim()[1] - maskAx.get_xlim()[0]) / \
                        (maskAx.get_ylim()[1] - maskAx.get_ylim()[0])
        print(ratio_default)
        errMapAx.set_aspect('equal')
        maskAx.set_aspect('equal')
        cb2 = plt.colorbar(mappable=cont, ax=maskAx)
        cb2.set_clim(vmin=0, vmax=errMax)
        cb2.set_ticks(V)
        cb2.draw_all()
        cb2.remove()
        i = 0
        for coord, turbine in self.WF.turbines:
            x_0 = coord.x1 + np.sin(self.params.y0) * turbine.rotor_radius
            x_1 = coord.x1 - np.sin(self.params.y0) * turbine.rotor_radius
            y_0 = coord.x2 - np.cos(self.params.y0) * turbine.rotor_radius
            y_1 = coord.x2 + np.cos(self.params.y0) * turbine.rotor_radius
            X_0 = coord.x1 + np.sin(self.params.yBar) * turbine.rotor_radius
            X_1 = coord.x1 - np.sin(self.params.yBar) * turbine.rotor_radius
            Y_0 = coord.x2 - np.cos(self.params.yBar) * turbine.rotor_radius
            Y_1 = coord.x2 + np.cos(self.params.yBar) * turbine.rotor_radius

            errMapAx.plot([x_0, x_1], [y_0, y_1], color='lime', linewidth=1)
            maskAx.plot([x_0, x_1], [y_0, y_1], color='red', linewidth=1)
            errMapAx.plot([X_0, X_1], [Y_0, Y_1], color='white', linewidth=1)
            maskAx.plot([X_0, X_1], [Y_0, Y_1], color='black', linewidth=1)

            for tick in maskAx.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            for tick in maskAx.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            i += 1
        if (MAT):
            scio.savemat('mat/' + str(i) + 'turbs_' + str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '-' + str(np.rad2deg(self.params.d0)) + \
                         str(np.rad2deg(self.params.yBar)) + '-' + str(np.rad2deg(self.params.y0)) + \
                         '_' + str(self.params.mask_thresh) + '_spdERR.mat',
                         mdict={'spdErrTh' + str(int(self.params.mask_thresh * 10)): speedError})
            scio.savemat('mat/' + str(i) + 'turbs_' + str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '-' + str(np.rad2deg(self.params.d0)) + \
                         str(np.rad2deg(self.params.yBar)) + '-' + str(np.rad2deg(self.params.y0)) + \
                         '_' + str(self.params.mask_thresh) + '_yawERR.mat',
                         mdict={'yawErrTh' + str(int(self.params.mask_thresh * 10)): yawError})
            scio.savemat('mat/' + str(i) + 'turbs_' + str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '_refERRnorm.mat',
                         mdict={'refERRnorm' + str(int(self.params.mask_thresh * 10)): ref_err_field_norm})
            scio.savemat('mat/' + str(i) + 'turbs_' + str(self.params.vBar) + '-' + str(self.params.v0) + '_' +
                         str(np.rad2deg(self.params.dBar)) + '_ERRnorm.mat',
                         mdict={'ERRnorm' + str(int(self.params.mask_thresh * 10)): err_field_norm})
        print('max coverage: ' + str(np.amax(mask_per)) + '%')
        print('min coverage: ' + str(np.amin(mask_per)) + '%')
        print('average coverage: ' + str(np.average(mask_per)) + '%')
        plt.subplots_adjust(left=0.05,
                            bottom=0.15,
                            right=0.95,
                            top=0.83,
                            wspace=0.27,
                            hspace=0.19)

        def update_plot(val):

            idx = int(round(sld.val))
            sld.valtext.set_text('iteration ' + '{}'.format(idx))
            spdErrAx.clear()
            yawErrAx.clear()
            errMapAx.clear()
            maskAx.clear()
            errMapAx.contourf(X, Y, err_field[idx], v, cmap='gnuplot2')
            maskAx.scatter(x, y, masks[idx], color='black')
            maskAx.text(-220, 290, str(round(V_k[idx], 2)) + ' m/s, ' + \
                        str(round(np.rad2deg(Y_k[idx]), 2)) + '\N{DEGREE SIGN}', fontsize=fontsize)
            maskAx.arrow(-220, 225, 15 * V_k[idx] * np.cos(Y_k[idx]),
                         -15 * V_k[idx] * np.sin(Y_k[idx]), head_width=50, head_length=50)
            maskAx.set_xlabel('meters', fontsize=fontsize)
            for tick in maskAx.xaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            for tick in maskAx.yaxis.get_major_ticks():
                tick.label.set_fontsize(fontsize)
            spdErrAx.plot(iterations, speedError)
            spdErrAx.axvline(idx, color='red')
            spdErrAx.set_title('speed error')
            yawErrAx.plot(iterations, yawError)
            yawErrAx.axvline(idx, color='red')
            yawErrAx.set_title('yaw error')

            for coord, turbine in self.WF.turbines:
                x_0 = coord.x1 + np.sin(Y_k[idx]) * turbine.rotor_radius
                x_1 = coord.x1 - np.sin(Y_k[idx]) * turbine.rotor_radius
                y_0 = coord.x2 - np.cos(Y_k[idx]) * turbine.rotor_radius
                y_1 = coord.x2 + np.cos(Y_k[idx]) * turbine.rotor_radius
                X_0 = coord.x1 + np.sin(self.params.yBar) * turbine.rotor_radius
                X_1 = coord.x1 - np.sin(self.params.yBar) * turbine.rotor_radius
                Y_0 = coord.x2 - np.cos(self.params.yBar) * turbine.rotor_radius
                Y_1 = coord.x2 + np.cos(self.params.yBar) * turbine.rotor_radius

                errMapAx.plot([x_0, x_1], [y_0, y_1], color='lime', linewidth=1)
                maskAx.plot([x_0, x_1], [y_0, y_1], color='red', linewidth=1)
                errMapAx.plot([X_0, X_1], [Y_0, Y_1], color='white', linewidth=1)
                maskAx.plot([X_0, X_1], [Y_0, Y_1], color='black', linewidth=1)
                plt.draw()
            print(sum(masks[idx]))

        sld.on_changed(update_plot)

        def animate(frame, *fargs):
            # if it's not at the end, increment the slider value
            if sld.val < sld.valmax - 1:
                temp = sld.val
                sld.set_val(temp + 1)
            else:
                # if it's at the end, set it to the beginning
                sld.set_val(sld.valmin)

        # set the animate function to the FuncAnimation function for animation
        if ANIMATE:
            an = anim.FuncAnimation(f, animate, interval=100, frames=sld.valmax * 5)
            # render to video. to make it play faster, increase fps
            if FILE:
                an.save(FILE + '.mp4', fps=7, dpi=300)
        if SHOW:
            plt.show()

    # params
    ## @brief  A class used to pass parameters to various functions
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
        ## y0 serves as the initial yaw estimate
        y0 = np.deg2rad(0.0)
        ## speed epsilon for calculating the sensitivity matrix
        epSpeed = 0.001  # speed epsilon (ev)
        ## direction epsilon for calculating the sensitivity matrix
        epDir = 0.0001  # direction epsilon (ed)
        ## yaw epsilon for calculating sensitivity matrix
        epYaw = 0.0001
        ## speed error threshold for stopping iterations
        spErrMax = 0.1  # speed error threshold
        ## direction error threshold for stopping iterations
        dirErrMax = 0.01  # direction error threshold
        ## yaw error threshold for stopping iterations
        yawErrMax = 0.01
        ## iteration threshold for stopping iterations
        iterMax = 40  # iteration threshold
        ## vBar serves as the 'actual' wind speed
        vBar = 8.0
        ## dBar serves as the 'actual' wind direction
        dBar = np.deg2rad(0.0) # actual wind direction
        ## yBar serves as the 'actual' yaw direction
        yBar = np.deg2rad(0.0)
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
        ## a range of yaws: [min yaw err, max yaw err, step]
        Yrange = [-2.0, 2.0, 0.25]
