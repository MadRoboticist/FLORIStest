## \file pathPlan.py
# This file contains functions which handle UAV path planning using Floris data.
#
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from math import pi as PI
import math
from copy import deepcopy
from matplotlib import animation as anim

## \class d
# An enumeration of the cardinal directions
# starting from East and moving counterclockwise.
#
# This enumeration can be used to generate compass
# angles, e.g. NW = 3x45 degrees or 3xPI/4
class d(enumerate):
    ## @var E
    # East = 0
    E = 0
    ## @var NE
    # Northeast = 1
    NE = 1
    ## @var N
    # North = 2
    N = 2
    ## @var NW
    # Northwest = 3
    NW = 3
    ## @var W
    # West = 4
    W = 4
    ## @var SW
    # Southwest = 5
    SW = 5
    ## @var S
    # South = 6
    S = 6
    ## @var SE
    # Southeast = 7
    SE = 7


## a structure which holds information for each node on
# a discretized transition map
class _X_map_node:
    ## Class constructor
    def __init__(self):
        ## @var Xitions
        # A list of the 8 possible transitions
        # from the corresponding node
        self.Xitions = [self._Transition() for i in range(8)]
        ## @var GPS
        # The GPS coordinates [lat,lon] of the corresponding node
        self.GPS = [0, 0]
        ## @var idx
        # The coordinates [x, y] of the node
        self.idx = [0, 0]
        # set aspect ratio of planner plot, initialize to 1
        self._aspect = 1

    ## a structure which holds information for each possible transition from
    #   one position on the map to another
    class _Transition:
        ## Class constructor
        def __init__(self):
            ## @var thetaCost
            # A cost value
            # for wind direction vs. transition heading
            self.thetaCost = 0.0
            ## @var dSscore
            # A score/cost value
            # for the rate of change of the total score
            # resulting from the corresponding transition
            self.dSscore = 0.0

## \class PathPlanner
#   @brief Contains all path planning functions and parameters
class PathPlanner:

    ## @var minDist2turbine
    # sets a minimum distance from the center of a turbine
    # which serves as a boundary for the UAV
    minDist2turbine = 160

    ## Class constructor
    #   @param vman a VisualManager object
    def __init__(self, vman, Xbar=None, YAW=False):
        ## @var YAW
        #  A boolean for if yaw estimates are being calculated
        self.YAW = YAW
        ## @var vman
        #   A VisualizationManager object
        self.vman = deepcopy(vman)
        ## @var params
        #   A local copy of vman.params
        self.params = self.vman.params
        ## @var WF
        #  a local copy of vman.WF
        self.WF = deepcopy(self.vman.WF)
        # time to calculate Xbar
        if not self.YAW:
            self.WF.set_incoming(self.params.vBar, np.rad2deg(self.params.dBar))
        else:
            self.WF.set_vy(self.params.vBar, np.rad2deg(self.params.yBar))
        ## @var Xbar
        #   this is the 'actual' u_field
        if Xbar is None:
            print("no xBar provided, using xBar from FLORIS")
            if not self.YAW:
                self.WF.set_incoming(self.params.vBar, np.rad2deg(self.params.dBar))
            else:
                self.WF.set_vy(self.params.vBar, np.rad2deg(self.params.yBar))
            ff_bar = deepcopy(self.WF.u_field)
            ff = deepcopy(ff_bar)
            self.Xbar = np.vstack(ff.flatten())  # 1xn
        elif Xbar.ndim == 2:
            print("Static xBar provided")
            ff_bar = deepcopy(Xbar)
            ff = deepcopy(ff_bar)
            self.Xbar = np.vstack(ff.flatten())
        else:
            print("Time varying xBar provided")
            TV = True
            self.Xbar = deepcopy(Xbar)

        ## @var X_map
        #   A map of all possible transitions on the map and their respective scores
        self.X_map = [[_X_map_node() for x in range(self.vman.grid_resolution[1])]
                 for y in range(self.vman.grid_resolution[0])]
        ## @var X
        #   An X grid mesh of the map's 'GPS' coordinates
        self.X = self.WF.x_mesh
        ## @var Y
        #   A Y grid mesh of the map's 'GPS' coordinates
        self.Y = self.WF.y_mesh
        # set the incoming speed and direction to the estimates
        ## @var error
        # holds a record of the direction and velocity errors for plotting
        self.error = list()

        if not self.YAW:
            self.WF.set_incoming(self.params.v0, np.rad2deg(self.params.d0))
            ## the local copy of the sensitivity matrix
            self.sens_mat = self.vman.calcSensitivityMatrix(self.params.v0, self.params.d0)
            # append the initial error
            self.error.append([self.params.vBar - self.params.v0,
                               self.params.dBar - self.params.d0])
        else:
            self.WF.set_vy(self.params.v0, np.rad2deg(self.params.y0))
            self.sens_mat = self.vman.calcSensitivityMatrix_vy(self.params.v0, self.params.y0)
            # append the initial error
            self.error.append([self.params.vBar-self.params.v0,
                               self.params.yBar-self.params.y0])
        # calculate the score map
        self.calcScoreMap()
        # calculate the cost map
        self._init_cost_map()
        self.hist = list()
        self.IDMhistory_max = list()
        self.IDMhistory_min = list()
        self.IDMhistory_avg = list()
        self.mask_size_history = list()
        self.learn_length = 1000
        self.learn_alpha = 0.5
        self.learn_gamma = -0.5
        self.learn_thresh = 0.5
        self.max_score_box = 2
        self.steps_planned = 0
        self.percent_plan = 0.75
        self.learn_path_length = 200
        self.MAX_BONUS = 100
        self.steps_taken = 0

    ## a method to adjust idx based on a cardinal direction
    def _shiftVals(self, val):
        value = {
            d.E: [1, 0],    # one right, zero up
            d.NE: [1, 1],   # one right, one up
            d.N: [0, 1],    # zero right, one up
            d.NW: [-1, 1],  # one left, one up
            d.W: [-1, 0],   # one left, zero up
            d.SW: [-1, -1], # one left, one down
            d.S: [0, -1],   # zero left, one down
            d.SE: [1, -1]   # one right, one down
        }
        return value.get(val)

    ## calculates the euclidean distance between two points
    def euclidDist(self, p1, p2):
        # sqrt((Y2-Y1)^2+(X2-X1)^2)
        return math.sqrt(((p2[0]-p1[0])**2)+((p2[1]-p1[1])**2))

    ## a method to update the transition/cost map
    def _init_cost_map(self):
        # iterate over the entire range of the discretized map
        for i in range(self.vman.grid_resolution[0]):
            for j in range(self.vman.grid_resolution[1]):
                # record the index of the node
                # so that it can easily be retrieved if non-indexed
                # search methods are used
                self.X_map[i][j].idx = [i, j]
                # record the 'GPS' location of the node
                self.X_map[i][j].GPS = [self.X[i][j], self.Y[i][j]]
                # iterate through the 8 possible transitions at each node
                for k in range(8):
                    # calculate the cost of moving through the wind at the transition heading
                    self.X_map[i][j].Xitions[k].thetaCost = k*np.deg2rad(45)-self.params.d0
                    # adjust the cost to a positive value 0-PI
                    if self.X_map[i][j].Xitions[k].thetaCost > PI:
                        self.X_map[i][j].Xitions[k].thetaCost = \
                            abs(k*np.deg2rad(45)-self.params.d0-2*PI)
                    # test the map indices after the transition
                    [I, J] = self._shiftVals(k)
                    # if the index is out of bounds, this will throw an exception
                    try:
                        # compute the first derivative of score with respect to transition
                        dS = self.score_map[i + I][j + J] - self.score_map[i][j]
                        # check to make sure the number is real
                        if math.isnan(dS):
                            # if not, then the move is not possible
                            self.X_map[i][j].Xitions[k].dSscore = None
                        else:
                            # otherwise, take a copy of the value
                            self.X_map[i][j].Xitions[k].dSscore = deepcopy(dS)
                    # if the index is out of bounds...
                    except:
                        # then the move is not possible
                        self.X_map[i][j].Xitions[k].dSscore = None

    def _update_cost_map(self, UAV):
        # iterate over the entire range of the discretized map
        for i in range(self.vman.grid_resolution[0]):
            for j in range(self.vman.grid_resolution[1]):
                # record the index of the node
                # so that it can easily be retrieved if non-indexed
                # search methods are used
                self.X_map[i][j].idx = [i, j]
                # record the 'GPS' location of the node
                self.X_map[i][j].GPS = [self.X[i][j], self.Y[i][j]]
                # iterate through the 8 possible transitions at each node
                for k in range(8):
                    # calculate the cost of moving through the wind at the transition heading
                    self.X_map[i][j].Xitions[k].thetaCost = k*np.deg2rad(45)-UAV.d0
                    # adjust the cost to a positive value 0-PI
                    if self.X_map[i][j].Xitions[k].thetaCost > PI:
                        self.X_map[i][j].Xitions[k].thetaCost = \
                            abs(k*np.deg2rad(45)-UAV.d0-2*PI)
                    # test the map indices after the transition
                    [I, J] = self._shiftVals(k)
                    # if the index is out of bounds, this will throw an exception
                    try:
                        # compute the first derivative of score with respect to transition
                        dS = self.score_map[i + I][j + J] - self.score_map[i][j]
                        # check to make sure the number is real
                        if math.isnan(dS):
                            # if not, then the move is not possible
                            self.X_map[i][j].Xitions[k].dSscore = None
                        else:
                            # otherwise, take a copy of the value
                            self.X_map[i][j].Xitions[k].dSscore = deepcopy(dS)
                    # if the index is out of bounds...
                    except:
                        # then the move is not possible
                        self.X_map[i][j].Xitions[k].dSscore = None

    ## a method to calculate the score map from the current sensitivity matrix
    def calcScoreMap(self):
        # compute the normalized df/dd column of the sensitivity matrix ||df/dd||
        Zd = abs(self.sens_mat[1]) / np.amax(abs(self.sens_mat[1]))
        # compute the normalized df/dv column of the sensitivity matrix ||df/dv||
        Zv = abs(self.sens_mat[0]) / np.amax(abs(self.sens_mat[0]))
        ## @var score_map
        # the score map is computed as ||df/dd||*||df/dv||
        self.score_map = Zd * Zv
        #self.score_map = Zd
        # iterate through the entire score map
        for i in range(self.vman.grid_resolution[0]):
            for j in range(self.vman.grid_resolution[1]):
                # create a no-fly zone around each turbine
                for coord, turbine in self.WF.turbines:
                    # compute if the euclidean distance is less
                    # than the minimum distance threshold minDist2turbine
                    if self.euclidDist([self.X[i][j],
                                         self.Y[i][j]],
                                        [coord.x1, coord.x2]) < \
                            self.minDist2turbine:
                        # if it is, make it untraversable
                        self.score_map[i][j] = None
        # calculate a wave map from this score map

    ## calculates a wave map using the highest score on the score map
    def _calcWaveMap(self, UAV, plan_mask=None):
        # set the max unobtainably low (since our range is 0-1)
        if plan_mask is None:
            plan_mask = UAV.plan_mask
        max = -100
        idx = [0,0]
        # set the wave_map to all ones
        UAV.wave_map = [[1 for i in range(self.vman.grid_resolution[1])]
                        for j in range(self.vman.grid_resolution[0])]
        # iterate through the entire map
        for i in range(self.vman.grid_resolution[0]):
            for j in range(self.vman.grid_resolution[1]):
                # if we aren't using a UAV for this calculation
                # (such as for the initial wave map calculation)
                # the program will throw an exception
                try:
                    # we test the score using the UAV's path mask
                    if self.score_map[i][j]*plan_mask[i][j]*UAV.path_mask[i][j] > max:
                        max = self.score_map[i][j]
                        idx = [i, j]

                except:
                    if self.score_map[i][j] > max:
                        max = self.score_map[i][j]
                        idx = [i, j]
        # just for shorthand
        x = idx[0]
        y = idx[1]
        # iterate through the map again
        for i in range(self.vman.grid_resolution[0]):
            for j in range(self.vman.grid_resolution[1]):
                # set the value of the wave map to the euclidean
                # distance from the maximum score
                try:
                    UAV.wave_map[i][j] = 1/(self.euclidDist([x, y], [i, j])*self.euclidDist([x, y], [i, j]))
                except:
                    UAV.wave_map[i][j] = 1
        UAV.max_wave_idx = idx

    ## A method to convert an index to a GPS point
    def _findGPSindex(self, GPS, check=True):
        # Figure out how 'wide' each range is
        XleftSpan = np.amax(self.X) - np.amin(self.X)
        XrightSpan = self.vman.grid_resolution[0]
        YleftSpan = np.amax(self.Y) - np.amin(self.Y)
        YrightSpan = self.vman.grid_resolution[1]

        # Convert the left range into a 0-1 range (float)
        idx = int((GPS[0] - np.amin(self.X))*XrightSpan / float(XleftSpan))
        idy = int((GPS[1] - np.amin(self.Y))*YrightSpan / float(YleftSpan))
        # Is it a valid GPS point?
        if (math.isnan(self.score_map[idx][idy]) or idx < 0 or idy < 0) and check:
            # if not, let the user know
            print("Error: GPS point is out of bounds: "+str(GPS))
            return [None, None]
        else:
            # otherwise, send back the coordinates
            return [idx, idy]

    ## a method to convert a GPS point to an index
    def _getGPSfromIDX(self, idx):
        # we already have indexed lists of these values
        X = self.X[idx[0]][idx[1]]
        Y = self.Y[idx[0]][idx[1]]
        return [X, Y]

    # greedyPath
    ## A function which populates a UAV's list with a greedy path
    # with the given number of steps in the UAV's planning horizon
    #
    # @param UAV The UAV whose path is to be generated
    # @param UAV.plan_horizon the number of steps to plan ahead
    # @param recalc (optional) the number of moves before the score map is recalculated
    # @return UAV.IDXplan a path of indices through which the UAV will travel
    # @return UAV.GPSplan a path of GPS points throught which the UAV will travel
    def greedyPath(self, UAV):
        plan_horizon = self._init_planner(UAV)
        # if the GPS point was valid, the indices should exist
        if UAV.idx:
            # calculate the current score map
            self.calcScoreMap()
            self._update_cost_map(UAV)
            # let the iteration begin
            for i in range(plan_horizon):
                # just for shorthand notation
                x = UAV.IDXplan[i][0]
                y = UAV.IDXplan[i][1]

                # if the indices are invalid, an exception will be thrown
                try:
                    # start by updating the UAV's score mask for the current position
                    UAV.update_mask(UAV.plan_mask, [x, y])
                    # calculate the new wave map with the updated mask
                    self._calcWaveMap(UAV)
                    # calculate the next step
                    max, dir = self._greedyStep(UAV)
                    UAV.planner.steps_planned = UAV.planner.steps_planned + 1

                except:
                    # if we run off the map, or into a no-fly zone, end the run
                    print("problem occurred")
                    return
                # add the max to the UAV's score
                UAV.score += max
                # add the current dS to the list for d2S calculation
                UAV.dSplan.append(self.X_map[x][y].Xitions[dir].dSscore)
                # update the UAV's heading
                UAV.plan_heading.append(dir)
                # update idx based on the direction of travel
                [I, J] = self._shiftVals(dir)
                # UAV.idx = [x+I,y+J]
                # add the index to the IDXplan
                UAV.IDXplan.append([UAV.IDXplan[i][0]+I, UAV.IDXplan[i][1]+J])
                # add the GPS point to the UAV's path
                UAV.GPSplan.append(self._getGPSfromIDX(UAV.IDXplan[i+1]))
                # uncomment the following to plot each step of the planner
            # print(UAV.GPSplan)
        # the starting point was invalid
        else:
            print("Aborting. Invalid start point.")

    ## a method which calculates a single step in a greedy path
    def _greedyStep(self, UAV):
        # just for shorthand notation
        x = UAV.IDXplan[len(UAV.IDXplan)-1][0]
        y = UAV.IDXplan[len(UAV.IDXplan)-1][1]

        # first check to make sure the transition is valid
        if x >= self.vman.grid_resolution[0] or \
            y >= self.vman.grid_resolution[1] or \
            x < 0 or y < 0:
            return None
        # set an unobtainably low maximum value
        max = -10000000
        # set the direction to zero just for starters
        dir = 0
        # if there isn't a transition list for this node
        # then the node doesn't exist, and any move is invalid
        try:
            # just for shorthand
            Xition = self.X_map[x][y].Xitions
        except:
            # if that doesn't work, something is wrong. end the run
            return None
        # iterate through the 8 possible transitions
        for i in range(8):
            # check the indices of the transition
            [I, J] = self._shiftVals(i)
            # if the transition index is invalid an exception will be thrown
            try:
                # pull the score from the score map for the transition node
                score = self.score_map[x+I][y+J]
                # pull the UAV's mask value for the transition node
                mask = UAV.plan_mask[x+I][y+J]
                # pull the wave map value for the transition node
                wave = UAV.wave_map[x+I][y+J]
                # calculate the cost of changing the UAV's heading
                headCost = abs(i * np.deg2rad(45) -
                               (UAV.plan_heading[len(UAV.plan_heading) - 1]) *
                               np.deg2rad(45))

                # normalize to a value between -PI and PI
                if headCost > PI:
                    headCost = abs(headCost - 2 * PI)
                # calculate the 2nd derivative of the score wrt the transition
                d2S = Xition[i].dSscore-UAV.dSplan[len(UAV.dSplan)-1]
                # calculate the reward of the move
                reward = (UAV.scoreWt * score +
                          UAV.dSwt * Xition[i].dSscore +
                          UAV.d2Swt * d2S -
                          UAV.windWt * Xition[i].thetaCost / PI -
                          UAV.headWt * headCost / PI +
                          UAV.waveWt * wave +
                          UAV.timeWt * len(UAV.IDXplan)) * \
                         mask
                # we are trying to maximize the reward
                if reward > max and \
                    x+I < self.vman.grid_resolution[0] and \
                    y+J < self.vman.grid_resolution[1] and \
                    x+I >= 0 and y+J >= 0:

                    # if the transition is valid and has the highest reward
                    # then track the reward score and the direction of the transition
                    max = reward
                    dir = i
            # if the transition is invalid, skip it
            except:
                pass
        # return the decision and the resulting score
        return max, dir

    def printAVGs(self):
        print("IDM max: "+str(np.average(self.IDMhistory_max)))
        print("IDM min: " + str(np.average(self.IDMhistory_min)))
        print("IDM avg: "+str(np.average(self.IDMhistory_avg)))
        print("mask max: "+str(np.amax(self.mask_size_history)))
        print("mask min: "+str(np.amin(self.mask_size_history)))
        print("mask avg: "+str(np.average(self.mask_size_history)))

    # plotHistory
    ## A function which plots a UAV's history
    #
    # @param UAV a UAV object or list of UAV objects to plot the history of
    # @param filename a filename for the video WITHOUT THE FILE EXTENSION
    #       omit for no output video
    # @param plot boolean- to show the plot or not to show the plot
    #
    def plotHistory(self, UAV, filename=None, plot=True):
        if self.YAW:
            self.plotHistory_vy(UAV, filename, plot)
            return
        self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'lime', 'orange','purple']
        # variable which tells us if we are looking at the Wave or the Score map
        self._map_sel = 0
        # boolean which tells us if the animation is playing
        if filename is None:
            self._play = False
        else:
            self._play = True
        fontsize = 14
        # create a figure
        f = plt.figure(figsize=(10, 9))
        # separate it into a grid
        gs = f.add_gridspec(2, 3)
        # put the error plots in the leftmost figures
        ax_v = f.add_subplot(gs[0, 0], title='speed error')
        ax_d = f.add_subplot(gs[1, 0])
        # the path plot takes up the right 2/3 of the figure
        axbig = f.add_subplot(gs[0:, 1:])
        axbig.set_aspect('equal')
        # plot the velocity error
        ax_v.plot([i for i in range(len(self.error))], [i[0] for i in self.error])
        # format x,y and c for scatterplotting
        x = deepcopy(self.X)
        y = deepcopy(self.Y)
        if type(UAV) == list:
            c = deepcopy(UAV[0].planner.hist[1][0])
        else:
            c = deepcopy(UAV.planner.hist[1][0])
        # plot the score map
        cont = axbig.scatter(x=x.flatten(), y=y.flatten(), c=c.flatten(), s=15, cmap='gnuplot2')
        divider = make_axes_locatable(axbig)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = plt.colorbar(cont, cax=cax)
        cb.set_label('IDM score',fontsize=fontsize)
        cb.ax.tick_params(labelsize=fontsize)
        for tick in axbig.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        for tick in axbig.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        axbig.set_xlabel('meters', fontsize=fontsize)
        # plot the direction error
        axbig.clear()
        if type(UAV) == list:
            for ind, uav in enumerate(UAV):
                ax_v.plot([i for i in range(len(uav.planner.error))],
                          [i[0] for i in uav.planner.error],
                          color=self.colors[ind % len(self.colors)])
                ax_d.plot([i for i in range(len(uav.planner.error))],
                          [np.rad2deg(i[1]) for i in uav.planner.error],
                          color=self.colors[ind % len(self.colors)])
                axbig.scatter(x=x[uav.minX:uav.maxX, uav.minY:uav.maxY].flatten(),
                              y=y[uav.minX:uav.maxX, uav.minY:uav.maxY].flatten(),
                              c=uav.planner.hist[0][self._map_sel][uav.minX:uav.maxX, uav.minY:uav.maxY].flatten(),
                              s=15, cmap='gnuplot2')
        else:
            ax_v.plot([i for i in range(len(UAV.planner.error))],
                      [i[0] for i in UAV.planner.error], color='black')
            ax_d.plot([i for i in range(len(UAV.planner.error))],
                      [np.rad2deg(i[1]) for i in UAV.planner.error],color='black')
            axbig.scatter(x=x[UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                          y=y[UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                          c=UAV.planner.hist[0][self._map_sel][UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                          s=15, cmap='gnuplot2')
        for coord, turbine in self.WF.turbines:
            x_0 = coord.x1 + np.sin(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            x_1 = coord.x1 - np.sin(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            y_0 = coord.x2 - np.cos(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            y_1 = coord.x2 + np.cos(np.deg2rad(turbine.yaw_angle)) * turbine.rotor_radius
            axbig.plot([x_0, x_1], [y_0, y_1], linewidth=1, color='black')
        # the plot's title holds a lot of info
        if type(UAV) == list:
            plt.suptitle("Wind Speed and Direction Estimates with a UAV swarm for Sensing\n" +
                         'Actual Speed: ' + str(self.params.vBar) +
                         ', Actual Direction: ' + str(np.rad2deg(self.params.dBar)) + '\N{DEGREE SIGN}' +
                         '\nPlanning Horizon: ' + str(UAV[0].plan_horizon) +
                         '          Moves until recalculation: ' + str(UAV[0].moves2recalc) +
                         '          Memory: ' + str(UAV[0].patrolMax) + ' nodes')
        else:
            plt.suptitle("Wind Speed and Direction Estimates with a UAV swarm for Sensing\n" +
                         'Actual Speed: ' + str(self.params.vBar) +
                         ', Actual Direction: ' + str(np.rad2deg(self.params.dBar)) + '\N{DEGREE SIGN}' +
                         '\nPlanning Horizon: ' + str(UAV.plan_horizon) +
                         '          Moves until recalculation: ' + str(UAV.moves2recalc) +
                         '          Memory: ' + str(UAV.patrolMax) + ' nodes')
        # set plot and axis titles for the error plots
        ax_d.set_ylabel('direction error (\N{DEGREE SIGN})')
        ax_d.set_xlabel('# of recalculations')
        ax_v.set_ylabel('speed error (m/s)')
        # adjust how the plots fill the figure
        plt.subplots_adjust(left=0.05,
                            bottom=0.15,
                            right=0.95,
                            top=0.83,
                            wspace=0.27,
                            hspace=0.19)
        # slider axis
        sld_ax = plt.axes((0.2, 0.02, 0.56, 0.02))
        # create the slider
        LEN = 1000000000000000
        if type(UAV)==list:
            for uav in UAV:
                LEN = np.amin([len(uav.planner.hist), LEN])
        else:
            LEN = len(UAV.planner.hist)
        sld = Slider(sld_ax, 'moves', 0, LEN - 1, valinit=0)
        # set the initial slider text
        sld.valtext.set_text('move 0')
        # button axis
        btn_ax = plt.axes([0.85, 0.925, 0.125, 0.05])
        # create a button which toggles between score map and wave map
        btn = Button(btn_ax, 'Show Wave')
        # button axis
        btn2_ax = plt.axes([0.85, 0.01, 0.125, 0.05])
        # create a button which plays/pauses the animation
        if filename is None:
            btn2 = Button(btn2_ax, 'Play')
        else:
            btn2 = Button(btn2_ax, 'Pause')

        # function for the map button
        def map_btn(event):
            if self._map_sel == 0:  # currently showing the score map
                self._map_sel = 3  # change to showing the wave map
                btn.label.set_text('Show Score')  # update button text
            else:
                self._map_sel = 0  # otherwise change to showing the score map
                btn.label.set_text('Show Wave')  # update button text
            update_plot(0)  # and update the plot
        # set the above function for the map button
        btn.on_clicked(map_btn)

        # function for the play button
        def play_btn(event):
            if self._play:  # currently playing
                self._play = False  # pause it
                btn2.label.set_text("Play")  # update button text
            else:
                self._play = True  # otherwise, play the animation
                btn2.label.set_text("Pause")  # update button text
        # set the above function for the play button
        btn2.on_clicked(play_btn)

        # function which updates the plot
        def update_plot(val):
            # discretize the slider to integer values
            idx = int(round(sld.val))
            # set the slider text
            sld.valtext.set_text('move ' + '{}'.format(idx))
            # clear the plots
            axbig.clear()
            ax_v.clear()
            ax_d.clear()
            # check map selection
            if self._map_sel == 0:  # it's the score map
                btn.label.set_text('Show Wave')  # set the button text
            else:  # it's the wave map
                btn.label.set_text('Show Score')  # set the button text\
            # set the subplot labels
            ax_d.set_ylabel('direction error (\N{DEGREE SIGN})')
            ax_d.set_xlabel('# of recalculations')
            ax_v.set_ylabel('speed error (m/s)')
            # plot the map
            if type(UAV) == list:
                for inx, uv in enumerate(UAV):
                    axbig.scatter(x=x[uv.minX:uv.maxX, uv.minY:uv.maxY].flatten(),
                                  y=y[uv.minX:uv.maxX, uv.minY:uv.maxY].flatten(),
                                  c=uv.planner.hist[idx][self._map_sel][uv.minX:uv.maxX,
                                                                        uv.minY:uv.maxY].flatten(),
                                  s=15, cmap='gnuplot2')
                    # plot the plan
                    axbig.plot([i[0] for i in uv.planner.hist[idx][2]],
                               [i[1] for i in uv.planner.hist[idx][2]],
                               linewidth=1.0, color=self.colors[inx % len(self.colors)])
                    # plot the actual path
                    axbig.plot([i[0] for i in uv.planner.hist[idx][1]],
                               [i[1] for i in uv.planner.hist[idx][1]],
                               linewidth=1.0, color=self.colors[inx % len(self.colors)])
                    for tic in axbig.xaxis.get_major_ticks():
                        tic.label.set_fontsize(fontsize)
                    for tic in axbig.yaxis.get_major_ticks():
                        tic.label.set_fontsize(fontsize)
                    # for shorthand, UAV's location on the map
                    try:
                        _x = uv.planner.hist[idx][1][len(uv.planner.hist[idx][1])-1][0]
                        _y = uv.planner.hist[idx][1][len(uv.planner.hist[idx][1])-1][1]
                        # plot the UAV, a big red circle
                        axbig.plot(_x, _y, marker='o', markersize=10, color=self.colors[inx % len(self.colors)])
                    except:
                        print("couldn't draw UAV")
                        pass

                    # plot the velocity error
                    ax_v.plot([i for i in range(len(uv.planner.error))], [i[0] for i in uv.planner.error],
                              color=self.colors[inx % len(self.colors)])
                    # plot the direction error
                    ax_d.plot([i for i in range(len(uv.planner.error))], [np.rad2deg(i[1]) for i in uv.planner.error],
                              color=self.colors[inx % len(self.colors)])
            else:
                axbig.scatter(x=x[UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                              y=y[UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                              c=UAV.planner.hist[idx][self._map_sel][UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                              s=15, cmap='gnuplot2')
                # plot the plan
                axbig.plot([i[0] for i in UAV.planner.hist[idx][2]],
                           [i[1] for i in UAV.planner.hist[idx][2]],
                           linewidth=1.0, color='lime')
                # plot the actual path
                axbig.plot([i[0] for i in UAV.planner.hist[idx][1]],
                           [i[1] for i in UAV.planner.hist[idx][1]],
                           linewidth=1.0, color='red')
                for tic in axbig.xaxis.get_major_ticks():
                    tic.label.set_fontsize(fontsize)
                for tic in axbig.yaxis.get_major_ticks():
                    tic.label.set_fontsize(fontsize)
                # for shorthand, UAV's location on the map
                try:
                    _x = UAV.planner.hist[idx][1][len(UAV.planner.hist[idx][1]) - 1][0]
                    _y = UAV.planner.hist[idx][1][len(UAV.planner.hist[idx][1]) - 1][1]
                    # plot the UAV, a big red circle
                    axbig.plot(_x, _y, marker='o', markersize=10, color='red')
                except:
                    print("couldn't draw UAV")
                    pass
                # plot the velocity error
                ax_v.plot([i for i in range(len(UAV.planner.error))],
                          [i[0] for i in UAV.planner.error], color='black')
                # plot the direction error
                ax_d.plot([i for i in range(len(UAV.planner.error))],
                          [np.rad2deg(i[1]) for i in UAV.planner.error], color='black')
            # plot the turbines
            for coor, turbin in self.WF.turbines:
                x1_0 = coor.x1 + np.sin(np.deg2rad(turbin.yaw_angle)) * turbin.rotor_radius
                x1_1 = coor.x1 - np.sin(np.deg2rad(turbin.yaw_angle)) * turbin.rotor_radius
                y1_0 = coor.x2 - np.cos(np.deg2rad(turbin.yaw_angle)) * turbin.rotor_radius
                y1_1 = coor.x2 + np.cos(np.deg2rad(turbin.yaw_angle)) * turbin.rotor_radius
                axbig.plot([x1_0, x1_1], [y1_0, y1_1], linewidth=1, color='black')
            # show the plot
            plt.draw()
        # set the update function for what happens when the slider value changes
        sld.on_changed(update_plot)
        # adjust how the plots fill the figure
        plt.subplots_adjust(left=0.1,
                            bottom=0.14,
                            right=1.0,
                            top=0.84,
                            wspace=0.2,
                            hspace=0.24)

        # function to animate the plot
        def animate(frame, *fargs):
            # check if it is playing
            if self._play:
                # if it's not at the end, increment the slider value
                if sld.val < sld.valmax-1:
                    temp = sld.val
                    sld.set_val(temp + 1)
                else:  # if it's at the end, set it to the beginning
                    sld.set_val(sld.valmin)

        # set the animate function to the FuncAnimation function for animation
        an = anim.FuncAnimation(f, animate, interval=100, frames=LEN)
        # render to video. to make it play faster, increase fps
        if filename is not None:
            an.save(filename+'.mp4', fps=15, dpi=300)
        # show the plot
        if plot:
            plt.show()

    def learnPath(self, UAV):
        # if not plan_horizon:
        plan_horizon = self._init_planner(UAV)
        print("learn horizon: "+str(self.learn_path_length))
        # if the GPS point was valid, the indices should exist
        if UAV.idx:
            # calculate the current score map
            switches = 0
            self.calcScoreMap()
            Qmap = [[[0 for i in range(8)] for x in range(self.vman.grid_res[1])] \
                              for y in range(self.vman.grid_res[0])]
            print("Learning...")
            plan_mask = deepcopy(UAV.path_mask)
            self._calcWaveMap(UAV, plan_mask)
            print("initial max coordinate: "+str(UAV.max_wave_idx))
            for j in range(self.learn_length):
                self.scoreTOT = 0
                plan_mask = deepcopy(UAV.path_mask)
                X = UAV.idx[0]
                Y = UAV.idx[1]
                self.head = 0
                self.action = 0
                reward = 0
                SWITCH = False

                for i in range(self.learn_path_length):
                    # start with a fresh Qmap
                    step = False
                    while not step:
                        maxQ = -1000000000000
                        for m in range(8):
                            if Qmap[X][Y][m] > maxQ:
                                maxQ = Qmap[X][Y][m]
                                self.action = m
                        if (np.random.random_sample()) < self.learn_thresh:
                            self.action = np.random.randint(0, 8)
                        [I, J] = self._shiftVals(self.action)
                        if 0 <= X+I < self.vman.grid_res[0] and 0 <= Y+J < self.vman.grid_res[1] \
                                and not math.isnan(self.score_map[X+I][Y+J]):
                                maxQ = np.amax(Qmap[X+I][Y+J])
                                Xition = self.X_map[X][Y].Xitions
                                # pull the score from the score map for the transition node
                                score = self.score_map[X + I][Y + J]
                                self.scoreTOT = self.scoreTOT + score
                                # pull the UAV's mask value for the transition node
                                mask = plan_mask[X + I][Y + J]*UAV.path_mask[X + I][Y + J]
                                # pull the wave map value for the transition node
                                wave = UAV.wave_map[X + I][Y + J]
                                # calculate the cost of changing the UAV's heading
                                headCost = abs(self.action * np.deg2rad(45) - \
                                               (self.head * np.deg2rad(45)))
                                if X + I == UAV.max_wave_idx[0] and Y + J == UAV.max_wave_idx[1]:
                                    self.scoreTOT = self.scoreTOT + self.MAX_BONUS
                                # normalize to a value between -PI and PI
                                if headCost > PI:
                                    headCost = abs(headCost - 2 * PI)
                                # calculate the reward of the move
                                '''
                                temp_reward = (UAV.scoreWt * score - \
                                               UAV.windWt * Xition[self.action].thetaCost / PI - \
                                               UAV.headWt * headCost / PI + \
                                               UAV.waveWt * wave) -(mask)*UAV.maskSUB
                                               '''
                                temp_reward = self.scoreTOT
                                if not math.isnan(temp_reward):
                                    reward = temp_reward
                                    if X+I == UAV.max_wave_idx[0] and Y+J == UAV.max_wave_idx[1]:
                                        SWITCH = True
                                    tempQmap = Qmap[X][Y][self.action] + \
                                               self.learn_alpha * (reward +
                                                                   self.learn_gamma * maxQ - Qmap[X][Y][self.action])

                                    if not math.isnan(tempQmap):

                                        Qmap[X][Y][self.action] = tempQmap
                                        X = X + I
                                        Y = Y + J
                                        UAV.update_mask(plan_mask, [X, Y])
                                        step = True
                                        self.head = self.action
                                        # start by updating the UAV's score mask for the current position

                                        if SWITCH:
                                            switch_mask = deepcopy(UAV.path_mask)
                                            for boxi in range(self.max_score_box):
                                                for boxj in range(self.max_score_box):
                                                    if np.sqrt(boxi*boxi+boxj*boxj) < self.max_score_box:
                                                        try:
                                                            UAV.update_mask(switch_mask, [X-boxi, Y-boxj])
                                                        except:
                                                            pass
                                                        try:
                                                            UAV.update_mask(switch_mask, [X + boxi, Y + boxj])
                                                        except:
                                                            pass

                                            temp_mask = np.array(switch_mask)*np.array(plan_mask)
                                            self._calcWaveMap(UAV,temp_mask.tolist())
                                            SWITCH = False
                        else:
                            Qmap[X][Y][self.action] = -100000000000

            print("learning complete.")
            print("plan horizon: " + str(plan_horizon))
            print(" Planning...")
            X = UAV.idx[0]
            Y = UAV.idx[1]
            switch=0
            prevMaxQ = np.amin(Qmap)
            prevCoords=[X, Y]
            if len(UAV.GPSpath) == 0:
                plan_len = UAV.init_plan_horizon
            else:
                plan_len = UAV.plan_horizon
            for n in range(plan_len):
                step = False
                while not step and n<plan_len:
                    maxQ = np.amin(Qmap)
                    index = 0
                    for j in range(8):

                        if Qmap[X][Y][j] > maxQ:
                            [I, J] = self._shiftVals(j)
                            if 0 <= X + I < self.vman.grid_res[0] and \
                               0 <= Y + J < self.vman.grid_res[1] and \
                                    X + I != prevCoords[0] and \
                                    Y + J != prevCoords[1]:
                                        if self.euclidDist([X+I, Y+J], prevCoords) < np.sqrt(2):
                                            pass
                                        else:
                                            maxQ = Qmap[X][Y][j]
                                            Qmap[X][Y][j] = Qmap[X][Y][j] - 100
                                            index = j
                    if maxQ < prevMaxQ:
                        plan_len=plan_len-1
                        step = True
                    elif maxQ == prevMaxQ or maxQ == -100000000000:
                        break
                    else:
                        step = True
                if maxQ == prevMaxQ or maxQ == -100000000000:
                    break
                print([maxQ, prevMaxQ, [X, Y]])
                prevMaxQ = maxQ
                # update the UAV's heading
                UAV.plan_heading.append(index)
                # update idx based on the direction of travel
                prevCoords = [X, Y]
                [I, J] = self._shiftVals(index)
                X = X + I
                Y = Y + J
                # UAV.idx = [x+I,y+J]
                # add the index to the IDXplan
                UAV.IDXplan.append([UAV.IDXplan[n][0] + I, UAV.IDXplan[n][1] + J])
                # add the GPS point to the UAV's path
                UAV.GPSplan.append(self._getGPSfromIDX(UAV.IDXplan[n + 1]))
                # uncomment the following to plot each step of the planner
            print(str(len(UAV.IDXplan))+' steps planned')
            self.steps_planned = len(UAV.IDXplan)

    def _init_planner(self, UAV):
        if len(UAV.GPSpath) == 0:
            print("GPS minimum coords")
            GPSmin = self._getGPSfromIDX([UAV.minX, UAV.minY])
            UAV.minX_GPS = GPSmin[0]
            UAV.minY_GPS = GPSmin[1]
            print(GPSmin)
            print("GPS maximum coords")
            GPSmax = self._getGPSfromIDX([UAV.maxX-1, UAV.maxY-1])
            UAV.maxX_GPS = GPSmax[0]
            UAV.maxY_GPS = GPSmax[1]
            print(GPSmax)
            self._aspect = float(GPSmax[1]-GPSmin[1])/(GPSmax[0]-GPSmin[0])
            print("y/x aspect ratio:")
            print(self._aspect)
        # start the planner fresh
        UAV.reset_planner()
        # figure out what the index of the starting point is
        if len(self.hist) < UAV.init_mask_size:
            plan_horizon = UAV.init_plan_horizon
        else:
            plan_horizon = UAV.plan_horizon

        if not UAV.idx:
            # if there isn't already an index, get it from the 'GPS' value
            UAV.idx = self._findGPSindex(UAV.GPS)
            # and append the index to the planned path
            UAV.IDXplan.append(UAV.idx)
        else:
            # otherwise just add the index to the planned path
            UAV.IDXplan.append(UAV.idx)
        return plan_horizon

    # Cardinal-Ordinal Straight Planner
    def COSPlan(self, UAV):
        plan_horizon = self._init_planner(UAV)
        print("plan horizon: "+str(plan_horizon))
        print("planning...")
        if UAV.idx:
            [X, Y]=UAV.idx
            #print([X, Y])
            # calculate the current score map
            self.calcScoreMap()
            dirplan = list()
            UAV.update_mask(UAV.plan_mask, [X, Y])
            max_sum=-10000000000000
            for i in range(8):
                total = self._total_line(UAV, X, Y, i)
                # print([total,i])
                if total == 'nan':
                    pass
                elif total > max_sum:
                    max_sum = total
                    heading = i

            [I, J] = self._shiftVals(heading)
            #print("first heading: "+str(heading))
            step = 0
            while step < plan_horizon:
                max_sum = -1000000000000
                x = X+I
                y = Y+J
                max_coords=[x, y]
                while UAV.minX <= x < UAV.maxX and \
                      UAV.minY <= y < UAV.maxY:
                    if not np.isnan(self.score_map[x][y]):
                        for i in range(8):
                            if i != heading-4 and i != heading+4 \
                              and i != heading:
                                total = self._total_line(UAV, deepcopy(x), deepcopy(y), i)
                                if total == 'nan':
                                    pass
                                elif total > max_sum:
                                    max_sum = total
                                    max_dir = i
                                    max_coords = [x, y]
                            else:
                                total = 0
                        x = x+I
                        y = y+J
                    else:
                        break
                if max_coords == [X, Y]:
                    print("BROKEN: "+str([X, Y, heading, step]))
                    break
                if heading == max_dir:
                    print("broken")
                    break
                while X != max_coords[0] or Y != max_coords[1]:
                    dirplan.append(heading)
                    X = X + I
                    Y = Y + J
                    UAV.IDXplan.append([X, Y])
                    UAV.GPSplan.append(self._getGPSfromIDX([X, Y]))
                    UAV.plan_heading.append(heading)
                    UAV.update_mask(UAV.plan_mask, [X, Y])
                    step = step+1
                heading = max_dir
                [I, J] = self._shiftVals(max_dir)
            print("steps planned: "+str(len(UAV.IDXplan)))
            self.steps_planned = len(UAV.IDXplan)

    def _total_line(self, UAV, X, Y, d):
        [I, J] = self._shiftVals(d)
        X = X+I
        Y = Y+J
        i = 1
        total = 0
        while UAV.minX <= X < UAV.maxX and \
              UAV.minY <= Y < UAV.maxY:
            if not np.isnan(self.score_map[X][Y]):
                total = total + self.checkScore(UAV, [X, Y], None)
                i = i+1
            else:
                if i == 1:
                    return 'nan'
                else:
                    return total
            X = X + I
            Y = Y + J
        if i > 1:
            return total
        else:
            return 'nan'

    def checkScore(self, UAV, coord, path_mask=None):
        X = coord[0]
        Y = coord[1]
        if path_mask is None:
            path_mask = UAV.plan_mask
        try:
            return self.score_map[X][Y] * path_mask[X][Y] - UAV.maskSUB * (1 - UAV.plan_mask[X][Y])
        except:
            return None

    def quadSearchPlan(self, UAV):
        plan_horizon = self._init_planner(UAV)
        print("plan horizon: " + str(plan_horizon))
        print("planning...")
        if UAV.idx:
            [X, Y] = UAV.idx
            # calculate the current score map
            self.calcScoreMap()
            UAV.update_mask(UAV.plan_mask, [X, Y])
            self.avgQuad([X, Y])

    def avgQuad(self, XY):
        X = XY[0]
        Y = XY[1]
        Q1 = 0
        Q2 = 0
        Q3 = 0
        Q4 = 0
        for i in range(X, self.vman.grid_res[0], 1):
            for j in range(Y, self.vman.grid_res[1], 1):
                if not np.isnan(self.score_map[i][j]):
                    Q1 = Q1+self.score_map[i][j]
        for i in range(0, X, 1):
            for j in range(Y, self.vman.grid_res[1], 1):
                if not np.isnan(self.score_map[i][j]):
                    Q2 = Q2+self.score_map[i][j]
        for i in range(0, X, 1):
            for j in range(0, Y, 1):
                if not np.isnan(self.score_map[i][j]):
                    Q3 = Q3+self.score_map[i][j]
        for i in range(X, self.vman.grid_res[0], 1):
            for j in range(0, Y, 1):
                if not np.isnan(self.score_map[i][j]):
                    Q4 = Q4+self.score_map[i][j]
        print(Q1)
        print(Q2)
        print(Q3)
        print(Q4)
        TOT = Q1+Q2+Q3+Q4
        I = (Q1-Q2-Q3+Q4)/TOT
        J = (Q1+Q2-Q3-Q4)/TOT
        if round(I) == 0 and round(J) == 0:
            g = 2
        if not np.isnan(self.score_map[X+round(I)][Y+round(J)]):
            return round(I), round(J)
        elif not np.isnan(self.score_map[X + round(I)][Y]):
            return

    def updateEstimates(self, UAV):
        if self.YAW:
            self.updateEstimates_vy(UAV)
            return
        print("Updating Estimates...")
        print("Prior estimate: " + str([UAV.v0, UAV.d0]))
        # current field estimate
        temp_Xk = deepcopy(self.WF.u_field).flatten()
        # current velocity and direction estimates
        vd_k = [[deepcopy(UAV.v0)],
                [deepcopy(UAV.d0)]]
        # current sensitivity matrix
        temp_sens_mat0 = self.sens_mat[0].flatten()
        temp_sens_mat1 = self.sens_mat[1].flatten()

        # field actual
        if self.Xbar.ndim == 2:
            temp_Xbar = deepcopy(self.Xbar).flatten()
        else: # time varying
            temp_Xbar = deepcopy(self.Xbar[self.steps_taken])
            for i in range(len(UAV.IDXpath)):
                XY = UAV.IDXpath[i]
                temp_Xbar[XY[0]][XY[1]] = UAV.measure[i]

        temp_mask = deepcopy(np.array(UAV.path_mask).flatten())
        # if the node is not masked by the UAV's path, use the estimate
        self.mask_size_history.append(1 - float(np.sum(UAV.path_mask)) /
                                      (self.vman.grid_resolution[0] * self.vman.grid_resolution[1]))
        sens_mat0 = list()
        sens_mat1 = list()
        Xbar = list()
        Xk = list()
        IDM = list()
        temp_score_map = deepcopy(self.score_map).flatten()
        # reduce the matrix
        for i in range(len(temp_mask)):
            if temp_mask[i] != 1:
                Xbar.append(temp_Xbar.flatten()[i])
                Xk.append(temp_Xk[i])
                sens_mat0.append(temp_sens_mat0[i])
                sens_mat1.append(temp_sens_mat1[i])
                IDM.append(temp_score_map[i])

        # convert Xk and Xbar to column vectors
        Xbar = np.vstack(np.array(Xbar))
        Xk = np.vstack(np.array(Xk))
        sens_mat0 = np.vstack(np.array(sens_mat0))
        sens_mat1 = np.vstack(np.array(sens_mat1))
        # calculate the pseudoinverse
        sens_mat_pinv = np.linalg.pinv(np.column_stack([sens_mat0, sens_mat1]))
        # calculate the next estimate: Vdk+1 = Vdk + sens_mat_pinv*(Xbar-Xk)
        vd_kp1 = vd_k + np.matmul(sens_mat_pinv, Xbar - Xk)
        # update the estimates
        UAV.v0 = vd_kp1[0][0]
        UAV.d0 = vd_kp1[1][0]
        # update the wind field estimate
        self.WF.set_incoming(UAV.v0, np.rad2deg(UAV.d0))
        # update the sensitivity matrix
        self.sens_mat = deepcopy(self.vman.calcSensitivityMatrix(UAV.v0, UAV.d0))
        # update the score map
        self.calcScoreMap()
        # keep track of the error signals for plotting
        error = [self.params.vBar - UAV.v0,
                 self.params.dBar - UAV.d0]
        print("Actual: "+str([self.params.vBar, self.params.dBar]))
        print("New estimate: " + str([UAV.v0, UAV.d0]))
        print("Error: [v,theta]=" + str(error))

    ## Update estimates for speed and yaw
    def updateEstimates_vy(self, UAV):
        if not self.YAW:
            self.updateEstimates(UAV)
        print("Updating Estimates...")
        print("Prior estimate: " + str([UAV.v0, UAV.y0]))
        # current field estimate
        temp_Xk = deepcopy(self.WF.u_field).flatten()
        # current velocity and direction estimates
        vy_k = [[deepcopy(UAV.v0)],
                [deepcopy(UAV.y0)]]
        # current sensitivity matrix
        temp_sens_mat0 = self.sens_mat[0].flatten()
        temp_sens_mat1 = self.sens_mat[1].flatten()
        # field actual
        if self.Xbar.ndim == 2:
            temp_Xbar = deepcopy(self.Xbar).flatten()
        else: # time varying
            temp_Xbar = deepcopy(self.Xbar[self.steps_taken])
            for i in range(len(UAV.IDXpath)):
                XY = UAV.IDXpath[i]
                temp_Xbar[XY[0]][XY[1]] = UAV.measure[i]
        temp_mask = deepcopy(np.array(UAV.path_mask).flatten())
        # if the node is not masked by the UAV's path, use the estimate
        self.mask_size_history.append(1 - float(np.sum(UAV.path_mask)) /
                                      (self.vman.grid_resolution[0] * self.vman.grid_resolution[1]))
        sens_mat0 = list()
        sens_mat1 = list()
        Xbar = list()
        Xk = list()
        IDM = list()
        temp_score_map = deepcopy(self.score_map).flatten()
        # reduce the matrix
        for i in range(len(temp_mask)):
            if temp_mask[i] != 1:
                Xbar.append(temp_Xbar.flatten()[i])
                Xk.append(temp_Xk[i])
                sens_mat0.append(temp_sens_mat0[i])
                sens_mat1.append(temp_sens_mat1[i])
                IDM.append(temp_score_map[i])

        # convert Xk and Xbar to column vectors
        Xbar = np.vstack(np.array(Xbar))
        Xk = np.vstack(np.array(Xk))
        sens_mat0 = np.vstack(np.array(sens_mat0))
        sens_mat1 = np.vstack(np.array(sens_mat1))
        # calculate the pseudoinverse
        sens_mat_pinv = np.linalg.pinv(np.column_stack([sens_mat0, sens_mat1]))
        # calculate the next estimate: Vdk+1 = Vdk + sens_mat_pinv*(Xbar-Xk)
        vy_kp1 = vy_k + np.matmul(sens_mat_pinv, Xbar - Xk)
        # update the estimates
        UAV.v0 = vy_kp1[0][0]
        UAV.y0 = vy_kp1[1][0]
        # update the wind field estimate
        self.WF.set_vy(UAV.v0, np.rad2deg(UAV.y0))
        # update the sensitivity matrix
        self.sens_mat = deepcopy(self.vman.calcSensitivityMatrix_vy(UAV.v0, UAV.y0))
        # update the score map
        self.calcScoreMap()
        # keep track of the error signals for plotting
        error = [self.params.vBar - UAV.v0,
                 self.params.yBar - UAV.y0]
        print("Actual: "+str([self.params.vBar, self.params.yBar]))
        print("New estimate: " + str([UAV.v0, UAV.y0]))
        print("Error: [v,theta]=" + str(error))

    ## Gets measurement from time-varying xBar (SOWFA)
    def getMeasure(self, XY):
        print("getting measurement")
        if self.Xbar.ndim == 2:
            measure = self.Xbar[XY[0]][XY[1]]
        else:
            measure = self.Xbar[self.steps_taken][XY[0]][XY[1]]
        print("measurement taken")
        self.steps_taken = self.steps_taken+1
        if self.steps_taken >= self.Xbar.shape[0]:
            print("looping to beginning of data")
            self.steps_taken = 0
        return measure

    # plotHistory_yv
    ## A function which plots a UAV's history
    #
    # @param UAV a UAV object or list of UAV objects to plot the history of
    # @param filename a filename for the video WITHOUT THE FILE EXTENSION
    #       omit for no output video
    # @param plot boolean- to show the plot or not to show the plot
    #
    def plotHistory_vy(self, UAV, filename=None, plot=True):
        if not self.YAW:
            self.plotHistory(UAV, filename, plot)
            return
        self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'lime', 'orange', 'purple']
        # variable which tells us if we are looking at the Wave or the Score map
        self._map_sel = 0
        # boolean which tells us if the animation is playing
        if filename is None:
            self._play = False
        else:
            self._play = True
        fontsize = 14
        # create a figure
        f = plt.figure(figsize=(10, 9))
        # separate it into a grid
        gs = f.add_gridspec(2, 3)
        # put the error plots in the leftmost figures
        ax_v = f.add_subplot(gs[0, 0], title='speed error')
        ax_y = f.add_subplot(gs[1, 0])
        # the path plot takes up the right 2/3 of the figure
        axbig = f.add_subplot(gs[0:, 1:])
        axbig.set_aspect('equal')
        # plot the velocity error
        ax_v.plot([i for i in range(len(self.error))], [i[0] for i in self.error])
        # format x,y and c for scatterplotting
        x = deepcopy(self.X)
        y = deepcopy(self.Y)
        if type(UAV) == list:
            c = deepcopy(UAV[0].planner.hist[1][0])
        else:
            c = deepcopy(UAV.planner.hist[1][0])
        # plot the score map
        cont = axbig.scatter(x=x.flatten(), y=y.flatten(), c=c.flatten(), s=15, cmap='gnuplot2')
        divider = make_axes_locatable(axbig)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cb = plt.colorbar(cont, cax=cax)
        cb.set_label('IDM score', fontsize=fontsize)
        cb.ax.tick_params(labelsize=fontsize)
        for tick in axbig.xaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        for tick in axbig.yaxis.get_major_ticks():
            tick.label.set_fontsize(fontsize)
        axbig.set_xlabel('meters', fontsize=fontsize)
        # plot the direction error
        if type(UAV) == list:
            print("Plotting "+str(len(UAV))+" UAVs")
            for ind, uav in enumerate(UAV):
                ax_v.plot([i for i in range(len(uav.planner.error))],
                          [i[0] for i in uav.planner.error],
                          color=self.colors[ind % len(self.colors)])
                ax_y.plot([i for i in range(len(uav.planner.error))],
                          [np.rad2deg(i[1]) for i in uav.planner.error],
                          color=self.colors[ind % len(self.colors)])
                axbig.scatter(x=x[uav.minX:uav.maxX, uav.minY:uav.maxY].flatten(),
                              y=y[uav.minX:uav.maxX, uav.minY:uav.maxY].flatten(),
                              c=uav.planner.hist[0][self._map_sel][uav.minX:uav.maxX,
                                                                   uav.minY:uav.maxY].flatten(),
                              s=15, cmap='gnuplot2')
                for coord, turbine in self.WF.turbines:
                    turbidx = self._findGPSindex([coord.x1, coord.x2], False)
                    if uav.minX <= turbidx[0] < uav.maxX and \
                            uav.minY <= turbidx[1] <= uav.maxY:
                        # x_0 = coord.x1 + np.sin(self.params.yBar) * turbine.rotor_radius
                        # x_1 = coord.x1 - np.sin(self.params.yBar) * turbine.rotor_radius
                        # y_0 = coord.x2 - np.cos(self.params.yBar) * turbine.rotor_radius
                        # y_1 = coord.x2 + np.cos(self.params.yBar) * turbine.rotor_radius
                        x_0 = coord.x1 + np.sin(uav.planner.hist[0][6]) * turbine.rotor_radius
                        x_1 = coord.x1 - np.sin(uav.planner.hist[0][6]) * turbine.rotor_radius
                        y_0 = coord.x2 - np.cos(uav.planner.hist[0][6]) * turbine.rotor_radius
                        y_1 = coord.x2 + np.cos(uav.planner.hist[0][6]) * turbine.rotor_radius
                        # axbig.plot([x_0, x_1], [y_0, y_1], linewidth=1, color='black')
                        axbig.plot([x_0, x_1], [y_0, y_1], linewidth=1, color=self.colors[ind % len(self.colors)])

        else:
            print("Plotting single UAV")
            ax_v.plot([i for i in range(len(UAV.planner.error))],
                      [i[0] for i in UAV.planner.error], color='black')
            ax_y.plot([i for i in range(len(UAV.planner.error))],
                      [np.rad2deg(i[1]) for i in UAV.planner.error], color='black')
            axbig.scatter(x=x[UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                          y=y[UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                          c=UAV.planner.hist[0][self._map_sel][UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                          s=15, cmap='gnuplot2')

            for coord, turbine in self.WF.turbines:
                x_0 = coord.x1 + np.sin(self.params.yBar) * turbine.rotor_radius
                x_1 = coord.x1 - np.sin(self.params.yBar) * turbine.rotor_radius
                y_0 = coord.x2 - np.cos(self.params.yBar) * turbine.rotor_radius
                y_1 = coord.x2 + np.cos(self.params.yBar) * turbine.rotor_radius
                x0_0 = coord.x1 + np.sin(UAV.planner.hist[0][6]) * turbine.rotor_radius
                x0_1 = coord.x1 - np.sin(UAV.planner.hist[0][6]) * turbine.rotor_radius
                y0_0 = coord.x2 - np.cos(UAV.planner.hist[0][6]) * turbine.rotor_radius
                y0_1 = coord.x2 + np.cos(UAV.planner.hist[0][6]) * turbine.rotor_radius
                axbig.plot([x_0, x_1], [y_0, y_1], linewidth=1, color='black')
                axbig.plot([x0_0, x0_1], [y0_0, y0_1], linewidth=1, color='black')

        # the plot's title holds a lot of info
        if type(UAV) == list:
            plt.suptitle("Wind Speed and Direction Estimates with a UAV swarm for Sensing\n" +
                         'Actual Speed: ' + str(self.params.vBar) +
                         ', Actual Direction: ' + str(np.rad2deg(self.params.dBar)) + '\N{DEGREE SIGN}' +
                         '\nPlanning Horizon: ' + str(UAV[0].plan_horizon) +
                         '          Moves until recalculation: ' + str(UAV[0].moves2recalc) +
                         '          Memory: ' + str(UAV[0].patrolMax) + ' nodes')
        else:
            plt.suptitle("Wind Speed and Direction Estimates with a UAV swarm for Sensing\n" +
                         'Actual Speed: ' + str(self.params.vBar) +
                         ', Actual Direction: ' + str(np.rad2deg(self.params.dBar)) + '\N{DEGREE SIGN}' +
                         '\nPlanning Horizon: ' + str(UAV.plan_horizon) +
                         '          Moves until recalculation: ' + str(UAV.moves2recalc) +
                         '          Memory: ' + str(UAV.patrolMax) + ' nodes')
        # set plot and axis titles for the error plots
        ax_y.set_ylabel('yaw error (\N{DEGREE SIGN})')
        ax_y.set_xlabel('# of recalculations')
        ax_v.set_ylabel('speed error (m/s)')
        # adjust how the plots fill the figure
        plt.subplots_adjust(left=0.05,
                            bottom=0.15,
                            right=0.95,
                            top=0.83,
                            wspace=0.27,
                            hspace=0.19)
        # slider axis
        sld_ax = plt.axes((0.2, 0.02, 0.56, 0.02))
        # create the slider
        LEN = 1000000000000000
        if type(UAV) == list:
            for uav in UAV:
                LEN = np.amin([len(uav.planner.hist), LEN])
        else:
            LEN = len(UAV.planner.hist)
        sld = Slider(sld_ax,
                     'moves',
                     0, LEN - 1, valinit=0)
        # set the initial slider text
        sld.valtext.set_text('move 0')
        # button axis
        btn_ax = plt.axes([0.85, 0.925, 0.125, 0.05])
        # create a button which toggles between score map and wave map
        btn = Button(btn_ax, 'Show Wave')
        # button axis
        btn2_ax = plt.axes([0.85, 0.01, 0.125, 0.05])
        # create a button which plays/pauses the animation
        if filename is None:
            btn2 = Button(btn2_ax, 'Play')
        else:
            btn2 = Button(btn2_ax, 'Pause')

        # function for the map button
        def map_btn(event):
            if self._map_sel == 0:  # currently showing the score map
                self._map_sel = 3  # change to showing the wave map
                btn.label.set_text('Show Score')  # update button text
            else:
                self._map_sel = 0  # otherwise change to showing the score map
                btn.label.set_text('Show Wave')  # update button text
            update_plot(0)  # and update the plot
        # set the above function for the map button
        btn.on_clicked(map_btn)

        # function for the play button
        def play_btn(event):
            if self._play:  # currently playing
                self._play = False  # pause it
                btn2.label.set_text("Play")  # update button text
            else:
                self._play = True  # otherwise, play the animation
                btn2.label.set_text("Pause")  # update button text
        # set the above function for the play button
        btn2.on_clicked(play_btn)

        # function which updates the plot
        def update_plot(val):
            # discretize the slider to integer values
            idx = int(round(sld.val))
            # set the slider text
            sld.valtext.set_text('move ' + '{}'.format(idx))
            # clear the plots
            axbig.clear()
            ax_v.clear()
            ax_y.clear()
            # check map selection
            if self._map_sel == 0:  # it's the score map
                btn.label.set_text('Show Wave')  # set the button text
            else:  # it's the wave map
                btn.label.set_text('Show Score')  # set the button text\
            # set the subplot labels
            ax_y.set_ylabel('direction error (\N{DEGREE SIGN})')
            ax_y.set_xlabel('# of recalculations')
            ax_v.set_ylabel('speed error (m/s)')
            # plot the map
            if type(UAV) == list:
                for ids, uv in enumerate(UAV):
                    axbig.scatter(x=x[uv.minX:uv.maxX, uv.minY:uv.maxY].flatten(),
                                  y=y[uv.minX:uv.maxX, uv.minY:uv.maxY].flatten(),
                                  c=uv.planner.hist[idx][self._map_sel][uv.minX:uv.maxX,
                                                                        uv.minY:uv.maxY].flatten(),
                                  s=15, cmap='gnuplot2')
                    # plot the plan
                    axbig.plot([i[0] for i in uv.planner.hist[idx][2]],
                               [i[1] for i in uv.planner.hist[idx][2]], linewidth=1.0,
                               color=self.colors[ids % len(self.colors)])
                    # plot the actual path
                    axbig.plot([i[0] for i in uv.planner.hist[idx][1]],
                               [i[1] for i in uv.planner.hist[idx][1]], linewidth=1.0,
                               color=self.colors[ids % len(self.colors)])
                    print(ids, np.rad2deg(uv.planner.hist[idx][6]))
                    for coor, turb in self.WF.turbines:
                        turbix = self._findGPSindex([coor.x1, coor.x2], False)
                        if uv.minX <= turbix[0] < uv.maxX and \
                                uv.minY <= turbix[1] <= uv.maxY:
                            # x1_0 = coor.x1 + np.sin(self.params.yBar) * turb.rotor_radius
                            # x1_1 = coor.x1 - np.sin(self.params.yBar) * turb.rotor_radius
                            # y1_0 = coor.x2 - np.cos(self.params.yBar) * turb.rotor_radius
                            # y1_1 = coor.x2 + np.cos(self.params.yBar) * turb.rotor_radius
                            x2_0 = coor.x1 + np.sin(uv.planner.hist[idx][6]) * turb.rotor_radius
                            x2_1 = coor.x1 - np.sin(uv.planner.hist[idx][6]) * turb.rotor_radius
                            y2_0 = coor.x2 - np.cos(uv.planner.hist[idx][6]) * turb.rotor_radius
                            y2_1 = coor.x2 + np.cos(uv.planner.hist[idx][6]) * turb.rotor_radius

                            # axbig.plot([x1_0, x1_1], [y1_0, y1_1], linewidth=1, color='black')
                            axbig.plot([x2_0, x2_1], [y2_0, y2_1], linewidth=1,
                                       color=self.colors[ids % len(self.colors)])

                    for tic in axbig.xaxis.get_major_ticks():
                        tic.label.set_fontsize(fontsize)
                    for tic in axbig.yaxis.get_major_ticks():
                        tic.label.set_fontsize(fontsize)
                    # for shorthand, UAV's location on the map
                    try:
                        _x = uv.planner.hist[idx][1][len(uv.planner.hist[idx][1])-1][0]
                        _y = uv.planner.hist[idx][1][len(uv.planner.hist[idx][1])-1][1]
                        # plot the UAV, a big red circle
                        axbig.plot(_x, _y, marker='o', markersize=10, color=self.colors[ids % len(self.colors)])
                    except:
                        print("couldn't draw UAV")
                        pass

                    # plot the velocity error
                    ax_v.plot([i for i in range(len(uv.planner.error))],
                              [i[0] for i in uv.planner.error],
                              color=self.colors[ids % len(self.colors)])

                    # plot the direction error
                    ax_y.plot([i for i in range(len(uv.planner.error))],
                              [np.rad2deg(i[1]) for i in uv.planner.error],
                              color=self.colors[ids % len(self.colors)])

            else:
                # print("single UAV")
                axbig.scatter(x=x[UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                              y=y[UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                              c=UAV.planner.hist[idx][self._map_sel][UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                              s=15,
                              cmap='gnuplot2')
                # plot the plan
                axbig.plot([i[0] for i in UAV.planner.hist[idx][2]],
                           [i[1] for i in UAV.planner.hist[idx][2]],
                           linewidth=1.0, color='lime')
                # plot the actual path
                axbig.plot([i[0] for i in UAV.planner.hist[idx][1]],
                           [i[1] for i in UAV.planner.hist[idx][1]],
                           linewidth=1.0, color='red')
                for tic in axbig.xaxis.get_major_ticks():
                    tic.label.set_fontsize(fontsize)
                for tic in axbig.yaxis.get_major_ticks():
                    tic.label.set_fontsize(fontsize)
                # for shorthand, UAV's location on the map
                try:
                    _x = UAV.planner.hist[idx][1][len(UAV.planner.hist[idx][1]) - 1][0]
                    _y = UAV.planner.hist[idx][1][len(UAV.planner.hist[idx][1]) - 1][1]
                    # plot the UAV, a big red circle
                    axbig.plot(_x, _y, marker='o', markersize=10, color='red')
                except:
                    print("couldn't draw UAV")
                    pass

                # plot the velocity error
                ax_v.plot([i for i in range(len(UAV.planner.error))],
                          [i[0] for i in UAV.planner.error], color='black')
                # plot the direction error
                ax_y.plot([i for i in range(len(UAV.planner.error))],
                          [np.rad2deg(i[1]) for i in UAV.planner.error], color='black')
                # plot the turbines
                for coor, turb in self.WF.turbines:
                    x2_0 = coor.x1 + np.sin(self.params.yBar) * turb.rotor_radius
                    x2_1 = coor.x1 - np.sin(self.params.yBar) * turb.rotor_radius
                    y2_0 = coor.x2 - np.cos(self.params.yBar) * turb.rotor_radius
                    y2_1 = coor.x2 + np.cos(self.params.yBar) * turb.rotor_radius
                    x3_0 = coor.x1 + np.sin(UAV.planner.hist[idx][6]) * turb.rotor_radius
                    x3_1 = coor.x1 - np.sin(UAV.planner.hist[idx][6]) * turb.rotor_radius
                    y3_0 = coor.x2 - np.cos(UAV.planner.hist[idx][6]) * turb.rotor_radius
                    y3_1 = coor.x2 + np.cos(UAV.planner.hist[idx][6]) * turb.rotor_radius
                    axbig.plot([x2_0, x2_1], [y2_0, y2_1], linewidth=1, color='black')
                    axbig.plot([x3_0, x3_1], [y3_0, y3_1], linewidth=1, color='red')

            # show the plot
            plt.draw()
        # set the update function for what happens when the slider value changes
        sld.on_changed(update_plot)
        # adjust how the plots fill the figure
        plt.subplots_adjust(left=0.1,
                            bottom=0.14,
                            right=1.0,
                            top=0.84,
                            wspace=0.2,
                            hspace=0.24)

        # function to animate the plot
        def animate(frame, *fargs):
            # check if it is playing
            if self._play:
                # if it's not at the end, increment the slider value
                if sld.val < sld.valmax-1:
                    temp = sld.val
                    sld.set_val(temp + 1)
                else:
                    # if it's at the end, set it to the beginning
                    sld.set_val(sld.valmin)

        # set the animate function to the FuncAnimation function for animation
        an = anim.FuncAnimation(f, animate, interval=100, frames=LEN)
        # render to video. to make it play faster, increase fps
        if filename is not None:
            an.save(filename+'.mp4', fps=15, dpi=300)
        # show the plot
        if plot:
            plt.show()
