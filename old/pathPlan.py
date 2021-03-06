## \file pathPlan.py
# This file contains functions which handle UAV path planning using Floris data.
#
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from floris.coordinate import Coordinate
from visualization_manager_DJ import VisualizationManager
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
    ## @var hist
    # keeps track of everything that happens so it can be plotted afterward.
    #
    # This is how history is populated:
    # self.hist.append(deepcopy([self.score_map,
    #                                UAV.GPSpath,
    #                                UAV.GPSplan,
    #                                UAV.wave_map,
    #                               self.error[len(self.error)-1],
    #                               self.params.v0,
    #                               self.params.d0]))


    ## Class constructor
    #   @param vman a VisualManager object
    def __init__(self, vman):
        ## @var vman
        #   A VisualizationManager object
        self.vman = deepcopy(vman)
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            turbine.yaw_angle = self.vman.flowfield.wind_direction
        self.vman.flowfield.calculate_wake()
        ## @var params
        #   A local copy of vman.params
        self.params = deepcopy(vman.params)
        ## @var plane
        # the plane where the computations and visualizations are taking place
        self.plane = int(self.vman.flowfield.grid_resolution.z * self.params.percent_height)
        ## @var WF
        #  a local copy of vman.WF
        self.WF = deepcopy(vman.WF)
        # time to calculate Xbar
        # set the wind_speed in the JSON object to vBar
        vman.WF['farm']['properties']['wind_speed'] = self.vman.params.vBar
        # create a new VisualizationManager object with vBar as the velocity
        vman_bar = VisualizationManager(vman.WF, vman.grid_res)
        # set the direction to dBar
        vman_bar.flowfield.wind_direction = vman_bar.params.dBar
        for coord, turbine in vman_bar.flowfield.turbine_map.items():
            turbine.yaw_angle = vman_bar.flowfield.wind_direction
        # recalculate the wake with the new direction
        vman_bar.flowfield.calculate_wake()
        ## @var Xbar
        #   this is the 'actual' u_field
        self.Xbar = deepcopy(vman_bar.flowfield.u_field[:, :, self.plane])
        ## @var X_map
        #   A map of all possible transitions on the map and their respective scores
        self.X_map = [[_X_map_node() for x in range(self.vman.grid_res[1])]
                 for y in range(self.vman.grid_res[0])]
        ## @var X
        #   An X grid mesh of the map's 'GPS' coordinates
        self.X = deepcopy(self.vman.flowfield.x[:, :, self.plane])
        ## @var Y
        #   A Y grid mesh of the map's 'GPS' coordinates
        self.Y = deepcopy(self.vman.flowfield.y[:, :, self.plane])
        ## the local copy of the sensitivity matrix
        self.sens_mat = self.vman.calcSensitivityMatrix(self.params.v0,self.params.d0)
        ## @var error
        # holds a record of the direction and velocity errors for plotting
        self.error = list()
        # append the initial error
        self.error.append([self.params.vBar-self.params.v0,
                           self.params.dBar-self.params.d0])
        # delete the old vman objects
        del vman, vman_bar
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
        self.VTKrange = {0, 0}
        self.VTKpath = None
        self.VTKfilename = None
        self.VTKincrement = 5
        self.VTKlist = list()
        self.VTK_FULL = False

    ## a method to adjust idx based on a cardinal direction
    def _shiftVals(self, val):
        value = {
            d.E: [1,0],   # one right, zero up
            d.NE: [1,1],  # one right, one up
            d.N: [0,1],   # zero right, one up
            d.NW: [-1,1], # one left, one up
            d.W: [-1,0],  # one left, zero up
            d.SW: [-1,-1],# one left, one down
            d.S: [0,-1],  # zero left, one down
            d.SE: [1,-1]  # one right, one down
        }
        return value.get(val)

    ## calculates the euclidean distance between two points
    def euclidDist(self, P1, P2):
        # sqrt((Y2-Y1)^2+(X2-X1)^2)
        return math.sqrt(((P2[0]-P1[0])**2)+((P2[1]-P1[1])**2))

    ## a method to update the transition/cost map
    def _init_cost_map(self):
        # iterate over the entire range of the discretized map
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                # record the index of the node
                # so that it can easily be retrieved if non-indexed
                # search methods are used
                self.X_map[i][j].idx = [i,j]
                # record the 'GPS' location of the node
                self.X_map[i][j].GPS = [self.X[i][j],self.Y[i][j]]
                # iterate through the 8 possible transitions at each node
                for k in range(8):
                    # calculate the cost of moving through the wind at the transition heading
                    self.X_map[i][j].Xitions[k].thetaCost = k*np.deg2rad(45)-self.params.d0
                    # adjust the cost to a positive value 0-PI
                    if self.X_map[i][j].Xitions[k].thetaCost > PI:
                        self.X_map[i][j].Xitions[k].thetaCost = \
                            abs(k*np.deg2rad(45)-self.params.d0-2*PI)
                    # test the map indices after the transition
                    [I,J] = self._shiftVals(k)
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
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                # record the index of the node
                # so that it can easily be retrieved if non-indexed
                # search methods are used
                self.X_map[i][j].idx = [i,j]
                # record the 'GPS' location of the node
                self.X_map[i][j].GPS = [self.X[i][j],self.Y[i][j]]
                # iterate through the 8 possible transitions at each node
                for k in range(8):
                    # calculate the cost of moving through the wind at the transition heading
                    self.X_map[i][j].Xitions[k].thetaCost = k*np.deg2rad(45)-UAV.d0
                    # adjust the cost to a positive value 0-PI
                    if self.X_map[i][j].Xitions[k].thetaCost > PI:
                        self.X_map[i][j].Xitions[k].thetaCost = \
                            abs(k*np.deg2rad(45)-UAV.d0-2*PI)
                    # test the map indices after the transition
                    [I,J] = self._shiftVals(k)
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
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                # create a no-fly zone around each turbine
                for coord, turbine in self.vman.flowfield.turbine_map.items():
                    # compute if the euclidean distance is less
                    # than the minimum distance threshold minDist2turbine
                    if self.euclidDist([self.X[i][j],
                                         self.Y[i][j]],
                                        [coord.x, coord.y]) < \
                            self.minDist2turbine:
                        # if it is, make it untraversable
                        self.score_map[i][j] = None
        # calculate a wave map from this score map

    ## calculates a wave map using the highest score on the score map
    def _calcWaveMap(self, UAV, plan_mask=None):
        # set the max unobtainably low (since our range is 0-1)
        if(plan_mask==None):
            plan_mask=UAV.plan_mask
        max = -100
        idx = [0,0]
        # set the wave_map to all ones
        UAV.wave_map = [[1 for i in range(self.vman.grid_res[1])] \
                          for j in range(self.vman.grid_res[0])]
        # iterate through the entire map
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                # if we aren't using a UAV for this calculation
                # (such as for the initial wave map calculation)
                # the program will throw an exception
                try:
                    # we test the score using the UAV's path mask
                    if self.score_map[i][j]*plan_mask[i][j]*UAV.path_mask[i][j]>max:
                        max=self.score_map[i][j]
                        idx = [i,j]

                except:
                    if self.score_map[i][j]>max:
                        max=self.score_map[i][j]
                        idx = [i,j]
        # just for shorthand
        x = idx[0]
        y = idx[1]
        # iterate through the map again
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                # set the value of the wave map to the euclidean
                # distance from the maximum score
                try:
                    UAV.wave_map[i][j] = 1/(self.euclidDist([x,y],[i,j])*self.euclidDist([x,y],[i,j]))
                except:
                    UAV.wave_map[i][j] = 1
        # normalize the wave values to 0-1,
        # but with increasing values as they approach the maximum
        #UAV.wave_map = 1-UAV.wave_map/np.amax(UAV.wave_map)
        UAV.max_wave_idx = idx
        #print(idx)

    ## Simply calls the plt.show() method to show any figures
    #   which have been generated.
    def show(self):
        plt.show()

    ## Generates a plot of the current wave map
    #   @param ax can be passed to the function if the plot
    #           is to be in a subplot
    def plotWaveMap(self, UAV, ax=None):
        # a scatter plot of the wave map
        if(ax):
            # using a subplot axis that was generated elsewhere
            axx = ax.scatter(x=self.X, y=self.Y, c=UAV.wave_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx, ax=ax)
        else:
            # standalone plot
            axx = plt.scatter(x=self.X, y=self.Y, c=UAV.wave_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx)

        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle-self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle-self.vman.flowfield.wind_direction, coord.as_tuple())
            if ax:
                ax.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')
            else:
                plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')

    ## Generates a plot of the current score map
    #   @param ax can be passed to the function if the plot
    #           is to be in a subplot
    def plotScoreMap(self, ax=None):
        # a scatter plot of the score map
        if ax:
            # using a subplot axis that was generated elsewhere
            axx = ax.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx, ax=ax)
        else:
            # standalone plot
            axx = plt.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx)

        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle-self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle-self.vman.flowfield.wind_direction, coord.as_tuple())
            if ax:
                ax.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')
            else:
                plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')

    ## Generates a plot of a UAV's score mask (planned or actual)
    #   @param mask the mask to be plotted (planned or actual)
    #   @param ax Can be passed to the function if the plot
    #           is to be in a subplot
    def plotScoreMask(self, mask, ax=None):
        # a scatter plot of the UAV's score mask
        if ax:
            # using a subplot axis that was generated elsewhere
            axx = ax.scatter(x=self.X, y=self.Y, c=mask)
            plt.colorbar(axx, ax=ax)
        else:
            # standalone plot
            axx = plt.scatter(x=self.X, y=self.Y, c=mask)
            plt.colorbar(axx)

        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle, coord.as_tuple())
            if ax:
                ax.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=0.2, color='lime')
            else:
                plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=0.2, color='lime')

    ## Generates a plot of a given UAV's path over the score map
    #   @param UAV the UAV whose path is to be plotted
    #   @param ax Can be passed to the function if the plot
    #           is to be in a subplot
    def plotScoreMapUAV(self, UAV, ax=None):
        # a score map with the UAV's planned path and actual path overlayed
        if ax:
            # using an axis that was generated elsewhere
            axx = ax.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx, ax=ax)
        else:
            # standalone plot
            axx = plt.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx)
        plt.axes().set_aspect('equal')
        # the title
        plt.suptitle("Greedy UAV path\n"+str(len(UAV.GPSplan)-1)+" moves")
        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle, coord.as_tuple())
            if ax:
                ax.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
            else:
                plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
        # plot the UAV's path
        if ax:
            # planned path
            ax.plot([i[0] for i in UAV.GPSplan],
                     [i[1] for i in UAV.GPSplan], color='lime')
            # actual path
            ax.plot([i[0] for i in UAV.GPSpath],
                     [i[1] for i in UAV.GPSpath], color='red')
        else:
            # planned path
            plt.plot([i[0] for i in UAV.GPSplan],
                    [i[1] for i in UAV.GPSplan], color='lime')
            # actual path
            plt.plot([i[0] for i in UAV.GPSpath],
                    [i[1] for i in UAV.GPSpath], color='red')
        # size the plot to fill the window
        if not ax:
            plt.subplots_adjust(left=0.07,
                                bottom=0.08,
                                right=1.0,
                                top=0.91,
                                wspace=0.27,
                                hspace=0.19)

    ## A method to convert an index to a GPS point
    def _findGPSindex(self, GPS):
        # Figure out how 'wide' each range is
        XleftSpan = np.amax(self.X) - np.amin(self.X)
        XrightSpan = self.vman.grid_res[0]
        YleftSpan = np.amax(self.Y) - np.amin(self.Y)
        YrightSpan = self.vman.grid_res[1]

        # Convert the left range into a 0-1 range (float)
        idx = int((GPS[0] -  np.amin(self.X))* XrightSpan / float(XleftSpan))
        idy = int((GPS[1] -  np.amin(self.Y))* YrightSpan / float(YleftSpan))
        # Is it a valid GPS point?
        if math.isnan(self.score_map[idx][idy]) or idx<0 or idy<0:
            # if not, let the user know
            print("Error: GPS point is out of bounds: "+str(GPS))
            return [None, None]
        else:
            # otherwise, send back the coordinates
            return [idx, idy]

    ## a method to convert a GPS point to an index
    def _getGPSfromIDX(self, IDX):
        # we already have indexed lists of these values
        X = self.X[IDX[0]][IDX[1]]
        Y = self.Y[IDX[0]][IDX[1]]
        return [X,Y]

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
        if (UAV.idx):
            # calculate the current score map
            self.calcScoreMap()
            self._update_cost_map(UAV)
            # let the iteration begin
            for i in range(plan_horizon):
                # just for shorthand notation
                x=UAV.IDXplan[i][0]
                y=UAV.IDXplan[i][1]

                # if the indices are invalid, an exception will be thrown
                try:
                    # start by updating the UAV's score mask for the current position
                    UAV.update_mask(UAV.plan_mask, [x,y])
                    #print("1")
                    # calculate the new wave map with the updated mask
                    self._calcWaveMap(UAV)
                    # print("wave map calculated")
                    # calculate the next step
                    max, dir = self._greedyStep(UAV)
                    UAV.planner.steps_planned = UAV.planner.steps_planned + 1
                    #print(str(max) +','+str(dir))
                    # print("step planned")
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
        if x>=self.vman.grid_res[0] or \
            y>=self.vman.grid_res[1] or \
            x<0 or y<0:
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
            #print("1")
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
                #print("1")
                # pull the UAV's mask value for the transition node
                mask = UAV.plan_mask[x+I][y+J]
                #print("2")
                # pull the wave map value for the transition node
                wave = UAV.wave_map[x+I][y+J]
                #print("3")
                # calculate the cost of changing the UAV's heading
                headCost = abs(i*np.deg2rad(45)- \
                                  (UAV.plan_heading[len(UAV.plan_heading)-1])*
                                     np.deg2rad(45))
                #print("4")
                # normalize to a value between -PI and PI
                if headCost > PI:
                    headCost = abs(headCost - 2 * PI)
                # calculate the 2nd derivative of the score wrt the transition
                d2S = Xition[i].dSscore-UAV.dSplan[len(UAV.dSplan)-1]
                #print("5")
                # calculate the reward of the move
                reward = (UAV.scoreWt*score + \
                       UAV.dSwt*Xition[i].dSscore + \
                       UAV.d2Swt*d2S - \
                       UAV.windWt*Xition[i].thetaCost/PI - \
                       UAV.headWt*headCost/PI + \
                       UAV.waveWt*wave + \
                       UAV.timeWt*len(UAV.IDXplan)) * \
                       mask
                #print("6")
                # we are trying to maximize the reward
                if reward > max and \
                    x+I<self.vman.grid_res[0] and \
                    y+J<self.vman.grid_res[1] and \
                    x+I>=0 and y+J>=0:

                    # if the transition is valid and has the highest reward
                    # then track the reward score and the direction of the transition
                    max = reward
                    #print("max: "+ str(max))
                    dir = i
                #print("7")
                #print("index" + str(i))
            # if the transition is invalid, skip it
            except:
                pass
        # return the decision and the resulting score
        #print("direction:" + str(dir))
        return max, dir

    # updateEstimates
    ## A function which updates the direction and speed estimates
    # based on the UAV's path mask over Xbar
    #
    # @param UAV The UAV whose path is being used to update the estimates
    def updateEstimates(self, UAV):
        print("Updating Estimates...")
        # current field estimate
        Xk = deepcopy(self.vman.flowfield.u_field[:, :, self.plane])
        # current velocity and direction estimates
        vd_k = [[deepcopy(UAV.v0)],
                [deepcopy(UAV.d0)]]
        # current sensitivity matrix
        sens_mat = np.column_stack([self.sens_mat[0].flatten(), self.sens_mat[1].flatten()])
        # apply the pseudoinverse
        sens_mat_pinv = np.linalg.pinv(sens_mat)
        # field actual
        if self.Xbar.ndim == 2:
            Xbar = deepcopy(self.Xbar)
        else:
            Xbar = deepcopy(self.Xbar[self.steps_taken].flatten())
        # if the node is not masked by the UAV's path, use the estimate
        self.mask_size_history.append(1-float(np.sum(UAV.path_mask))/ \
                    (self.vman.grid_res[0]*self.vman.grid_res[1]))
        IDM = list()
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                if UAV.path_mask[i][j] == 1:
                    Xbar[i][j] = Xk[i][j]
                else:
                    IDM.append(self.score_map[i][j])
        #self.IDMhistory_max.append(np.amax(IDM))
        #self.IDMhistory_min.append(np.amin(IDM))
        #self.IDMhistory_avg.append(np.average(IDM))
        # convert Xk and Xbar to column vectors
        Xbar = np.vstack(Xbar.flatten())
        Xk = np.vstack(Xk.flatten())
        # calculate the next estimate: Vdk+1 = Vdk + sens_mat_pinv*(Xbar-Xk)
        vd_kp1 = vd_k + np.matmul(sens_mat_pinv, Xbar - Xk)
        # update the estimates
        UAV.v0 = vd_kp1[0][0]
        UAV.d0 = vd_kp1[1][0]
        # update the wind speed estimate
        self.WF['farm']['properties']['wind_speed'] = UAV.v0
        # calculate the FLORIS model with the new speed
        vmanTemp = VisualizationManager(self.WF, self.vman.grid_res)
        # update the direction estimate
        vmanTemp.flowfield.wind_direction = UAV.d0
        # recalculate the FLORIS wake model with the new direction estimate
        turbines = [turbine for _, turbine in vmanTemp.flowfield.turbine_map.items()]
        for k, turbine in enumerate(turbines):
            turbine.yaw_angle = UAV.d0
        vmanTemp.flowfield.calculate_wake()
        # update the sensitivity matrix
        self.sens_mat = deepcopy(vmanTemp.calcSensitivityMatrix(UAV.v0,UAV.d0))
        # update the score map
        self.calcScoreMap()
        # output the error to the terminal

        # update the local VisualManager object
        self.vman = deepcopy(vmanTemp)
        # keep track of the error signals for plotting
        self.error.append([self.params.vBar - UAV.v0,
                           self.params.dBar - UAV.d0])
        print("Error: [v,theta]=" + str(self.error[len(self.error) - 1]))

    def updateEstimatesRM(self, UAV):
        print("Updating Estimates...")
        print("Prior estimate: " + str([UAV.v0, UAV.d0]))
        # current field estimate
        temp_Xk = deepcopy(self.vman.flowfield.u_field[:, :, self.plane]).flatten()
        # current velocity and direction estimates
        vd_k = [[deepcopy(UAV.v0)],
                [deepcopy(UAV.d0)]]
        # current sensitivity matrix
        temp_sens_mat0 = self.sens_mat[0].flatten()
        temp_sens_mat1 = self.sens_mat[1].flatten()

        # field actual
        if self.Xbar.ndim == 2:
            temp_Xbar = deepcopy(self.Xbar).flatten()
        else:
            temp_Xbar = deepcopy(self.Xbar[self.steps_taken]).flatten()
        temp_mask = deepcopy(np.array(UAV.path_mask).flatten())
        # if the node is not masked by the UAV's path, use the estimate
        self.mask_size_history.append(1-float(np.sum(UAV.path_mask))/ \
                    (self.vman.grid_res[0]*self.vman.grid_res[1]))

        sens_mat0 = list()
        sens_mat1 = list()
        Xbar = list()
        Xk = list()
        IDM = list()
        temp_score_map = deepcopy(self.score_map).flatten()
        for i in range(len(temp_mask)):
                if temp_mask[i] != 1:
                    Xbar.append(temp_Xbar[i])
                    Xk.append(temp_Xk[i])
                    sens_mat0.append(temp_sens_mat0[i])
                    sens_mat1.append(temp_sens_mat1[i])
                    IDM.append(temp_score_map[i])


        # convert Xk and Xbar to column vectors
        Xbar = np.vstack(np.array(Xbar))
        Xk = np.vstack(np.array(Xk))
        sens_mat0 = np.vstack(np.array(sens_mat0))
        sens_mat1 = np.vstack(np.array(sens_mat1))
        #print("Xbar shape")
        #print(Xbar.shape)
        #print("Xk shape")
        #print(Xk.shape)
        #print("sens_mat0 shape")
        #print(sens_mat0.shape)
        #print("sens_mat1 shape")
        #print(sens_mat1.shape)
        # apply the pseudoinverse
        sens_mat_pinv = np.linalg.pinv(np.column_stack([sens_mat0, sens_mat1]))
        #print("sens_mat_pinv shape")
        #print(sens_mat_pinv.shape)
        # calculate the next estimate: Vdk+1 = Vdk + sens_mat_pinv*(Xbar-Xk)
        vd_kp1 = vd_k + np.matmul(sens_mat_pinv, Xbar - Xk)
        # update the estimates
        UAV.v0 = vd_kp1[0][0]
        UAV.d0 = vd_kp1[1][0]
        # update the wind speed estimate
        self.WF['farm']['properties']['wind_speed'] = UAV.v0
        # calculate the FLORIS model with the new speed
        vmanTemp = VisualizationManager(self.WF, self.vman.grid_res)
        # update the direction estimate
        vmanTemp.flowfield.wind_direction = UAV.d0
        # recalculate the FLORIS wake model with the new direction estimate
        turbines = [turbine for _, turbine in vmanTemp.flowfield.turbine_map.items()]
        for k, turbine in enumerate(turbines):
            turbine.yaw_angle = UAV.d0
        vmanTemp.flowfield.calculate_wake()
        # update the sensitivity matrix
        self.sens_mat = deepcopy(vmanTemp.calcSensitivityMatrix(UAV.v0,UAV.d0))
        # update the score map
        self.calcScoreMap()
        # output the error to the terminal

        # update the local VisualManager object
        self.vman = deepcopy(vmanTemp)
        # keep track of the error signals for plotting
        error = [self.params.vBar - UAV.v0,
                           self.params.dBar - UAV.d0]
        print("Actual: "+str([self.params.vBar, self.params.dBar]))
        print("New estimate: " + str([UAV.v0, UAV.d0]))
        print("Error: [v,theta]=" + str(error))


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
        self.colors = ['b', 'g', 'r', 'c', 'm', 'y', 'lime', 'orange','purple']
        # variable which tells us if we are looking at the Wave or the Score map
        self._map_sel = 0
        # boolean which tells us if the animation is playing
        if filename==None:
            self._play = False
        else:
            self._play = True
        fontsize = 14
        # create a figure
        f = plt.figure(figsize=(10,9))
        # separate it into a grid
        gs = f.add_gridspec(2,3)
        # put the error plots in the leftmost figures
        ax_v = f.add_subplot(gs[0,0], title='speed error')
        ax_d = f.add_subplot(gs[1,0])
        # the path plot takes up the right 2/3 of the figure
        axbig = f.add_subplot(gs[0:,1:])
        axbig.set_aspect('equal')
        # plot the velocity error
        ax_v.plot([i for i in range(len(self.error))], [i[0] for i in self.error])
        # format x,y and c for scatterplotting
        x = deepcopy(self.X)
        y = deepcopy(self.Y)
        if type(UAV)==list:
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
        if type(UAV)==list:
            for ind, uav in enumerate(UAV):
                ax_v.plot([i for i in range(len(uav.planner.error))], [i[0] for i in uav.planner.error],color=self.colors[ind % len(self.colors)])
                ax_d.plot([i for i in range(len(uav.planner.error))], [np.rad2deg(i[1]) for i in uav.planner.error],color=self.colors[ind % len(self.colors)])
                axbig.scatter(x=x[uav.minX:uav.maxX,uav.minY:uav.maxY].flatten(), y=y[uav.minX:uav.maxX,uav.minY:uav.maxY].flatten(),
                              c=uav.planner.hist[0][self._map_sel][uav.minX:uav.maxX,uav.minY:uav.maxY].flatten(), s=15,
                              cmap='gnuplot2')
        else:
            ax_v.plot([i for i in range(len(UAV.planner.error))],
                      [i[0] for i in UAV.planner.error], color='black')
            ax_d.plot([i for i in range(len(UAV.planner.error))],
                      [np.rad2deg(i[1]) for i in UAV.planner.error],color='black')
            axbig.scatter(x=x[UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                          y=y[UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(),
                          c=UAV.planner.hist[0][self._map_sel][UAV.minX:UAV.maxX, UAV.minY:UAV.maxY].flatten(), s=15,
                          cmap='gnuplot2')
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            axbig.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
        # the plot's title holds a lot of info
        if type(UAV)==list:
            plt.suptitle("Wind Speed and Direction Estimates with a UAV swarm for Sensing\n" +
                         'Actual Speed: ' + str(self.params.vBar) +
                         ', Actual Direction: ' + str(np.rad2deg(self.params.dBar)) + '\N{DEGREE SIGN}' +
                         '\nPlanning Horizon: ' + str(UAV[0].plan_horizon) +
                         '          Moves until recalculation: ' + str(UAV[0].moves2recalc) +
                         '          Memory: ' + str(UAV[0].patrolMax) + ' nodes')
                                 # set plot and axis titles for the error plots
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
                LEN = np.amin([len(uav.planner.hist),LEN])
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
        if filename==None:
            btn2 = Button(btn2_ax, 'Play')
        else:
            btn2 = Button(btn2_ax, 'Pause')
        # function for the map button
        def map_btn(event):
            if self._map_sel == 0: # currently showing the score map
                self._map_sel = 3 # change to showing the wave map
                btn.label.set_text('Show Score') # update button text
            else:
                self._map_sel = 0 # otherwise change to showing the score map
                btn.label.set_text('Show Wave') # update button text
            update_plot(0) # and update the plot
        # set the above function for the map button
        btn.on_clicked(map_btn)
        # function for the play button
        def play_btn(event):
            if self._play: # currently playing
                self._play = False # pause it
                btn2.label.set_text("Play") # update button text
            else:
                self._play = True # otherwise, play the animation
                btn2.label.set_text("Pause") # update button text
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
            if self._map_sel == 0: # it's the score map
                btn.label.set_text('Show Wave') # set the button text
            else: # it's the wave map
                btn.label.set_text('Show Score') # set the button text\
            # set the subplot labels
            ax_d.set_ylabel('direction error (\N{DEGREE SIGN})')
            ax_d.set_xlabel('# of recalculations')
            ax_v.set_ylabel('speed error (m/s)')
            # plot the map
            if type(UAV)==list:
                for ind, uav in enumerate(UAV):
                    axbig.scatter(x=x[uav.minX:uav.maxX,uav.minY:uav.maxY].flatten(),
                                  y=y[uav.minX:uav.maxX,uav.minY:uav.maxY].flatten(),
                                  c = uav.planner.hist[idx][self._map_sel][uav.minX:uav.maxX,uav.minY:uav.maxY].flatten(), s = 15,
                                  cmap = 'gnuplot2')
                    # plot the plan
                    axbig.plot([i[0] for i in uav.planner.hist[idx][2]],
                               [i[1] for i in uav.planner.hist[idx][2]], linewidth=1.0, color=self.colors[ind % len(self.colors)])
                    # plot the actual path
                    axbig.plot([i[0] for i in uav.planner.hist[idx][1]],
                                     [i[1] for i in uav.planner.hist[idx][1]], linewidth=1.0, color = self.colors[ind % len(self.colors)])
                    for tick in axbig.xaxis.get_major_ticks():
                        tick.label.set_fontsize(fontsize)
                    for tick in axbig.yaxis.get_major_ticks():
                        tick.label.set_fontsize(fontsize)
                    # for shorthand, UAV's location on the map
                    try:
                        _x = uav.planner.hist[idx][1][len(uav.planner.hist[idx][1])-1][0]
                        _y = uav.planner.hist[idx][1][len(uav.planner.hist[idx][1])-1][1]
                        # plot the UAV, a big red circle
                        axbig.plot(_x,_y, marker='o', markersize=10, color=self.colors[ind % 8])
                    except:
                        print("couldn't draw UAV")
                        pass

                    # plot the velocity error
                    ax_v.plot([i for i in range(len(uav.planner.error))], [i[0] for i in uav.planner.error],color=self.colors[ind % 8])

                    # plot the direction error
                    ax_d.plot([i for i in range(len(uav.planner.error))], [np.rad2deg(i[1]) for i in uav.planner.error],color=self.colors[ind % 8])

            else:
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
                for tick in axbig.xaxis.get_major_ticks():
                    tick.label.set_fontsize(fontsize)
                for tick in axbig.yaxis.get_major_ticks():
                    tick.label.set_fontsize(fontsize)
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
                          [i[0] for i in UAV.planner.error],color='black')
                # add a vertical red line to show where we are in the iterations
                # if idx > uav.init_mask_size:
                #    ax_v.axvline((idx-uav.init_mask_size+1)/(uav.moves2recalc)+1, color='red')
                # else:
                #    ax_v.axvline(idx/uav.init_mask_size, color='red')
                # plot the direction error
                ax_d.plot([i for i in range(len(UAV.planner.error))],
                          [np.rad2deg(i[1]) for i in UAV.planner.error],color='black')
                # add a vertical red line to show where we are in the iterations
                # if idx > uav.init_mask_size:
                #    ax_d.axvline((idx-uav.init_mask_size+1) / (uav.moves2recalc)+1, color='red')
                # else:
                #    ax_d.axvline(idx/uav.init_mask_size, color='red')
            # plot the turbines
            for coord, turbine in self.vman.flowfield.turbine_map.items():
                a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
                b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
                a.rotate_z(turbine.yaw_angle-self.vman.flowfield.wind_direction, coord.as_tuple())
                b.rotate_z(turbine.yaw_angle-self.vman.flowfield.wind_direction, coord.as_tuple())
                axbig.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
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
        if type(UAV)==list:
            an = anim.FuncAnimation(f, animate, interval=100,frames=len(UAV[0].planner.hist))
        else:
            an = anim.FuncAnimation(f, animate, interval=100,frames=len(UAV.planner.hist))
        # render to video. to make it play faster, increase fps
        if filename:
            an.save(filename+'.mp4',fps=15,dpi=300)
        # show the plot
        if plot:
            plt.show()

    def learnPath(self, UAV):
        #if not plan_horizon:
        plan_horizon = self._init_planner(UAV)

        print("learn horizon: "+str(self.learn_path_length))
        # if the GPS point was valid, the indices should exist
        if (UAV.idx):
            # calculate the current score map
            switches=0
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
                #print("learning iteration " + str(j))
                self.head = 0
                self.action = 0
                reward = 0
                SWITCH = False

                for i in range(self.learn_path_length):
                    # start with a fresh Qmap
                    #print("iteration "+str(i))
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

                        if X+I<self.vman.grid_res[0] and Y+J<self.vman.grid_res[1] and \
                            X+I>=0 and Y+J>=0 and not math.isnan(self.score_map[X+I][Y+J]):
                                maxQ = np.amax(Qmap[X+I][Y+J])
                                Xition = self.X_map[X][Y].Xitions
                                # pull the score from the score map for the transition node
                                score = self.score_map[X + I][Y + J]
                                self.scoreTOT = self.scoreTOT + score
                                #print(self.scoreTOT)
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
                                # print("5")
                                # calculate the reward of the move
                                '''
                                temp_reward = (UAV.scoreWt * score - \
                                               UAV.windWt * Xition[self.action].thetaCost / PI - \
                                               UAV.headWt * headCost / PI + \
                                               UAV.waveWt * wave) -(mask)*UAV.maskSUB
                                               '''
                                temp_reward = self.scoreTOT
                                #print("temp_reward: "+str(temp_reward))
                                if not math.isnan(temp_reward):
                                    reward = temp_reward
                                    if X+I == UAV.max_wave_idx[0] and Y+J == UAV.max_wave_idx[1]:
                                        SWITCH = True
                                        #print(reward)
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
                                                    if(np.sqrt(boxi*boxi+boxj*boxj)<self.max_score_box):
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
                                            #plan_mask = deepcopy(UAV.path_mask)
                                            SWITCH = False
                        else:
                            Qmap[X][Y][self.action] = -100000000000

            print("learning complete.")
            print("plan horizon: " + str(plan_horizon))
            print(" Planning...")
            #print(Qmap)
            X = UAV.idx[0]
            Y = UAV.idx[1]
            switch=0
            prevMaxQ = np.amin(Qmap)
            prevCoords=[X,Y]
            if (len(UAV.GPSpath) == 0):
                plan_len = UAV.init_plan_horizon
            else:
                plan_len = UAV.plan_horizon
            for n in range(plan_len):
                step = False
                while not step and n<plan_len:
                    maxQ = np.amin(Qmap)
                    index = 0
                    for j in range(8):

                        if Qmap[X][Y][j]>maxQ:
                            [I, J] = self._shiftVals(j)
                            if X + I < self.vman.grid_res[0] and X + I >= 0 and \
                                    Y + J < self.vman.grid_res[1] and Y + J >= 0\
                                    and X+I != prevCoords[0] \
                                    and Y+J != prevCoords[1]:
                                        if self.euclidDist([X+I,Y+J],prevCoords)<np.sqrt(2):
                                            pass
                                        else:
                                            maxQ = Qmap[X][Y][j]
                                            Qmap[X][Y][j] = Qmap[X][Y][j] - 100
                                            #print(maxQ)
                                            index = j
                    if maxQ < prevMaxQ:
                        plan_len=plan_len-1
                        step = True
                    elif maxQ==prevMaxQ or maxQ==-100000000000:
                    #if maxQ == prevMaxQ:
                        break

                    else:
                        step = True
                if maxQ == prevMaxQ or maxQ==-100000000000:
                    break
                print([maxQ,prevMaxQ,[X,Y]])
                prevMaxQ = maxQ
                # update the UAV's heading
                UAV.plan_heading.append(index)
                # update idx based on the direction of travel
                prevCoords=[X,Y]
                [I, J] = self._shiftVals(index)
                X = X + I
                Y = Y + J

                #print([index, I, J])
                # UAV.idx = [x+I,y+J]
                # add the index to the IDXplan
                UAV.IDXplan.append([UAV.IDXplan[n][0] + I, UAV.IDXplan[n][1] + J])
                #print([X, Y])
                # add the GPS point to the UAV's path
                UAV.GPSplan.append(self._getGPSfromIDX(UAV.IDXplan[n + 1]))
                # uncomment the following to plot each step of the planner
            #print(UAV.IDXplan)
            print(str(len(UAV.IDXplan))+' steps planned')
            self.steps_planned = len(UAV.IDXplan)


    def _init_planner(self, UAV):
        if len(UAV.GPSpath) == 0:
            print("GPS minimum coords")
            GPSmin = self._getGPSfromIDX([UAV.minX,UAV.minY])
            print(GPSmin)
            print("GPS maximum coords")
            GPSmax = self._getGPSfromIDX([UAV.maxX-1,UAV.maxY-1])
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

        if not (UAV.idx):
            # if there isn't already an index, get it from the 'GPS' value
            UAV.idx = self._findGPSindex(UAV.GPS)
            # and append the index to the planned path
            UAV.IDXplan.append(UAV.idx)
        else:
            # otherwise just add the index to the planned path
            UAV.IDXplan.append(UAV.idx)
        return plan_horizon


    def COSPlan(self,UAV):
        plan_horizon = self._init_planner(UAV)
        print("plan horizon: "+str(plan_horizon))
        print("planning...")
        if (UAV.idx):
            [X,Y]=UAV.idx
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


            step=0
            while step < plan_horizon:
                max_sum = -1000000000000
                x=X+I
                y=Y+J
                max_coords=[x,y]
                while x < UAV.maxX and x >= UAV.minX \
                    and y < UAV.maxY and y >= UAV.minY:
                    #print(".")
                    if not np.isnan(self.score_map[x][y]):
                        #print(",")
                        for i in range(8):
                            #print(":")
                            if i != heading-4 and i != heading+4 \
                                and i != heading:
                                #print("direction: "+str(i))
                                total=self._total_line(UAV, deepcopy(x),deepcopy(y),i)

                                if total=='nan':
                                    pass
                                elif total>max_sum:
                                    #print("_")
                                    max_sum=total
                                    max_dir=i
                                    max_coords=[x,y]

                            else:
                                total = 0
                        x=x+I
                        y=y+J
                    else:
                        #print("nope1: "+str([x,y]))
                        break
                #print("tested coords: "+str(max_coords))
                #print("current coords: "+str([X,Y]))
                if(max_coords==[X,Y]):
                    print("BROKEN: "+str([X,Y,heading,step]))
                    break
                #print("plan horizon: "+str(plan_horizon))
                #print("dir1, dir2, pos, next, step, sum:"+str([heading, max_dir, [X,Y], max_coords, step, max_sum]))
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
                    step=step+1
                #print("new coords: "+str([X, Y]))
                heading = max_dir
                [I, J] = self._shiftVals(max_dir)
            #print(UAV.IDXplan)
            print("steps planned: "+str(len(UAV.IDXplan)))
            self.steps_planned=len(UAV.IDXplan)






    def _total_line(self, UAV, X,Y, d):
        [I, J]=self._shiftVals(d)
        X=X+I
        Y=Y+J
        i=1
        total = 0
        while X<UAV.maxX and X>=UAV.minX\
                and Y<UAV.maxY and Y>=UAV.minY:
            if not np.isnan(self.score_map[X][Y]):
                total = total + self.checkScore(UAV,[X, Y], None)
                i=i+1
            else:
                if i==1:
                    #print("nope2: "+str([X, Y, d]))

                    return 'nan'
                else:
                    #print("EOL: "+str(total))
                    return total
            X = X + I
            Y = Y + J
        #print("EOL: " + str(total))
        #print([X, Y])
        if i>1:
            return total
        else:
            return 'nan'

    def checkScore(self,UAV,coord, path_mask=None):
        X=coord[0]
        Y=coord[1]
        if path_mask == None:
            path_mask = UAV.plan_mask
        try:
            #print(path_mask[X][Y])
            return self.score_map[X][Y] * path_mask[X][Y] - UAV.maskSUB * (1 - UAV.plan_mask[X][Y])

        except:
            return None

    def quadSearchPlan(self, UAV):
        plan_horizon = self._init_planner(UAV)
        print("plan horizon: " + str(plan_horizon))
        print("planning...")
        if (UAV.idx):
            [X,Y]=UAV.idx
            #print([X, Y])
            # calculate the current score map
            self.calcScoreMap()
            UAV.update_mask(UAV.plan_mask, [X, Y])
            self.avgQuad([X,Y])



    def avgQuad(self,XY):
        X = XY[0]
        Y = XY[1]
        Q1=0
        Q2=0
        Q3=0
        Q4=0
        for i in range(X,self.vman.grid_res[0],1):
            for j in range(Y,self.vman.grid_res[1],1):
                if not np.isnan(self.score_map[i][j]):
                    Q1=Q1+self.score_map[i][j]
        for i in range(0,X,1):
            for j in range(Y,self.vman.grid_res[1],1):
                if not np.isnan(self.score_map[i][j]):
                    Q2=Q2+self.score_map[i][j]
        for i in range(0,X,1):
            for j in range(0,Y,1):
                if not np.isnan(self.score_map[i][j]):
                    Q3=Q3+self.score_map[i][j]
        for i in range(X,self.vman.grid_res[0],1):
            for j in range(0,Y,1):
                if not np.isnan(self.score_map[i][j]):
                    Q4=Q4+self.score_map[i][j]
        print(Q1)
        print(Q2)
        print(Q3)
        print(Q4)
        TOT=Q1+Q2+Q3+Q4
        I = (Q1-Q2-Q3+Q4)/TOT
        J = (Q1+Q2-Q3-Q4)/TOT
        if round(I)==0 and round(J)==0:
            g=2
        if not np.isnan(self.score_map[X+round(I)][Y+round(J)]):
            return round(I),round(J)
        elif not np.isnan(self.score_map[X + round(I)][Y]):
            return


    def updateEstimatesTV(self, UAV):
        print("Updating Estimates...")
        # current field estimate
        Xk = deepcopy(self.vman.flowfield.u_field[:, :, self.plane])
        # current velocity and direction estimates
        vd_k = [[deepcopy(UAV.v0)],
                [deepcopy(UAV.d0)]]
        # current sensitivity matrix
        sens_mat = np.column_stack([self.sens_mat[0].flatten(), self.sens_mat[1].flatten()])
        # apply the pseudoinverse
        sens_mat_pinv = np.linalg.pinv(sens_mat)
        # field actual

        # if the node is not masked by the UAV's path, use the estimate
        self.mask_size_history.append(1-float(np.sum(UAV.path_mask))/ \
                    (self.vman.grid_res[0]*self.vman.grid_res[1]))
        temp_Xbar = deepcopy(self.Xbar[self.steps_taken])
        for i in range(len(UAV.IDXpath)):
            XY=UAV.IDXpath[i]
            temp_Xbar[XY[0]][XY[1]]=UAV.measure[i]

        # convert Xk and Xbar to column vectors

        print("Updating Estimates...")
        print("Prior estimate: " + str([UAV.v0, UAV.d0]))
        # current field estimate
        temp_Xk = deepcopy(self.vman.flowfield.u_field[:, :, self.plane]).flatten()
        # current velocity and direction estimates
        vd_k = [[deepcopy(UAV.v0)],
                [deepcopy(UAV.d0)]]
        # current sensitivity matrix
        temp_sens_mat0 = self.sens_mat[0].flatten()
        temp_sens_mat1 = self.sens_mat[1].flatten()

        temp_mask = deepcopy(np.array(UAV.path_mask).flatten())
        # if the node is not masked by the UAV's path, use the estimate
        self.mask_size_history.append(1 - float(np.sum(UAV.path_mask)) / \
                                      (self.vman.grid_res[0] * self.vman.grid_res[1]))
        Xbar = list()
        sens_mat0 = list()
        sens_mat1 = list()
        Xk = list()
        IDM = list()
        temp_score_map = deepcopy(self.score_map).flatten()
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
        # print("Xbar shape")
        # print(Xbar.shape)
        # print("Xk shape")
        # print(Xk.shape)
        # print("sens_mat0 shape")
        # print(sens_mat0.shape)
        # print("sens_mat1 shape")
        # print(sens_mat1.shape)
        # apply the pseudoinverse
        sens_mat_pinv = np.linalg.pinv(np.column_stack([sens_mat0, sens_mat1]))
        # print("sens_mat_pinv shape")
        # print(sens_mat_pinv.shape)
        # calculate the next estimate: Vdk+1 = Vdk + sens_mat_pinv*(Xbar-Xk)
        vd_kp1 = vd_k + np.matmul(sens_mat_pinv, Xbar - Xk)
        # update the estimates
        UAV.v0 = vd_kp1[0][0]
        UAV.d0 = vd_kp1[1][0]
        # update the wind speed estimate
        self.WF['farm']['properties']['wind_speed'] = UAV.v0
        # calculate the FLORIS model with the new speed
        vmanTemp = VisualizationManager(self.WF, self.vman.grid_res)
        # update the direction estimate
        vmanTemp.flowfield.wind_direction = UAV.d0
        # recalculate the FLORIS wake model with the new direction estimate
        turbines = [turbine for _, turbine in vmanTemp.flowfield.turbine_map.items()]
        for k, turbine in enumerate(turbines):
            turbine.yaw_angle = UAV.d0
        vmanTemp.flowfield.calculate_wake()
        # update the sensitivity matrix
        self.sens_mat = deepcopy(vmanTemp.calcSensitivityMatrix(UAV.v0,UAV.d0))
        # update the score map
        self.calcScoreMap()
        # output the error to the terminal
        print("Error: [v,theta]="+str(self.error[len(self.error) - 1]))
        # update the local VisualManager object
        self.vman = deepcopy(vmanTemp)
        # keep track of the error signals for plotting
        self.error.append([self.params.vBar - UAV.v0,
                           self.params.dBar - UAV.d0])

    def getMeasure(self, XY):
        measure = self.Xbar[self.steps_taken][XY[0]][XY[1]]
        self.steps_taken=self.steps_taken+1
        if self.steps_taken >= self.Xbar.shape[0]:
            print("looping to beginning of data")
            self.steps_taken = 0
        return measure
