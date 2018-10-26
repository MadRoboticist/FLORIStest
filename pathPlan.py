## \file pathPlan.py
# This file contains functions which handle UAV path planning using Floris data.
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from floris.coordinate import Coordinate
from visualization_manager_DJ import VisualizationManager
from math import pi as PI
import math
from copy import deepcopy

## \class dir
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
    history = list()

    ## Class constructor
    #   @param vman a VisualManager object
    def __init__(self, vman):
        ## @var vman
        #   A VisualizationManager object
        self.vman = deepcopy(vman)
        ## @var params
        #   A local copy of vman.params
        self.params = deepcopy(vman.params)
        ## @var plane
        # the plane where the computations and visualizations are taking place
        self.plane = int(self.vman.flowfield.grid_resolution.z * self.params.percent_height)
        # time to calculate Xbar
        vman.WF['farm']['properties']['wind_speed'] = self.vman.params.vBar
        vman_bar = VisualizationManager(vman.WF, vman.grid_res)
        ## @var Xbar
        #   this is the 'actual' u_field
        self.Xbar = deepcopy(vman_bar.flowfield.u_field[:, :, self.plane])
        ## @var WF
        #  a local copy of vman.WF
        self.WF = deepcopy(vman.WF)
        ## @var X_map
        #   A map of all possible transitions on the map and their respective scores
        self.X_map = [[_X_map_node() for x in range(self.vman.grid_res[0])]
                 for y in range(vman.grid_res[1])]
        ## @var X
        #   An X grid mesh of the map's 'GPS' coordinates
        self.X = deepcopy(vman.flowfield.x[:, :, self.plane])
        ## @var Y
        #   A Y grid mesh of the map's 'GPS' coordinates
        self.Y = deepcopy(vman.flowfield.y[:, :, self.plane])
        ## the initial sensitivity matrix
        self._sens_mat = self.vman.calcSensitivityMatrix()
        self.sens_mat = list()
        self.sens_mat.append(self._sens_mat)
        # delete the old vman objects
        self.error = list()
        self.error.append([self.params.vBar-self.params.v0,
                           self.params.dBar-self.params.d0])
        del vman, vman_bar
        # calculate the score map
        self._calcScoreMap()
        # calculate the cost map
        self._update_cost_map();

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
    def _euclidDist(self, P1, P2):
        # sqrt((Y2-Y1)^2+(X2-X1)^2)
        return math.sqrt(((P2[0]-P1[0])**2)+((P2[1]-P1[1])**2))

    ## a method to update the transition/cost map
    def _update_cost_map(self):
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

    ## a method to calculate the score map from the current sensitivity matrix
    def _calcScoreMap(self):
        # compute the normalized df/dd column of the sensitivity matrix ||df/dd||
        Zd = abs(self._sens_mat[1]) / np.amax(abs(self._sens_mat[1]))
        # compute the normalized df/dv column of the sensitivity matrix ||df/dv||
        Zv = abs(self._sens_mat[0]) / np.amax(abs(self._sens_mat[0]))
        # the score is computed as ||df/dd||*||df/dv||
        self.score_map = Zd * Zv
        # iterate through the entire score map
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                # create a no-fly zone around each turbine
                for coord, turbine in self.vman.flowfield.turbine_map.items():
                    # compute if the euclidean distance is less
                    # than the minimum distance threshold minDist2turbine
                    if self._euclidDist([self.X[i][j],
                                         self.Y[i][j]],
                                        [coord.x, coord.y]) < \
                            self.minDist2turbine:
                        # if it is, make it untraversable
                        self.score_map[i][j] = None
        # calculate a wave map from this score map

    ## calculates a wave map using the highest score on the score map
    def _calcWaveMap(self, UAV):
        # set the max unobtainably low (since our range is 0-1)
        max = -100
        idx = [0,0]
        # set the wave_map to all ones
        UAV.wave_map = [[1 for i in range(self.vman.grid_res[0])] \
                          for j in range(self.vman.grid_res[1])]
        # iterate through the entire map
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                # if we aren't using a UAV for this calculation
                # (such as for the initial wave map calculation)
                # the program will throw an exception
                try:
                    # we test the score using the UAV's mask
                    if self.score_map[i][j]*UAV.plan_mask[i][j]>max:
                        max=self.score_map[i][j]
                        idx = [i,j]
                except:
                    if self.score_map[i][j]>max:
                        max=self.score_map[i][j]
                        idx = [i,j]
        x = idx[0]
        y = idx[1]
        # iterate through the map again
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                # set the value of the wave map to the euclidean
                # distance from the maximum score
                UAV.wave_map[i][j] = self._euclidDist([x,y],[i,j])
        # normalize the wave values to 0-1,
        # but with increasing values as they approach the maximum
        UAV.wave_map = 1-UAV.wave_map/np.amax(UAV.wave_map)

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
            axx = ax.scatter(x=self.X, y=self.Y, c=UAV.wave_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx, ax=ax)
        else:
            axx = plt.scatter(x=self.X, y=self.Y, c=UAV.wave_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx)

        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
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
            axx = ax.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx, ax=ax)
        else:
            axx = plt.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx)
        # the color bar

        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
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
            axx = ax.scatter(x=self.X, y=self.Y, c=mask)
            plt.colorbar(axx, ax=ax)
        else:
            axx = plt.scatter(x=self.X, y=self.Y, c=mask)
            plt.colorbar(axx)
        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            if ax:
                ax.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')
            else:
                plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')

    ## Generates a plot of a given UAV's path over the score map
    #   @param UAV the UAV whose path is to be plotted
    #   @param ax Can be passed to the function if the plot
    #           is to be in a subplot
    def plotScoreMapUAV(self, UAV, ax=None):
        if ax:
            axx = ax.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx, ax=ax)
        # a scatter plot of the score map
        else:
            axx = plt.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
            plt.colorbar(axx)
        # the color bar

        # the title
        plt.suptitle("Greedy UAV path\n"+str(len(UAV.GPSplan)-1)+" moves")
        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            if ax:
                ax.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
            else:
                plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
        # plot the UAV's path
        if ax:
            ax.plot([i[0] for i in UAV.GPSplan],
                     [i[1] for i in UAV.GPSplan], color='lime')
            ax.plot([i[0] for i in UAV.GPSpath],
                     [i[1] for i in UAV.GPSpath], color='red')
        else:
            plt.plot([i[0] for i in UAV.GPSplan],
                    [i[1] for i in UAV.GPSplan], color='lime')
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
            print("Error: GPS point is out of bounds")
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

        UAV.reset_planner()
        # figure out what the index of the starting point is
        if not (UAV.idx):
            UAV.idx = self._findGPSindex(UAV.GPS)
            # and append the index to the new path
            UAV.IDXplan.append(UAV.idx)
        else:
            UAV.IDXplan.append(UAV.idx)
        # if the GPS point was valid, the indices should exist
        if (UAV.idx[0] and UAV.idx[1]):
            # calculate the current score map
            self._calcScoreMap()
            # let the iteration begin
            for i in range(UAV.plan_horizon):
                # just for shorthand notation
                x=UAV.IDXplan[i][0]
                y=UAV.IDXplan[i][1]

                # if the indices are invalid, an exception will be thrown
                try:
                    # start by updating the UAV's score mask for the current position
                    UAV.update_mask(UAV.plan_mask, [x,y])
                    # calculate the new wave map with the updated mask
                    self._calcWaveMap(UAV)
                    # calculate the next step
                    max, dir = self._greedyStep(UAV)
                except:
                    # if we run off the map, or into a no-fly zone, end the run
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
                self.history.append(deepcopy([self.score_map,
                                              UAV.GPSpath,
                                              UAV.GPSplan,
                                              UAV.wave_map,
                                              self.error[len(self.error)-1],
                                              self.params.v0,
                                              self.params.d0]))
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
            y>=self.vman.grid_res[1]:
            return None
        # set an unobtainably low maximum value
        max = -10000000
        # set the direction to zero just for starters
        dir = 0
        # if there isn't a transition list for this node
        # then the node doesn't exist, and any move is invalid
        try:
            Xition = self.X_map[x][y].Xitions
        except:
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
                headCost = abs(i*np.deg2rad(45)- \
                                  (UAV.plan_heading[len(UAV.plan_heading)-1])*
                                     np.deg2rad(45))
                if headCost > PI:
                    headCost = abs(headCost - 2 * PI)
                # calculate the 2nd derivative of the score wrt the transition
                d2S = Xition[i].dSscore-UAV.dSplan[len(UAV.dSplan)-1]
                # calculate the reward of the move
                reward = (UAV.scoreWt*score + \
                       UAV.dSwt*Xition[i].dSscore + \
                       UAV.d2Swt*d2S - \
                       UAV.windWt*Xition[i].thetaCost/PI - \
                       UAV.headWt*headCost/PI + \
                       UAV.waveWt*wave + \
                       UAV.timeWt*len(UAV.IDXplan)) * \
                       mask

                # we are trying to maximize the reward
                if reward > max and \
                    x+I<self.vman.grid_res[0] and \
                    y+J<self.vman.grid_res[1] and \
                    x+I>=0 and y+J>=0:
                    # if the transition is valid and has the highest reward
                    # then track the reward score and the direction of the transition
                    max = reward
                    dir = i
            # if the transition is invalid, skip it
            except:
                pass
        return max, dir

    def updateEstimates(self, UAV):

        Xk = deepcopy(self.vman.flowfield.u_field[:,:,self.plane])
        vd_k = [[self.params.v0],
                [self.params.d0]]
        sens_mat = np.column_stack([self._sens_mat[0].flatten(), self._sens_mat[1].flatten()])
        sens_mat_pinv = np.linalg.pinv(sens_mat)
        Xbar = deepcopy(self.Xbar)
        counter=0
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                if UAV.path_mask[i][j]==1:
                    counter += 1
                    Xbar[i][j]=Xk[i][j]
        Xbar = np.vstack(Xbar.flatten())
        Xk = np.vstack(Xk.flatten())
        # Vdk+1 = Vdk + sens_mat_pinv*(Xbar-Xk)
        [[self.params.v0],
         [self.params.d0]]=vd_k+np.matmul(sens_mat_pinv,Xbar-Xk)

        self.WF['farm']['properties']['wind_speed'] = self.params.v0
        vmanTemp = VisualizationManager(self.WF, self.vman.grid_res)
        vmanTemp.params.d0 = self.params.d0
        vmanTemp.flowfield.calculate_wake()
        self._sens_mat = deepcopy(vmanTemp.calcSensitivityMatrix())
        self.sens_mat.append(self._sens_mat)
        self._calcScoreMap()
        self.vman = deepcopy(vmanTemp)
        self.error.append([self.params.vBar-self.params.v0,
                           self.params.dBar-self.params.d0])
        print(self.error[len(self.error)-1])

    def plotHistory(self, UAV):
        steps = [i for i in range(len(self.history))]
        print(self.history[20][4])
        f, axarr = plt.subplots(2,2,gridspec_kw = {'width_ratios':[1, 3]})
        sld_ax = plt.axes([0.23, 0.02, 0.56, 0.02])
        sld = Slider(sld_ax,
                     'steps',
                     0, len(self.history) - 1, valinit=0)
        sld.valtext.set_text('step 0')
        f.suptitle('Speed estimate: ' + str(self.params.v0) +
                   ' m/s, Actual: ' + str(self.params.vBar) +
                   ' m/s\nDirection estimate: ' + str(np.rad2deg(self.params.d0)) + '\N{DEGREE SIGN}' +
                   ', Actual: ' + str(np.rad2deg(self.params.dBar)) + '\N{DEGREE SIGN}' +
                   '\nPlanning Horizon: ' + str(UAV.plan_horizon) +
                   ' Total steps planned: ' + str(len(self.history)) +
                   '\nFinal error: e\N{GREEK SMALL LETTER THETA} = ' +
                   str(self.error[len(self.error)-1][0]) +
                   ' ev = ' + str(self.error[len(self.error)-1][1]))
        axarr[0][0].plot([i for i in range(len(self.error))], [i[0] for i in self.error])
        axarr[0][0].set_title('speed error')
        v = np.linspace(0, np.amax([i[0] for i in self.history]), 100)
        V = np.linspace(0, np.amax([i[0] for i in self.history]), 5)
        x = deepcopy(self.X)
        y = deepcopy(self.Y)
        c = deepcopy(self.history[1][0])
        cont = axarr[0][1].scatter(x=x.flatten(), y=y.flatten(), c=c.flatten(), cmap='gnuplot2')
        # cont = axarr[0][2].scatter(x=X[:, :, plane].flatten(), y=Y[:, :, plane].flatten(), c=err_field[0].flatten(),
        #                    cmap='gnuplot2')
        plt.colorbar(cont, ax=axarr[0][1])
        # cb = plt.colorbar(cont, ax=axarr[0][2])
        # cb.set_clim(vmin=0, vmax=np.amax([i[0] for i in self.history]))
        # cb.set_ticks(V, True)
        # cb.set_label('||df/dd||*||df/dv||')
        # cb.draw_all()
        axarr[0][1].set_title('path over score map')
        axarr[1][0].plot([i for i in range(len(self.error))], [np.rad2deg(i[1]) for i in self.error])
        axarr[1][0].set_title('direction error')
        axarr[1][0].set_xlabel('steps')
        u = np.linspace(0, 1, 100)
        U = np.linspace(0, 1, 5)
        c = deepcopy(self.history[0][3])
        cont2 = axarr[1][1].scatter(x=x.flatten(), y=y.flatten(), c=c.flatten(), cmap='gnuplot2')
        plt.colorbar(cont2, ax=axarr[1][1])
        '''
        cb2 = plt.colorbar(mappable=cont, ax=axarr[1][2])
        cb2.set_clim(vmin=0, vmax=1)
        cb2.set_ticks(U)
        cb2.draw_all()'''
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            axarr[0][1].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
            axarr[1][1].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
        axarr[1][1].set_title('wave map')

        # f.tight_layout()
        plt.subplots_adjust(left=0.05,
                            bottom=0.15,
                            right=0.95,
                            top=0.83,
                            wspace=0.27,
                            hspace=0.19)

        def update_plot(val):

            idx = int(round(sld.val))
            sld.valtext.set_text('iteration ' + '{}'.format(idx))
            axarr[0][1].clear()
            axarr[0][1].scatter(x=x.flatten(), y=y.flatten(), c=self.history[idx][0].flatten(), cmap='gnuplot2')
            axarr[0][1].plot([i[0] for i in self.history[idx][2]],
                             [i[1] for i in self.history[idx][2]], color = 'lime')
            axarr[0][1].plot([i[0] for i in self.history[idx][1]],
                             [i[1] for i in self.history[idx][1]], color = 'red')
            # axarr[0][2].scatter(x=X[:, :, plane].flatten(), y=Y[:, :, plane].flatten(), c=err_field[idx].flatten(), cmap='gnuplot2')
            axarr[1][1].clear()
            axarr[1][1].scatter(x=x.flatten(), y=y.flatten(), c=self.history[idx][3].flatten(), cmap='gnuplot2')
            axarr[1][1].plot([i[0] for i in self.history[idx][2]],
                             [i[1] for i in self.history[idx][2]], color='lime')
            axarr[1][1].plot([i[0] for i in self.history[idx][1]],
                             [i[1] for i in self.history[idx][1]], color='red')
            for coord, turbine in self.vman.flowfield.turbine_map.items():
                a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
                b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
                a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
                b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
                axarr[0][1].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
                axarr[1][1].plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
            plt.draw()

        sld.on_changed(update_plot)
        plt.subplots_adjust(left=0.05,
                            bottom=0.11,
                            right=1.0,
                            top=0.84,
                            wspace=0.15,
                            hspace=0.24)
        plt.show()

