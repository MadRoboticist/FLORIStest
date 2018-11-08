## \file pathPlan.py
# This file contains functions which handle UAV path planning using Floris data.
#

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
    ## @var history
    # keeps track of everything that happens so it can be plotted afterward.
    #
    # This is how history is populated:
    # self.history.append(deepcopy([self.score_map,
    #                                UAV.GPSpath,
    #                                UAV.GPSplan,
    #                                UAV.wave_map,
    #                               self.error[len(self.error)-1],
    #                               self.params.v0,
    #                               self.params.d0]))
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
        # recalculate the wake with the new direction
        vman_bar.flowfield.calculate_wake()
        ## @var Xbar
        #   this is the 'actual' u_field
        self.Xbar = deepcopy(vman_bar.flowfield.u_field[:, :, self.plane])
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
        ## the local copy of the sensitivity matrix
        self._sens_mat = self.vman.calcSensitivityMatrix()
        ## @var error
        # holds a record of the direction and velocity errors for plotting
        self.error = list()
        # append the initial error
        self.error.append([self.params.vBar-self.params.v0,
                           self.params.dBar-self.params.d0])
        # delete the old vman objects
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
        ## @var score_map
        # the score map is computed as ||df/dd||*||df/dv||
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
                    # we test the score using the UAV's path mask
                    if self.score_map[i][j]*UAV.plan_mask[i][j]>max:
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
        # start the planner fresh
        UAV.reset_planner()
        # figure out what the index of the starting point is
        if len(self.history) < UAV.init_mask_size:
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
        # if the GPS point was valid, the indices should exist
        if (UAV.idx):
            # calculate the current score map
            self._calcScoreMap()
            # let the iteration begin
            for i in range(plan_horizon):
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
                # uncomment the following to plot each step of the planner
                '''
                self.history.append(deepcopy([self.score_map,
                                              UAV.GPSpath,
                                              UAV.GPSplan,
                                              UAV.wave_map,
                                              self.error[len(self.error)-1],
                                              self.params.v0,
                                              self.params.d0]))
                                              '''
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
                headCost = abs(i*np.deg2rad(45)- \
                                  (UAV.plan_heading[len(UAV.plan_heading)-1])*
                                     np.deg2rad(45))
                # normalize to a value between -PI and PI
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
        # return the decision and the resulting score
        return max, dir

    # updateEstimates
    ## A function which updates the direction and speed estimates
    # based on the UAV's path mask over Xbar
    #
    # @param UAV The UAV whose path is being used to update the estimates
    def updateEstimates(self, UAV):
        # current field estimate
        Xk = deepcopy(self.vman.flowfield.u_field[:,:,self.plane])
        # current velocity and direction estimates
        vd_k = [[deepcopy(self.params.v0)],
                [deepcopy(self.params.d0)]]
        # current sensitivity matrix
        sens_mat = np.column_stack([self._sens_mat[0].flatten(), self._sens_mat[1].flatten()])
        # apply the pseudoinverse
        sens_mat_pinv = np.linalg.pinv(sens_mat)
        # field actual
        Xbar = deepcopy(self.Xbar)
        # if the node is not masked by the UAV's path, use the estimate
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                if UAV.path_mask[i][j]==1:
                    Xbar[i][j]=Xk[i][j]
        # convert Xk and Xbar to column vectors
        Xbar = np.vstack(Xbar.flatten())
        Xk = np.vstack(Xk.flatten())
        # calculate the next estimate: Vdk+1 = Vdk + sens_mat_pinv*(Xbar-Xk)
        vd_kp1=vd_k+np.matmul(sens_mat_pinv,Xbar-Xk)
        # update the estimates
        self.params.v0 = vd_kp1[0][0]
        self.params.d0 = vd_kp1[1][0]
        # update the wind speed estimate
        self.WF['farm']['properties']['wind_speed'] = self.params.v0
        # calculate the FLORIS model with the new speed
        vmanTemp = VisualizationManager(self.WF, self.vman.grid_res)
        # update the direction estimate
        vmanTemp.flowfield.wind_direction = self.params.d0
        # recalculate the FLORIS wake model with the new direction estimate
        vmanTemp.flowfield.calculate_wake()
        # update the sensitivity matrix
        self._sens_mat = deepcopy(vmanTemp.calcSensitivityMatrix())
        # update the score map
        self._calcScoreMap()
        # output the error to the terminal
        print(self.error[len(self.error) - 1])
        # update the local VisualManager object
        self.vman = deepcopy(vmanTemp)
        # keep track of the error signals for plotting
        self.error.append([self.params.vBar-self.params.v0,
                           self.params.dBar-self.params.d0])

    # updateEstimates
    ## A function which updates the direction and speed estimates
    # based on the UAV's path mask over Xbar
    #
    # @param UAV The UAV whose path is being used to update the estimates
    def plotHistory(self, UAV, filename=None):
        # variable which tells us if we are looking at the Wave or the Score map
        self._map_sel = 0
        # boolean which tells us if the animation is playing
        self._play = True
        # create a figure
        f = plt.figure(figsize=(10,7.5))
        # separate it into a grid
        gs = f.add_gridspec(2,3)
        # put the error plots in the leftmost figures
        ax_v = f.add_subplot(gs[0,0], title='speed error')
        ax_d = f.add_subplot(gs[1,0])
        # the path plot takes up the right 2/3 of the figure
        axbig = f.add_subplot(gs[0:,1:])
        # plot the velocity error
        ax_v.plot([i for i in range(len(self.error))], [i[0] for i in self.error])
        # format x,y and c for scatterplotting
        x = deepcopy(self.X)
        y = deepcopy(self.Y)
        c = deepcopy(self.history[1][0])
        # plot the score map
        cont = axbig.scatter(x=x.flatten(), y=y.flatten(), c=c.flatten(), cmap='gnuplot2')
        plt.colorbar(cont, ax=axbig)
        # plot the direction error
        ax_d.plot([i for i in range(len(self.error))], [np.rad2deg(i[1]) for i in self.error])
        # the plot's title holds a lot of info
        plt.suptitle("Wind Speed and Direction Estimates with a UAV for Sensing\n" +
                     'Initial Speed estimate: ' + str(self.history[0][5]) +
                     ' m/s, Actual: ' + str(self.params.vBar) +
                     ' m/s\nInitial Direction estimate: ' + str(
            np.rad2deg(self.history[0][6])) + '\N{DEGREE SIGN}' +
                     ', Actual: ' + str(np.rad2deg(self.params.dBar)) + '\N{DEGREE SIGN}' +
                     '\nPlanning Horizon: ' + str(UAV.plan_horizon) +
                     '          Moves until recalculation: ' + str(UAV.moves2recalc) +
                     '          Memory: ' + str(UAV.patrolMax) + ' nodes'
                     '\nFinal error: e\N{GREEK SMALL LETTER THETA} = ' +
                     str(self.error[len(self.error) - 1][0]) + '\N{DEGREE SIGN}' +
                     '          ev = ' + str(self.error[len(self.error) - 1][1]) + ' m/s')
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
        sld = Slider(sld_ax,
                     'moves',
                     0, len(self.history) - 1, valinit=0)
        # set the initial slider text
        sld.valtext.set_text('move 0')
        # button axis
        btn_ax = plt.axes([0.85, 0.925, 0.125, 0.05])
        # create a button which toggles between score map and wave map
        btn = Button(btn_ax, 'Show Wave')
        # button axis
        btn2_ax = plt.axes([0.85, 0.01, 0.125, 0.05])
        # create a button which plays/pauses the animation
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
            axbig.scatter(x=x.flatten(), y=y.flatten(), c=self.history[idx][self._map_sel].flatten(), cmap='gnuplot2')
            # plot the plan
            axbig.plot([i[0] for i in self.history[idx][2]],
                       [i[1] for i in self.history[idx][2]], color='lime')
            # plot the actual path
            axbig.plot([i[0] for i in self.history[idx][1]],
                             [i[1] for i in self.history[idx][1]], color = 'red')
            # for shorthand, UAV's location on the map
            _x = self.history[idx][1][len(self.history[idx][1])-1][0]
            _y = self.history[idx][1][len(self.history[idx][1])-1][1]
            # plot the UAV, a big red circle
            axbig.plot(_x,_y, marker='o', markersize=10, color="red")
            # plot the velocity error
            ax_v.plot([i for i in range(len(self.error))], [i[0] for i in self.error])
            # add a vertical red line to show where we are in the iterations
            if idx > 100:
                ax_v.axvline((idx-UAV.init_mask_size)/(UAV.moves2recalc+1), color='red')
            else:
                ax_v.axvline(0, color='red')
            # plot the direction error
            ax_d.plot([i for i in range(len(self.error))], [np.rad2deg(i[1]) for i in self.error])
            # add a vertical red line to show where we are in the iterations
            if idx > 100:
                ax_d.axvline((idx-UAV.init_mask_size) / (UAV.moves2recalc + 1), color='red')
            else:
                ax_d.axvline(0, color='red')
            # plot the turbines
            for coord, turbine in self.vman.flowfield.turbine_map.items():
                a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
                b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
                a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
                b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
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
        an = anim.FuncAnimation(f, animate, interval=100,frames=len(self.history))
        # render to video. to make it play faster, increase fps
        if filename:
            an.save(filename+'.mp4',fps=15,dpi=300)
        # show the plot
        plt.show()

