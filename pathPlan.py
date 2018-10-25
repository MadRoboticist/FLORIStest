## \file pathPlan.py
# This file contains functions which handle UAV path planning using Floris data.
#

import numpy as np
import matplotlib.pyplot as plt
from floris.coordinate import Coordinate
from math import pi as PI
import math
import matplotlib as mpl
import time
from copy import deepcopy


## \class dir
# An enumeration of the cardinal directions
# starting from East and moving counterclockwise.
#
# This enumeration can be used to generate compass
# angles, e.g. NW = 3x45 degrees or 3xPI/4
class dir(enumerate):
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

## a structure which holds information for each possible transition from
#   one position on the map to another
class _Transition:
    ## Class constructor
    def __init__(self):
        ## @var thetaCost
        # A cost value
        # for wind direction vs. transition heading
        self.thetaCost = 0.0
        ## @var phiCost
        # A cost value
        # for transition heading vs. UAV heading
        self.phiCost = 0.0
        ## @var dSscore
        # A score/cost value
        # for the rate of change of the total score
        # resulting from the corresponding transition
        self.dSscore = 0.0
        ## @var d2Sscore
        # A score/cost value
        # for the 2nd derivative of score/transition
        self.d2Sscore = 0.0

## a structure which holds information for each node on
# a discretized transition map
class _X_map_node:
    ## Class constructor
    def __init__(self):
        ## @var Xitions
        # A list of the 8 possible transitions
        # from the corresponding node
        self.Xitions = [_Transition() for i in range(8)]
        ## @var GPS
        # The GPS coordinates [lat,lon] of the corresponding node
        self.GPS = [0, 0]
        ## @var idx
        # The coordinates [x, y] of the node
        self.idx = [0, 0]
        ## @var score
        # The value score of the node from the sensitivity matrix
        self.score = None

## \class UAV
# @brief A class for representing the state of a UAV on the windfarm
class UAV:
    ## Class constructor.
    def __init__(self):
        ## @var GPS
        # Holds the "GPS" location
        # of the UAV. Currently, the GPS location is an
        # ordered pair representing the XY position in meters
        # from the map's origin X=GPS[0], Y=GPS[1]
        self.GPS = [0, 0]
        ## @var idx
        # Holds the XY node index from the bottom left
        # corner of the map. X=idx[0], Y=idx[1]
        self.idx = [0, 0]
        ## @var heading
        # Holds the enumerated direction that the UAV
        # is currently traveling. see Class dir()
        self.heading = dir.E
        ## @var IDXplan
        # A list of idx points representing the UAV's planned path
        # through the nodes.
        self.IDXplan = list()
        ## @var GPSplan
        # A list of GPS points representing the UAV's planned path
        # through the map.
        self.GPSplan = list()
        ## @var IDXpath
        # A list of idx points representing the UAV's actual path
        # through the nodes.
        self.IDXpath = list()
        ## @var GPSpath
        # A list of GPS points representing the UAV's actual path
        # through the map.
        self.GPSpath = list()
        ## @var score
        # Holds the score of the UAV's path
        self.score = 0
        ## @var wave_map
        # Holds a map of values which can serve as
        # a potential field, drawing the UAV to
        # the node with the maximum score value
        self.wave_map = None
        ## @var score_mask
        # Holds a map of values which degrade the
        # score of the corresponding node after it has been visited.
        # (applied as calculated_score*mask_value)
        self.score_mask = None
        ## @var maskSUB
        #   sets a subtractive penalty value to be applied to a node's
        #   score mask each time it is visited
        #   (set to 0 to apply no subtractive penalty)
        self.maskSUB = 0
        ## @var maskMUL
        #   sets a multiplicative penalty value to be applied to a node's
        #   score mask each time it is visited
        #   (set to 1 to apply no multiplicative penalty)
        self.maskMUL = 0.5
        ## @var maskMULthenSUB
        #   boolean which decides the order maskSUB and maskMUL are carried out
        self.maskMULthenSUB = True

## \class PathPlanner
#   @brief Contains all path planning functions and parameters
class PathPlanner:

    ## \class params
    #   @brief A class through which parameters can be passed
    #   to the member functions
    class params:
        ## @var minDist2turbine
        # sets a minimum distance from the center of a turbine
        # which serves as a boundary for the UAV
        minDist2turbine = 160
        ## @var scoreWt
        #   sets the weight of the sensitivity score in path calculations
        scoreWt = 10
        ## @var dSwt
        #   sets the weight of the derivative of the sensitivity score
        #   in path calculations
        dSwt = 2
        ## @var windWt
        #   sets the weight of the wind direction vs. transition heading
        #   in path calculations
        windWt = 2
        ## @var waveWt
        #   sets the weight of the wave map in path calculations
        waveWt = 35
        ## @var timeWt
        #   sets the weight of time/iterations elapsed in path calculations
        timeWt = 0.25


    ## Class constructor
    #   @param vman a VisualManager object
    def __init__(self, vman):
        ## @var vman
        #   A VisualizationManager object
        self.vman = vman
        ## the plane where the computations and visualizations are taking place
        plane = int(self.vman.flowfield.grid_resolution.z * self.vman.params.percent_height)
        ## @var X_map
        #   A map of all possible transitions on the map and their respective scores
        self.X_map = [[_X_map_node() for x in range(self.vman.grid_res[0])]
                 for y in range(self.vman.grid_res[1])]
        ## @var X
        #   An X grid mesh of the map's 'GPS' coordinates
        self.X = deepcopy(self.vman.flowfield.x[:, :, plane])
        ## @var Y
        #   A Y grid mesh of the map's 'GPS' coordinates
        self.Y = deepcopy(self.vman.flowfield.y[:, :, plane])
        ## the initial sensitivity matrix
        self._sens_mat = vman.calcSensitivityMatrix()
        ## @var d0
        #   The initial wind direction estimate
        self.d0 = self.vman.params.d0
        ## @var v0
        #   The initial wind speed estimate
        self.v0 = self.vman.params.v0
        self._calcScoreMap()
        self._update_cost_map();

    ## a method to adjust idx based on a cardinal direction
    def _shiftVals(self, val):
        value = {
            dir.E: [1,0],   # one right, zero up
            dir.NE: [1,1],  # one right, one up
            dir.N: [0,1],   # zero right, one up
            dir.NW: [-1,1], # one left, one up
            dir.W: [-1,0],  # one left, zero up
            dir.SW: [-1,-1],# one left, one down
            dir.S: [0,-1],  # zero left, one down
            dir.SE: [1,-1]  # one right, one down
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
                # assign the node a score from the current score map
                self.X_map[i][j].score = self.score_map[i][j]
                # iterate through the 8 possible transitions at each node
                for k in range(8):
                    # calculate the cost of moving through the wind at the transition heading
                    self.X_map[i][j].Xitions[k].thetaCost = k*np.deg2rad(45)-self.d0
                    # adjust the cost to a positive value 0-PI
                    if self.X_map[i][j].Xitions[k].thetaCost > PI:
                        self.X_map[i][j].Xitions[k].thetaCost = abs(k*np.deg2rad(45)-self.d0-2*PI)
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
                            self.params.minDist2turbine:
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
                    if self.score_map[i][j]*UAV.score_mask[i][j]>max:
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
        ax = plt.scatter(x=self.X, y=self.Y, c=UAV.wave_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
        # the colorbar
        plt.colorbar(ax)
        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')

    ## Generates a plot of the current score map
    #   @param ax can be passed to the function if the plot
    #           is to be in a subplot
    def plotScoreMap(self, ax=None):
        # a scatter plot of the score map
        ax = plt.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
        # the color bar
        plt.colorbar(ax)
        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')

    ## Generates a plot of a given UAV's score mask
    #   @param UAV the UAV whose score mask is to be plotted
    #   @param ax Can be passed to the function if the plot
    #           is to be in a subplot
    def plotScoreMask(self, UAV, ax=None):
        # a scatter plot of the UAV's score mask
        ax = plt.scatter(x=self.X, y=self.Y, c=UAV.score_mask)
        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')

    ## Generates a plot of a given UAV's path over the score map
    #   @param UAV the UAV whose path is to be plotted
    #   @param ax Can be passed to the function if the plot
    #           is to be in a subplot
    def plotScoreMapUAV(self, UAV, ax=None):
        # a scatter plot of the score map
        ax = plt.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
        # the color bar
        plt.colorbar(ax)
        # the title
        plt.suptitle("Greedy UAV path\n"+str(len(UAV.GPSplan)-1)+" moves")
        # plot the turbines
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')
        # plot the UAV's path
        plt.plot([i[0] for i in UAV.GPSplan],
                 [i[1] for i in UAV.GPSplan], color='lime')
        # size the plot to fill the window
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
    # with the given number of steps
    #
    # @param UAV The UAV whose path is to be generated
    # @param steps (optional) the number of moves to be generated for the path
    # @param recalc (optional) the number of moves before the score map is recalculated
    # @return UAV.IDXplan a path of indices through which the UAV will travel
    # @return UAV.GPSplan a path of GPS points throught which the UAV will travel
    def greedyPath(self, UAV, steps=100, recalc=10):
        # start with a new GPS path
        UAV.GPSplan.clear()
        # append the starting GPS point
        UAV.GPSplan.append(UAV.GPS)
        # figure out what the index of the starting point is
        UAV.idx = self._findGPSindex(UAV.GPS)
        # also start with a new index path
        UAV.IDXplan.clear()
        # and append the index to the new path
        UAV.IDXplan.append(UAV.idx)
        # start the score at zero
        UAV.score=0
        # if the GPS point was valid, the indices should exist
        if (UAV.idx[0] and UAV.idx[1]):
            # set the UAV's score mask to all ones
            UAV.score_mask = [ [ 1 for i in range(len(self.score_map[0]))] \
                             for j in range(len(self.score_map))]
            # calculate the current score map
            self._calcScoreMap()
            # let the iteration begin
            for i in range(steps):
                    # just for shorthand notation
                    x=UAV.IDXplan[i][0]
                    y=UAV.IDXplan[i][1]

                    # if the indices are invalid, an exception will be thrown
                    try:
                        # start by updating the UAV's score mask for the current position
                        if(UAV.maskMULthenSUB):
                            # multiply then subtract
                            UAV.score_mask[x][y] *= UAV.maskMUL
                            UAV.score_mask[x][y] -= UAV.maskSUB
                        else:
                            # subtract then multiply
                            UAV.score_mask[x][y] -= UAV.maskSUB
                            UAV.score_mask[x][y] *= UAV.maskMUL
                        # calculate the new wave map with the updated mask
                        self._calcWaveMap(UAV)
                        # calculate the next step
                        max, dir = self._greedyStep(UAV)
                    except:
                        # if we run off the map, or into a no-fly zone, end the run
                        return
                    # add the max to the UAV's score
                    UAV.score += max
                    # update the UAV's heading
                    UAV.heading = dir
                    # update idx based on the direction of travel
                    [I, J] = self._shiftVals(dir)
                    # UAV.idx = [x+I,y+J]
                    # add the index to the IDXplan
                    UAV.IDXplan.append([UAV.IDXplan[i][0]+I, UAV.IDXplan[i][1]+J])
                    # set the UAV's new GPS point
                    # UAV.GPS = self._getGPSfromIDX(UAV.idx)
                    # add the GPS point to the UAV's path
                    UAV.GPSplan.append(self._getGPSfromIDX(UAV.IDXplan[i+1]))
        # the starting point was invalid
        else:
            print("Aborting.")

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
                mask = UAV.score_mask[x+I][y+J]
                # pull the wave map value for the transition node
                wave = UAV.wave_map[x+I][y+J]
                # calculate the reward of the move
                reward = (self.params.scoreWt*score + \
                       self.params.dSwt*Xition[i].dSscore - \
                       self.params.windWt*Xition[i].thetaCost/PI + \
                       self.params.waveWt*wave + \
                       self.params.timeWt*len(UAV.IDXplan)) * \
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
