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


## an enumeration of the cardinal directions
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
        # for wind direction vs. UAV heading
        self.thetaCost = 0.0
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

## A class for representing the state of a UAV on the windfarm
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
        ## @var IDXpath
        # A list of idx points representing the UAV's path
        # through the nodes.
        self.IDXpath = list()
        ## @var GPSpath
        # A list of GPS points representing the UAV's path
        # through the map.
        self.GPSpath = list()
        ## @var score
        # Holds the score of the UAV's path
        self.score = 0
        ## @var score_mask
        # Holds a map of values which degrade the
        # score of the corresponding node after it has been visited.
        self.score_mask = None


class PathPlanner:

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
        ## @var maskSUB
        #   sets a subtractive penalty value to be applied to a node each
        #   time it is visited (set to 0 to apply no subtractive penalty)
        maskSUB = 0
        ## @var maskMUL
        #   sets a multiplicative penalty value to be applied to a node each
        #   time it is visited (set to 1 to apply no multiplicative penalty)
        maskMUL = 0.5
        ## @var maskMULthenSUB
        #   boolean which decides the order maskSUB and maskMUL are carried out
        maskMULthenSUB = True

    ## class constructor
    def __init__(self, vman):
        self.vman = vman
        plane = int(self.vman.flowfield.grid_resolution.z * self.vman.params.percent_height)
        self.X_map = [[_X_map_node() for x in range(self.vman.grid_res[0])]
                 for y in range(self.vman.grid_res[1])]
        self.X = deepcopy(self.vman.flowfield.x[:, :, plane])
        self.Y = deepcopy(self.vman.flowfield.y[:, :, plane])
        self.sens_mat = vman.calcSensitivityMatrix()
        self.d0 = self.vman.params.d0
        self.v0 = self.vman.params.v0
        self._calcScoreMap()
        self.update_cost_map();

    def _shiftVals(self, val):
        value = {
            dir.E: [1,0],
            dir.NE: [1,1],
            dir.N: [0,1],
            dir.NW: [-1,1],
            dir.W: [-1,0],
            dir.SW: [-1,-1],
            dir.S: [0,-1],
            dir.SE: [1,-1]
        }
        return value.get(val)

    def _euclidDist(self, P1, P2):
        # print(math.sqrt(((P2[0]-P1[0])**2)+((P2[1]-P1[1])**2)))
        return math.sqrt(((P2[0]-P1[0])**2)+((P2[1]-P1[1])**2))

    def update_cost_map(self):

        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                self.X_map[i][j].idx = [i,j]
                self.X_map[i][j].GPS = [self.X[i][j],self.Y[i][j]]
                self.X_map[i][j].score = self.score_map[i][j]
                for k in range(8):
                    self.X_map[i][j].Xitions[k].thetaCost = k*np.deg2rad(45)-self.d0
                    if self.X_map[i][j].Xitions[k].thetaCost > PI:
                        self.X_map[i][j].Xitions[k].thetaCost = abs(k*np.deg2rad(45)-self.d0-2*PI)
                    [I,J] = self._shiftVals(k)


                    try:
                        dS = self.score_map[i + I][j + J] - self.score_map[i][j]
                        if math.isnan(dS):
                            self.X_map[i][j].Xitions[k].dSscore = None

                        else:
                            self.X_map[i][j].Xitions[k].dSscore = deepcopy(dS)
                            # print(self.X_map[i][j].Xitions[k].dSscore)
                            # print(str(i)+' '+str(j)+' '+str(I)+' '+str(J))
                            # print(str(self.score_map[i+I][j+J])+' '+str(self.score_map[i][j]))
                    except:
                        self.X_map[i][j].Xitions[k].dSscore = None

    def _calcScoreMap(self):

        Zd = abs(self.sens_mat[1]) / np.amax(abs(self.sens_mat[1]))
        # print(Zd)
        Zv = abs(self.sens_mat[0]) / np.amax(abs(self.sens_mat[0]))
        # print(Zv)
        self.score_map = Zd * Zv
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                for coord, turbine in self.vman.flowfield.turbine_map.items():
                    # print(Pts[0])
                    if self._euclidDist([self.X[i][j],
                                         self.Y[i][j]],
                                        [coord.x, coord.y]) < \
                            self.params.minDist2turbine:
                        self.score_map[i][j] = None
        self._calcWaveMap()

    def _calcWaveMap(self, UAV=None):
        max = -100
        self.wave_map = [[1 for i in range(self.vman.grid_res[0])] \
                          for j in range(self.vman.grid_res[1])]
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                try:
                    if self.score_map[i][j]*UAV.score_mask[i][j]>max:
                        max=self.score_map[i][j]
                except:
                    if self.score_map[i][j]>max:
                        max=self.score_map[i][j]
        spots = np.where(self.score_map==max)
        x = spots[0][0]
        y = spots[1][0]
        for i in range(self.vman.grid_res[0]):
            for j in range(self.vman.grid_res[1]):
                self.wave_map[i][j] = self._euclidDist([x,y],[i,j])
        self.wave_map = 1-self.wave_map/np.amax(self.wave_map)

    def show(self):
        plt.show()

    def plotWaveMap(self, ax=None):
        ax = plt.scatter(x=self.X, y=self.Y, c=self.wave_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
        plt.colorbar(ax)
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')

    def plotScoreMap(self, ax=None):
        ax = plt.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
        plt.colorbar(ax)
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')

    def plotScoreMask(self, UAV, ax=None):
        ax = plt.scatter(x=self.X, y=self.Y, c=UAV.score_mask,
                         vmin=0, vmax=1)
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='lime')

    def plotScoreMapUAV(self, UAV, ax=None):
        ax = plt.scatter(x=self.X, y=self.Y, c=self.score_map,
                         vmin=0, vmax=1, cmap='gnuplot2')
        plt.colorbar(ax)
        plt.suptitle("Greedy UAV path\n"+str(len(UAV.GPSpath)-1)+" moves")
        for coord, turbine in self.vman.flowfield.turbine_map.items():
            a = Coordinate(coord.x, coord.y - turbine.rotor_radius)
            b = Coordinate(coord.x, coord.y + turbine.rotor_radius)
            a.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            b.rotate_z(turbine.yaw_angle - self.vman.flowfield.wind_direction, coord.as_tuple())
            plt.plot([a.xprime, b.xprime], [a.yprime, b.yprime], 'k', linewidth=1, color='black')

        plt.plot([i[0] for i in UAV.GPSpath],
                 [i[1] for i in UAV.GPSpath], color='lime')
        plt.subplots_adjust(left=0.07,
                            bottom=0.08,
                            right=1.0,
                            top=0.91,
                            wspace=0.27,
                            hspace=0.19)

    def _findGPSindex(self, GPS):
        # Figure out how 'wide' each range is
        XleftSpan = np.amax(self.X) - np.amin(self.X)
        XrightSpan = self.vman.grid_res[0]
        YleftSpan = np.amax(self.Y) - np.amin(self.Y)
        YrightSpan = self.vman.grid_res[1]

        # Convert the left range into a 0-1 range (float)
        idx = int((GPS[0] -  np.amin(self.X))* XrightSpan / float(XleftSpan))
        idy = int((GPS[1] -  np.amin(self.Y))* YrightSpan / float(YleftSpan))
        if math.isnan(self.score_map[idx][idy]) or idx<0 or idy<0:
            print("Error: GPS point is out of bounds")
            return [None, None]
        else:
            return [idx, idy]

    def _getGPSfromIDX(self, IDX):
        # print(IDX)
        X = self.X[IDX[0]][IDX[1]]
        Y = self.Y[IDX[0]][IDX[1]]
        return [X,Y]

    def greedyPath(self, UAV, steps=100):
        UAV.GPSpath.clear()
        UAV.GPSpath.append(UAV.GPS)
        UAV.idx = self._findGPSindex(UAV.GPS)
        UAV.IDXpath.clear()
        UAV.IDXpath.append(UAV.idx)
        UAV.score=0
        self._calcScoreMap()
        self._calcWaveMap()
        # print(UAV.idx)
        if (UAV.idx[0] and UAV.idx[1]):
            UAV.score_mask = [ [ 1 for i in range(len(self.score_map[0]))] \
                             for j in range(len(self.score_map))]
            # print(UAV.idx)
            for i in range(steps):
                x=UAV.idx[0]
                # print(x)
                y=UAV.idx[1]
                # print(y)

                # print("iteration:"+str(i))
                try:
                    if(self.params.maskMULthenSUB):
                        UAV.score_mask[x][y] *= self.params.maskMUL
                        UAV.score_mask[x][y] -= self.params.maskSUB
                    else:
                        UAV.score_mask[x][y] -= self.params.maskSUB
                        UAV.score_mask[x][y] *= self.params.maskMUL
                    self._calcWaveMap(UAV)
                    max, dir = self._greedyStep(UAV)
                except:
                    return
                UAV.score += max
                UAV.heading = dir
                [I, J] = self._shiftVals(dir)
                UAV.idx = [x+I,y+J]
                UAV.IDXpath.append(UAV.idx)
                UAV.GPS = self._getGPSfromIDX(UAV.idx)
                UAV.GPSpath.append(UAV.GPS)

        else:
            print("Aborting.")

    def _greedyStep(self, UAV):
        x = UAV.idx[0]
        y = UAV.idx[1]
        if x>=self.vman.grid_res[0] or \
            y>=self.vman.grid_res[1]:
            return None
        max = -10000000
        ind = 0
        try:
            Xition = self.X_map[x][y].Xitions
        except:
            return None
        for i in range(8):
            [I, J] = self._shiftVals(i)
            try:
                score = self.score_map[x+I][y+J]
                mask = UAV.score_mask[x+I][y+J]
                wave = self.wave_map[x+I][y+J]
                # print(score)
                cost = (self.params.scoreWt*score + \
                       self.params.dSwt*Xition[i].dSscore - \
                       self.params.windWt*Xition[i].thetaCost/PI + \
                       self.params.waveWt*wave + \
                       self.params.timeWt*len(UAV.IDXpath)) * \
                       mask

                if cost > max and \
                    x+I<self.vman.grid_res[0] and \
                    y+J<self.vman.grid_res[1] and \
                    x+I>=0 and y+J>=0:
                    max = cost
                    ind = i
            except:
                pass
        return max, ind
