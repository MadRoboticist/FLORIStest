## \file UAV.py
# contains the definition of a UAV
from copy import deepcopy

## \class UAV
# @brief A class for representing the state of a UAV on the windfarm
class UAV:

    ## Class constructor.
    def __init__(self):
        '''Class constructor.'''

        ''' This section holds "ACTUAL" values '''
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
        self.heading = 0
        ## @var dS
        # holds the current dS value for the UAV
        self.dS = 0
        ## @var path_mask
        # Holds a map of values which degrade the
        # score of the corresponding node after it
        # has been visited by the UAV.
        # (applied as calculated_score*mask_value)
        self.path_mask = None

        ''' This section contains PLANNING values '''
        ## @var plan_heading
        # A list of the enumerated directions that the UAV
        # plans to travel. see Class dir()
        self.plan_heading = list()
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
        ## @var dS
        # holds the previous values of dS for calculation of d2S during planning
        self.dSplan = list()
        ## @var score
        # Holds the score of the UAV's path
        self.score = 0
        ## @var wave_map
        # Holds a map of values which can serve as
        # a potential field, drawing the UAV to
        # the node with the maximum score value
        self.wave_map = None
        ## @var plan_mask
        # Holds a map of values which degrade the
        # score of the corresponding node after it
        # has been visited during path planning.
        # (applied as calculated_score*mask_value)
        self.plan_mask = None
        ## @var maskSUB
        #   sets a subtractive penalty value to be applied to a node's
        #   score mask each time it is visited
        #   (set to 0 to apply no subtractive penalty)

        ''' This section contains WEIGHTS '''
        self.maskSUB = 0
        ## @var maskMUL
        #   sets a multiplicative penalty value to be applied to a node's
        #   score mask each time it is visited
        #   (set to 1 to apply no multiplicative penalty)
        self.maskMUL = 0.5
        ## @var maskMULthenSUB
        #   boolean which decides the order maskSUB and maskMUL are carried out
        self.maskMULthenSUB = True
        ## @var scoreWt
        #   sets the weight of the sensitivity score in path calculations
        self.scoreWt = 10
        ## @var dSwt
        #   sets the weight of the derivative of the sensitivity score
        #   wrt transition in path calculations
        self.dSwt = 2
        ## @var d2Swt
        #   sets the weight of the 2nd derivative of the sensitivity score
        #   wrt transition in path calculations
        self.d2Swt = 2
        ## @var windWt
        #   sets the weight of the wind direction vs. transition heading
        #   in path calculations
        self.windWt = 2
        ## @var headWt
        #   sets the weight of the UAV heading vs. transition heading
        #   in path calculations
        self.headWt = 2
        ## @var waveWt
        #   sets the weight of the wave map in path calculations
        self.waveWt = 35
        ## @var timeWt
        #   sets the weight of time/iterations elapsed in path calculations
        self.timeWt = 0.25

    # update mask
    ## a method to update either the plan or the actual mask
    def update_mask(self, mask, IDX):
        x = IDX[0]
        y = IDX[1]
        if (self.maskMULthenSUB):
            # multiply then subtract
            mask[x][y] *= self.maskMUL
            mask[x][y] -= self.maskSUB
        else:
            # subtract then multiply
            mask[x][y] -= self.maskSUB
            mask[x][y] *= self.maskMUL

    def reset_planner(self):
        # clear the heading list
        self.plan_heading.clear()
        # append the current heading
        self.plan_heading.append(self.heading)
        # start with a new GPS path
        self.GPSplan.clear()
        # append the starting GPS point
        self.GPSplan.append(self.GPS)
        # also start with a new index path
        self.IDXplan.clear()
        # start with a new list of dS values
        self.dSplan.clear()
        # append the current dS
        self.dSplan.append(self.dS)
        # start the score at zero
        self.score = 0
        # clear the plan mask
        self.plan_mask = None

    # move
    ## a method to advance the UAV through the current planned path
    #   @param step
    def move(self, steps):
        # just in case, make a backup copy of the current values
        tempIDXpath = deepcopy(self.IDXpath)
        tempGPSpath = deepcopy(self.GPSpath)
        tempdS = deepcopy(UAV.dS)
        tempHead = deepcopy(UAV.heading)
        # step through the current plan
        for i in range(steps):
            # if the moves haven't been planned yet,
            # an exception will be thrown
            try:
                self.IDXpath.append(self.IDXplan[i])
                self.GPSpath.append(self.GPSplan[i])
                self.IDX = self.IDXpath[len(self.IDXpath) - 1]
                self.GPS = self.GPSpath[len(self.GPSpath) - 1]
                self.dS = self.dSplan[i]
                self.heading = self.plan_heading[i]
            except:
                # set everything back to how it was
                self.IDXpath.clear()
                for nodes in tempIDXpath:
                    self.IDXpath.append(nodes)
                self.GPSpath.clear()
                for nodes in tempGPSpath:
                    self.GPSpath.append(nodes)
                self.dS = tempdS
                self.heading = tempHead
                print("Error: This number of steps has not been planned.")
                return

