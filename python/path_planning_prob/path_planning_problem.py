import numpy as np
from matplotlib import pyplot as plt

"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PathPlanningProblem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Represents an instance of a path planning problem to be solved. Handles 
% discretization over continuous space as well as obstacle checking and 
% performance metric, as all these should remain the same regardless of how
% the problem is solved.
"""

class PathPlanningProblem():

    def __init__(self, region, source, goal_region, obstacle_mod, cost_func):
        """
        %ObstacleModule - abstraction of obstacles. Must implement
        %CollisionFree(p1,p2)->{T,F}. Checks to see if the path from p1 to
        %p2 is collision free. May also implement a plot() method
        """
        self.obstacleMod = obstacle_mod

        #region = [d1_max d1_min d2_max d2_min ... dn_max dn_min]; Raw coordinates from problem
        self.region = region

        #source = [x_start, y_start]; Raw coordinates
        self.source = source
        
        #goalRegion - GoalRegion object. The set of points to which the
        #robot is heading. 
        self.goalRegion = goal_region
                
        ##costFunction - maps two points (n1, n2) to the cost between
        #them. n1, n2 are TreeNode objects.
        self.costFunction = cost_func

    """
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pathCost
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Creates an identical copy of the PathPlanningProblem but with a
        % different source
        % Input:
        % self - reference to the PathPlanningProblem object
        % n1 - starting node (TreeNode)
        % n2 - ending node (TreeNode)
        % path - path between n1 and n2
        %
        % Output:
        % cost of moving from n1 to n2
    """
    def pathCost(self, n1, n2):
        return self.costFunction(n1, n2)

    """
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % nodeInDest
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cheks if node is in the goal region
        % Input:
        % self - reference to the PathPlanningProblem object
        % n1 - node (TreeNode) to check
        %
        % Output:
        % true if the node is in the goal region. False otherwise
    """
    def nodeInDest(self, n1):
        return self.goalRegion.contains_point(n1.x)

    """
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % collisionFree
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cheks if path is collision free and trims if necessary
        % function.
        % Input:
        % self - reference to the PathPlanningProblem object
        % path - series of path waypoints
        %
        % Output:
        % viable_subpath - returns the portion of the path that is not
        %                   obstructed. For example if there is a collision
        %                   at path(i), will return path(1:i-1)
    """
    def collisionFree(self, path):
        terminate_early = False
        #TODO chop path into little pieces so that we can salvage some of the path
        p_len = len(path)
        for i in range(p_len-1):
            if self.goalRegion.contains_point(path[i]):
                terminate_early = True
                break
            p1 = path[i]
            p2 = path[i+1]
            if not self.obstacleMod.collisionFree(p1, p2):
                terminate_early = True
                break

        viable_subpath = path
        if terminate_early:
            viable_subpath = path[:i+1,:]
        return viable_subpath

    """
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ptInRegion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cheks if point is in the workspace
        % Input:
        % self - reference to the PathPlanningProblem object
        % pt - point [x,y] 
        %
        % Output:
        % tf - true if the point is in the workspace. False otherwise
    """
    def ptInRegion(self, pt):

        d_reg = len(self.region)//2
        d_pt = len(pt)
        if d_pt != d_reg:
            return False

        in_reg = True
        for i in range(d_reg):
            d_max = self.region[2*i]; d_min = self.region[2*i+1] 
            in_reg = in_reg and (d_min<pt[i]<d_max)
            if not in_reg:
                break
        return in_reg

    """
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotProb
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots the problem
        % Input:
        % self - reference to the PathPlanningProblem object
        %
    """
    def plotProb(self):
        self.obstacleMod.plotObstacles()
        plt.scatter(self.source[0], self.source[1])
        dests = self.goalRegion.getPoints()
        plt.scatter(dests[:,0], dests[:,1])
        plt.xlim(self.region[1], self.region[0])
        plt.ylim(self.region[3], self.region[2])