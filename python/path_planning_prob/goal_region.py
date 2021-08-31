
"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GoalRegion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Represents the destination of our path planning problem. The
% destinatation region is taken to be a set of points from the simulated
% communication channel, and will thus generally has the same resolution as
% the communication channel simulation. This may differ from the problem
% resolution.
"""

class GoalRegion():

    def __init__(self, point_cloud):
        self.pc = point_cloud
        self.pc.partition()

    def contains_point(self, point):
        contains = False
        position = point[:2]
        for ply in self.pc.polygons:
            if ply.contains_point(point):
                contains = True
                if len(point) > 2:
                    v = point[2:4]
                    contains = ( v[0]!=1 or v[1]!=0)
                    if len(point) >=5:
                        contains = point[4]==0
            if contains:
                break
        return contains

    def getPoints(self):
        return self.pc.points

class PtGoalRegion():

    def __init__(self, points):
        self.pts = points

    def contains_point(self, point):
        return point in self.pts

    def getPoints():
        return self.pts
