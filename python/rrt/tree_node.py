import numpy as np
import matplotlib
from matplotlib import pyplot as plt
class TreeNode():
    """
    %TreeNode - represents a node in the tree created by a number of related
    % Random Dense Tree algorithms. Also includes convenient printing
    % methods
    %
    %PROPERTIES
    % x - state of the node in the workspace
    % parent - parent node. If this node is the root, then parent is itself
    % children - array of child nodes
    % isRoot - Flag set high for root node, low otherwise
    % costToHere - the distance to this node through the tree
    """
    def __init__(self, state):
        self.x = state
        self.parent = None
        self.pathFromParent = None
        self.children = []
        self.isRoot = False
        self.costToHere = float('inf')
        self.isPropagating = False

    def __str__(self):
        return 'x: '+str(self.x)+'\ncost to here: '+ str(self.costToHere) + '\nNum Kids: %d'%(len(self.children))
        
    def distTo(self, other, norm_num=2):
        return np.linalg.norm(self.x - other.x, norm_num)

    def setParent(self, parent, cost, path_from_parent=None):

        if self == parent:
            error_string = "Setting self as parent. Cost: %.2f\nParent: %s\nChild: %s"%(cost, str(parent), str(self))
            raise ValueError(error_string)
        if self in parent._getAncestry():#check the node, not the value (rounding issues?)
            print(parent.pathToRoot())
            error_string = "Cycle Detected. Cost: %.2f\nParent: %s\nChild: %s"%(cost, str(parent), str(self))
            raise ValueError(error_string)

        #need to update cost to other nodes down the line, too
        self._updateDist(parent.costToHere, cost)
        if self.parent is not None:
            self.parent.removeChild(self)

        self.parent = parent;
        parent.children.append(self)


        #if the path_to_parent is actually passed
        if path_from_parent is not None:
            self.pathFromParent = path_from_parent
        else:
            #otherwise assume it's just the line between them
            self.pathFromParent = np.array([parent.x, self.x])
        
    def setDist(self, new_dist):

        if new_dist < 0:
            error("Distance must always be non-negative");
        else:
            self.costToHere = new_dist;
          
    def getRootNode(self):
        node = self
        while not node.isRoot:
            node = node.parent 
        return node
         
    def __eq__(self, other):
        if type(other) != TreeNode:
            return False
        return (np.linalg.norm(self.x - other.x) == 0)
        
    def removeChild(self, thatChild):
        for i in range(len(self.children)):
            thisChild = self.children[i]
            if  thisChild == thatChild:
                del self.children[i]
                break


    def pathToRoot(self, do_plot=False):
        current = self
        path = []
        while not current.isRoot:
            #only take the path up to but not including the parent
            #parent is first, so don't add it (will be added next
            #iteration
            path = [ *current.pathFromParent[0:,:], *path]
            current = current.parent
        #tack on the root
        path = [current.x, *path]
        path = np.array(path)
        if do_plot:
            plt.plot(path[:,0], path[:,1])
        return path
        
    def plotTree(self,style):
        queue = [self];
        while len(queue) >0:
            #pop the first off the queue
            current = queue.pop(0)
            for child in current.children:
                path_to_parent = child.pathFromParent
                plt.plot(path_to_parent[:,0], path_to_parent[:,1], style)
                queue = [*queue, child];

#Private Functions

    def _getAncestry(self):
        node = self
        ancestry = [node]
        while not node.isRoot:
            node = node.parent
            ancestry.append(node)
        return ancestry
            
    def _updateDist(self, parent_dist, dist_btwn):
        if (np.isnan(parent_dist) or np.isinf(parent_dist) or parent_dist < 0 or
           np.isnan(dist_btwn) or np.isinf(dist_btwn) or dist_btwn < 0):
            error_string = "Parent: %s\nChild: %s\nDistance must be non-negative, finite number."%(str(self.parent), str(self))
            raise Exception(error_string)
        
        if np.isinf(self.costToHere):
            #we've just connected, we have no children!
            self.setDist(parent_dist + dist_btwn);
        else:
            diff = self.costToHere - (parent_dist + dist_btwn);
            self._propogateDiff(diff);
            
    #if the distance to me gets updated, need to update distance to my
    #children
    #TODO - update to take into account the problem instance's distance
    #metric. Should only need to be calculated for rewire
    def _propogateDiff(self, diff):
        if self.isPropagating:
            raise StopIteration('Cycle detected in tree.\n'+str(self))
        self.isPropagating = True
        #avoid some very smoll rounding issues
        self.setDist( max(0,self.costToHere - diff) )

        for child in self.children:
            child._propogateDiff(diff);
        self.isPropagating = False