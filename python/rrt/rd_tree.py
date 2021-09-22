import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import time


#Import a few utility functions...
import sys  
from pathlib import Path
#sys.path.insert(0, "../comm_channel")
#sys.path.insert(0, "../polling_systems")
sys.path.insert(0, "../utils")

from pbin import PBTree 
from tree_node import TreeNode

Inf = float('inf')
class RDTSolver():
    """
    %RDTSolver - can solve instances of PathPlanningProblem using Randomly
    %Exploring Random Trees, including and RRT*. Flexible enough to
    %handle some variants 
    """
    def __init__(self, sampler, stop_criteria, do_rewire, steer_rad):
        #sampler - sampling module to use
        self.sampler = sampler

        #stopCriteria - function handle. Function determines whether or not
        #we should stop. Arguments are elapsed time, iteration number, and cost
        #of best solution so far
        self.stopCriteria = stop_criteria

        #rdTree - The tree being built out
        self.rdTree = RDTree(steer_rad, do_rewire)

    def solve(self, pppi):
        start_time = time.time()
        self._initializeTree(pppi)
        if not self.startInDest(pppi):
            self.sampler.setParams(pppi)

            iteration_count = 0
            series_delta = 0.2
            elapsed_time = time.time()-start_time
            time_for_next_recording = series_delta

            while not self.stopCriteria.stop(not np.isinf(self.getBSF().costToHere), 
                iteration_count, elapsed_time):
                x_rand = self.sampler.sample(iteration_count)
                #print(x_rand)
                improved = self.rdTree.executeIteration(x_rand, pppi)
                iteration_count += 1
                elapsed_time = time.time() - start_time

                if elapsed_time > time_for_next_recording:
                    time_for_next_recording += series_delta
                    self.rdTree.recordBSFCost(elapsed_time)

                if iteration_count %5== 0:
                    print(iteration_count)
                if (iteration_count in [100, 1000, 2000, 5000]) or improved:
                    print('Iteration count: %d'%(iteration_count))
                    n_nodes = self.rdTree.treeNodes.count
                    rut = self.getRoot()
                    rut.plotTree('k')
                    bsf = self.getBSF()
                    if not np.isinf(bsf.costToHere):
                        bsf.pathToRoot(do_plot=True)
                    pppi.plotProb()
                    plt.show()

        
    def getRoot(self):
        return self.rdTree.root 
            
    def getBSF(self):
        return self.rdTree.BSF

        
    def getBSFTimeSeries(self):
        return self.rdTree.solCosts

    def minCostSoFar(self):
        return self.getBSF().costToHere

    #Private methods

    def _initializeTree(self, pppi):
        #see Karaman & Frazzoli, 2011, page 29. Here, d=2 
        root_pos = pppi.source
        region = pppi.region
        lb_ms = 1
        d = len(region)//2
        lb_unit_ball = (np.pi**(d/2))/_euler_gamma(1+d/2)
        for i in range(d):
            lb_ms *=(region[(i*2)+1] - region[i*2])
        gamma = (2*(1+1/d)*(lb_ms/lb_unit_ball))**(1/d) + 1#need to make it greater

        self.rdTree.initialize(gamma, root_pos, pppi)

    def startInDest(self, pppi):
        root_pos = pppi.source
        in_dest = pppi.goalRegion.contains_point(root_pos)
        if in_dest:
            self.rdTree.BSF = self.getRoot()
        return in_dest
    

def _euler_gamma(t):
    term = t-1
    prod = 1
    while term >=1:
        prod *= term
        term -=1
    if term == 0.5:
        prod *=0.5*np.sqrt(np.pi)
    return prod

def steer(steer_rad, dest, start, pppi):
    diff = dest - start
    #Don't extend beyond the newly sampled point
    nrm = np.linalg.norm(diff)
    if nrm > steer_rad:
        new_diff = steer_rad*diff/nrm
        new = start + new_diff
        #TODO - do we need to snap to grid?
    else:
        new = dest

    return np.array([start, new])


class RDTree():

    def __init__(self, steer_rad, do_rewire=True):

        #doRewire - if true, will rewire, that is, will run RRT* rather
        #than RRT
        self.doRewire = do_rewire

        #steerRad - maximum distance the tree can grow at a time, that is
        #maximum distance between p_nearest and p_new
        self.steerRad = steer_rad

        #BSF = TreeNode. Path from BSF to root traces the best path found
        #so far
        self.BSF = None

        #treeNodes = [root, n1, n2, ....], the list of nodes in the tree,
        #in order in which they are added.
        self.treeNodes = None
        self.root = None

        #see "Sampling-based Algorithms for Optimal Motion Planning",
        #s. Karaman, E. Frazzoli, 2011 
        self.gamma= -1

        self.solCosts = []

        #neighborBoxes - data structure presented in "Minimising 
        #computational complexity of the RRT algorithm a practical
        #approach" M Sventrup et. al, 2011
        #neighborBoxes;

        self.theta = -1
    
        #setup asspects of the tree before beginning the algorithm
    def initialize(self, gamma, root_pos, pppi, theta=1):
        self.theta = theta
        self.gamma = gamma

        self.root = TreeNode(root_pos)
        self.root.isRoot = 1
        self.root.costToHere = 0
        self.treeNodes = PBTree(pppi.region)
        self.treeNodes.addTNode(self.root)

        #create a placeholder BSF with infinite cost
        self.BSF = TreeNode([-Inf,-Inf])
        self.BSF.costToHere = Inf

        
    def executeIteration(self, x_rand, pppi):
        improved = False
        success, new_node = self.tryAddNode(x_rand, pppi)
        # print(new_node.x)
        # print(success)
        #check if the new_node has given us a new BSF
        if success and pppi.nodeInDest(new_node) and self.BSF.costToHere > new_node.costToHere: 
            self.BSF = new_node
            improved = True
        if not success:
            print('new node not added')
        return improved
        
    def tryAddNode(self, x_rand, pppi):
        success = False
        #initialize dummy new node
        n_new = TreeNode([float('inf'), float('inf')])
        nearest, min_dist = self.nearest(x_rand)
        #check to see if we already have this node added to the tree
        #should be able to just check if zero again, but let's validate
        if  min_dist != 0:
            print('sampled new')
            path = self.steer(x_rand, nearest, pppi);
            #check if the path is collision free, trimming if necessary
            viable_path = pppi.collisionFree(path)
            if len(viable_path) > 1: # more than one point in path
                success = True
                n_new = self.addNode(nearest, viable_path, pppi)
            else:
                print(path)
                print(viable_path)
        else:
            print('sampled duplicate')
            #handle the case where we resample an already sampled
            #point. Need to rewire. Look into literature
            self.rewire(nearest, pppi)
        return success, n_new

     
    def recordBSFCost(self, time):
        if not np.isinf(self.BSF.costToHere):
            self.solCosts.append((time, self.BSF.costToHere))

    def addNode(self, nearest, path, pppi):
        #add a node to the tree
        new_node =TreeNode(path[-1, :])
        #print(path)
        cost = pppi.pathCost(nearest, new_node)
        new_node.setParent(nearest, cost, path)
        if self.doRewire:
            self.rewire(new_node, pppi)
        self.treeNodes.addTNode(new_node)

        return new_node

    """
    rewire
        % TODO - rewire also includes finding shortest path (not just
        % closest TO the new node)
        % TODO - implement with balanced box decomposition
        %can be implemented in O(log n) time using Balanced Box Decomposition
        %currently is O(n) (brute force)
    """
    def rewire(self, new_node, pppi):
            #See Karaman & Frazzoli, 2011
            cardV = self.treeNodes.count
            radius = min(self.steerRad, self.gamma*np.sqrt(np.log(cardV)/cardV))
            neighbors = self.near(new_node.x, radius)

            cost_btwn = -1*np.ones(len(neighbors))
            x_current = new_node.x
            if new_node.isRoot:
                min_cost = 0
                min_cost_btwn = 0
            else:
                best_neighbor = new_node.parent;
                min_cost = new_node.costToHere;
                min_cost_btwn = min_cost - best_neighbor.costToHere;
                best_path = new_node.pathFromParent;
            
            #first, find the least cost path to current
            for i in range(len(neighbors)):
                neighbor = neighbors[i]
                if (not new_node.isRoot) and (neighbor == new_node.parent):
                    continue;
                #want path from from  neighbor to current
                path = np.array([neighbor.x, x_current])
                
                #check to make sure there's nothing obstructing
                viable_path = pppi.collisionFree(path);
                if len(viable_path) == len(path):
                    #full path works!
                    cost_btwn[i] = self.theta*pppi.pathCost(new_node, neighbor)
                    cost_to_new_via_neighbor = cost_btwn[i] + neighbor.costToHere
                    if cost_to_new_via_neighbor < min_cost:            
                           best_path = path;
                           best_neighbor = neighbor;
                           min_cost = cost_to_new_via_neighbor;
                           min_cost_btwn = cost_btwn[i]
            
            if not new_node.isRoot:
                new_node.setParent(best_neighbor, min_cost_btwn, best_path); 
           
            #now check if there are any nodes that we can reach with less
            #cost from our new_node
            for j in range(len(neighbors)):
               neighbor = neighbors[j]
               cost = cost_btwn[j]
               if cost >= 0:
                    if (neighbor.costToHere > (new_node.costToHere + cost)): 
                       #we've already seen that there's nothing
                       #obstructing from the first time we looped through
                       #path from new to neighbor
                       path = np.array([x_current, neighbor.x])
                       #now add the new parent
                       neighbor.setParent(new_node, cost, path)

    def steer(self, x_rand, v_nearest, pppi):
        #do a little massaging to increase sample efficiency
        if v_nearest.x[4] < x_rand[4]:
            x_rand[4] = v_nearest.x[4]
        print(v_nearest.x)
        print(x_rand)
        return steer(self.steerRad, x_rand, v_nearest.x, pppi)
        
    def near(self, pt, radius):
        return self.treeNodes.findNeighborhood(pt, radius)

    def nearest(self, x_rand):
        return self.treeNodes.findNearest(x_rand)
