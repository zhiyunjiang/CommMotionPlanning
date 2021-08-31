import numpy as np
from matplotlib import pyplot as plt

#Import a few utility functions...
import sys  
from pathlib import Path
#sys.path.insert(0, "../comm_channel")
#sys.path.insert(0, "../polling_systems")
sys.path.insert(0, "../geometry")
#sys.path.insert(0, "../utils")
import pointcloud as PC

"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ObstacleMod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Containing class for the collection of obstacles the workspace of our
% path palnning problem

	% Properties (public):
    % obstacles - list of all obstacles
    
    % Methods (public):
    % (constructor) - create new ObstacleMod object
    % collisionFree - checks if the straightline path form p1 to p2 is
    %                   obstacle free
    % plotObstacles - plots all obstacles
"""
class ObstacleMod():

	def __init__(self, obstacles):
		self.obstacles = obstacles
	"""
	 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % collisionFree
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Checks if the straightline path from p1 to p2 is obstacle free
        %
        % Inputs:
        % this - reference to this ObstacleMod object
        % p1 - first waypoint
        % p2 - second waypoiny
        %
        % Output:
        % True if the path is collision free, False otherwise
	"""
	def collisionFree(self, p1, p2):
		c_free = True
		for ob in self.obstacles:
			c_free = not ob.obstructs(p1, p2)
			if not c_free:
				break
		return c_free

	"""
	 	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotObstacles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots all obstacles
        %
        % Inputs:
        % self - reference to this ObstacleMod object
	"""
	def plotObstacles(self):
		for ob in self.obstacles:
			ob.plot()

class Obstacle():

	def __init__(self):
		pass

	def obstructs(self, p1, p2):
		print('To be overwritten by implementing classes')

	def plot(self):
		print('To be overwritten by implementing classes')

class CircleObstacle(Obstacle):

	def __init__(self, c, r):
		self.c = c
		self.r = r

	"""
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % obstructs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % checks if the CircObs obstructs the path between two
        % waypoints.
        % Input:
        % self - reference to the CircObs object
        % p1 - first waypoint [x,y]
        % p2 - second waypoint [x,y]
        %
        % Output:
        % does_obstruct - true if obstacle blocks path. False otherwise.
	"""
	def obstructs(self, p1, p2):
		#only look at first two coordinates (x and y)
		p1 = p1[:2]
		p2 = p2[:2]
		does_obstruct = True
		# calc the min dist to the line on which the line segment a-b lies
		a_to_c = self.c - p1
		a_to_b = p1-p2

		#find angle between ac and ab

		theta = np.real(np.arccos(max(min( (a_to_c @ a_to_b)/(np.linalg.norm(a_to_c)*np.linalg.norm(a_to_b)),1), -1)))
		min_dist = np.linalg.norm(a_to_c)*np.sin(theta % np.pi)
		#If the min distance from the circle center to the line is
		# greater than the circle radius, we're set
		if min_dist > self.r:
			return False
		else:
			#otherwise, check if min occurs within the line segment
			theta_c = np.arctan2(a_to_c[1], a_to_c[0])
			x_delta = abs(min_dist*np.cos(theta_c))*np.sign(p1[0]-self.c[0])
			x_min = x_delta + self.c[0]

			if (p1[0]>x_min and p2[0] > x_min) or (p1[0]<x_min and p2[0]<x_min):
				#if not, then just check the end points
				does_obstruct = np.linalg.norm(a_to_c) <= self.r or np.linalg.norm(self.c - p2) <= self.r
		return does_obstruct

		"""
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% plotObstacle
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% plots the obstacle
		% Input:
		% this - reference to the CircObs object
		"""
	def plot(self):
		th = np.linspace(0, 2*np.pi,6.5*50)
		x = r*np.cos(th)+self.c[0]
		y = r*np.sin(th)+self.c[1]
		#plt.plot(x, y)
		plt.fill(x, y)
     
class PolyObstacle(Obstacle):

	def __init__(self, vertices):
		self.poly = PC.Poly(vertices)

	def obstructs(self, p1, p2):
		does_obstruct = False
		for edge in self.poly.edges():
			if edge.intersects_segment(p1, p2):
				does_obstruct = True
				break
		return does_obstruct

	def plot(self):
		pts = self.poly.points
		plt.fill(pts[:,0],pts[:,1])

class AccelerationCap(Obstacle):

	def __init__(self, a_max):
		self.a_max = a_max

	def obstructs(self, p1, p2):
		dist = np.linalg.norm(p1[:2] - p2[:2])

		v1 = p1[2:4]
		v2 = p2[2:4]
		time = 2*dist/np.linalg.norm(v1+v2) #assume linear, constant acceleration
		a = np.linalg.norm(v1 - v2)/time
		return (a > self.a_max) or a<=0

	def plot(self):
		pass

class SpectralEfficiencyCap(Obstacle):

	def __init__(self, r_max, BW, channel=None, power_max = float('inf')):
		self.r_max = r_max
		self.BW = BW
		self.channel = channel
		self.power_max = power_max


	def obstructs(self, p1, p2):
		dist = np.linalg.norm(p1[:2] - p2[:2])

		v1 = p1[2:4]
		v2 = p2[2:4]
		time = 2*dist/np.linalg.norm(v1+v2) #assume linear, constant acceleration
		data = p1[4] - p2[4]#data has to be increasing
		if data < 0:
			return True
		r_avg = data/(self.BW*time)
		return r_avg > self.r_max
		#req_power = self.channel.getReqTXPowerAtPoint(pt, qos)

	def plot(self):
		pass

