import numpy as np
from tree_node import TreeNode

#Import a few utility functions...
import sys  
from pathlib import Path
sys.path.insert(0, "../geometry")
from rrange import RRange

class PBTree():

	class PBin():
		P_MAX = 100

		def __init__(self, ranges, d_split):
			self.ranges = ranges
			self.d = len(ranges)
			self.d_split = d_split
			self.left = None
			self.right = None
			self.tnodes = []
			self.isLeaf = True

		def __str__(self):
			out_str = 'isLeaf: ' + str(self.isLeaf)
			if len(self.tnodes) >0:
				out_str+=str([n.x for n in self.tnodes])
			return out_str


		def containsPoint(self, point):
			contains = True
			for i in range(self.d):
				if not (self.ranges[i].includes(point[i])):
					contains = False
					break
			return contains

		def addTNode(self, tnode):
			if self.isLeaf:
				self.tnodes.append(tnode)
				if len(self.tnodes) > PBTree.PBin.P_MAX:
					self.split()
			elif self.left.containsPoint(tnode.x):
				self.left.addTNode(tnode)
			elif self.right.containsPoint(tnode.x):
				self.right.addTNode(tnode)
			else:
				print('Point out of range.')				

		def split(self):
			next_d_split = (self.d_split + 1) % self.d
			split_range = self.ranges[self.d_split]
			l_min = split_range.min
			#always split in half by range values, since we're assuming uniform sampling
			l_max = (split_range.max - split_range.min)/2 + l_min
			r_min = l_max
			r_max = split_range.max
			#build left pbin
			l_ranges = self.ranges
			l_ranges[self.d_split] = RRange(l_max, l_min)
			self.left = PBTree.PBin(l_ranges, next_d_split)
			#build right 
			r_ranges = self.ranges

			r_ranges[self.d_split] = RRange(r_max, r_min)
			self.right = PBTree.PBin(r_ranges, next_d_split)

			for tn in self.tnodes:
				if self.left.containsPoint(tn.x):
					self.left.addTNode(tn)
				else:
					self.right.addTNode(tn)

			self.tnodes = []
			self.isLeaf = False

		def nearestInBin(self, point):
			min_dist = float('inf')
			nearest = None
			for neighbor in self.tnodes:
				dist = np.linalg.norm(point - neighbor.x)
				if dist < min_dist:
					min_dist = dist
					nearest = neighbor
			return min_dist, nearest


	def __init__(self, region, p_max = 100):
		self.d = len(region)//2#number of dimensions
		self.bounds = []
		self.p_max = p_max #number of elements allowed to be in the same bin
		for i in range(self.d):
			self.bounds.append(RRange(region[2*i], region[2*i + 1]))

		self.root = PBTree.PBin(self.bounds, 0)
		self.count = 0

	def addTNode(self, tnode):
		self.root.addTNode(tnode)
		self.count += 1

	def _getContainingLeafPBin(self, point):
		pbin = self.root
		while not pbin.isLeaf:
			if pbin.left.containsPoint(point):
				pbin = pbin.left
			elif pbin.right.containsPoint(point):
				pbin = pbin.right
		return pbin

	def findNearest(self, point):
		pbin = self._getContainingLeafPBin(point)
		min_dist, nearest = pbin.nearestInBin(point)
		#now check to see if it's possible the nearest actually lies in a neighboring bin
		#TODO - this only checks points orthogonal, would want to check all
		neighbor_pbins = self._getNeighborhoodPBins(pbin, point, min_dist)
		for pbin in neighbor_pbins:
			test_min, test_neighbor = pbin.nearestInBin(point)
			if test_min < min_dist:
				min_dist = test_min
				nearest = test_neighbor
		return nearest, min_dist

	def findNeighborhood(self, point, radius):
		neighbors = []
		#only search neighboring, should be sufficient
		pbin = self._getContainingLeafPBin(point)
		neighbor_pbins = self._getNeighborhoodPBins(pbin, point, radius)
		all_pbins = [pbin, *neighbor_pbins]
		for pbin in all_pbins:
			for tn in pbin.tnodes:
				if np.linalg.norm(point - tn.x) <= radius:
					neighbors.append(tn)

		return neighbors

	def _getNeighborhoodPBins(self, pbin, point, max_dist):
		neighbor_pbins = []
		for i in range(self.d):
			if point[i] - pbin.ranges[i].min < max_dist:
				#check over that way
				test_point = [x for x in point]#create a deep copy
				test_point[i] = max(pbin.ranges[i].min - 0.1, self.root.ranges[i].min)
				pbin = self._getContainingLeafPBin(test_point)
				if not pbin in neighbor_pbins:
					neighbor_pbins.append(pbin)

			if pbin.ranges[i].max - point[i] < max_dist:
				#check over this way
				test_point = [x for x in point]#create a deep copy
				test_point[i] = min(pbin.ranges[i].min + 0.1, self.root.ranges[i].max)
				pbin = self._getContainingLeafPBin(test_point)
				if not pbin in neighbor_pbins:
					neighbor_pbins.append(pbin)

		return neighbor_pbins

