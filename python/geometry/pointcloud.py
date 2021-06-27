#the usual suspects
import numpy as np
import matplotlib.pyplot as plt

#a few auxillary imports needed for cgal_partition
import sys
sys.path.insert(0, "../cython")#Needed for cgal partitioning

#CGAL imports
from CGAL import CGAL_Alpha_shape_2 as CGAL_alphashape
from CGAL.CGAL_Kernel import Point_2
import cgal_partition

#Ramer-Douglas-Peucker algo
from rdp import rdp

#triangulation that can handle holes
import triangle


class PointCloud():

	def __init__(self, points):
		self.points = points
		self.polygons = None
		
	def partition(self, alpha=0.01, algo = 0):
		if self.polygons is None:
			self._find_polys(alpha, algo)
		
		return self.polygons
		
	def plot(self, frmt = ''):
		plt.plot(self.points[:,0], self.points[:,1], frmt)
		
	def plot_polys(self, show_partition=True):
		if self.polygons is not None:
			for poly in self.polygons:
				poly.plot(show_partition)
		
	def _find_polys(self, alpha, algo):
		#find the polygons using alphashapes, then extract
		CGAL_pts = to_CGAL_points(self.points)
		shapes = CGAL_alphashape.Alpha_shape_2(CGAL_alphashape.REGULARIZED)
		shapes.make_alpha_shape(CGAL_pts)
		shapes.set_alpha(alpha)
		
		edges = PointCloud._edges_from_alphashape(shapes)
		print('Extracted edges...')
		self.polygons, incompletes = PointCloud._build_polys(edges)
		print('Polygons constructed, looking for holes...')
		#Find holes
		n = len(self.polygons)
		mask = []
		for i in range(n):
			polyi = self.polygons[i]
			was_hole = False
			for j in range(n):
				polyj = self.polygons[j]
				if (i != j) and not polyj.is_interior:

					interior_tf_vals = [polyj.contains_point(pt) for pt in polyi.points]
					if all(interior_tf_vals):
						polyj.add_hole(polyi)
						print('Found hole')
				if polyi.is_interior:
					break
			if not polyi.is_interior:
				mask.append(i)
				
		self.polygons = [self.polygons[m] for m in mask]
		
		print('Constructed %d possibly non-convex polygons'%(len(self.polygons)))
		#now decompose each polygon
		n_subregions = 0
		for poly in self.polygons:
			poly.make_cc()
			poly.reduce()
			poly.partition(algo = algo)
			n_subregions += len(poly.cnvx_partition)

		print('%d total subregions'%(n_subregions))
		
	def _build_polys(incompletes):
		completes = []
		ni = len(incompletes)
		terminated = (ni==0)
		while not terminated:
			temp_incompletes = []

			for i in range(len(incompletes)):
				if not incompletes[i].merged:#if it was merged earlier, just ignore
					for j in range(i+1,len(incompletes)):
			    			#try to merge
						if not incompletes[j].merged:
							incompletes[i].try_merge(incompletes[j])
						
							#If we've completed this one, stop working on it
							if incompletes[i].complete:
								completes.append(Poly(incompletes[i].points))
								break
			
					if not incompletes[i].complete:
						temp_incompletes.append(incompletes[i])

			incompletes = temp_incompletes

			if ni == len(incompletes):
				#check for self loops
				found_loop = False
				for inc in incompletes:
					loop = inc.check_for_loop()
					if loop is not None:
						completes.append(loop)
						found_loop = True
				if not found_loop:
					if ni > 0:
						print('No Improvement, orphaned edges exist')
						for inc in incompletes:
							print(inc)
					terminated = True
			else:
				ni = len(incompletes)

		return completes, incompletes
	
	def _edges_from_alphashape(cgala):
		cgal_edges = [cgala.segment(it) for it in cgala.alpha_shape_edges()]
		edges = [PolyConstructor(np.array([[edge.point(0).x(), edge.point(0).y()], [edge.point(1).x(), edge.point(1).y()]])) for edge in cgal_edges]
		#there is almost certainly a better way to do this
		all_edges_good = False
		all_points = []
		#list out all the points, including duplicates
		for edge in edges:
			all_points += edge.points.tolist()

		while not all_edges_good:
			good_edges = []
			n = len(edges)
			for edge in edges:
				if all_points.count(edge.points[0].tolist()) >= 2 and all_points.count(edge.points[1].tolist()) >= 2:
					good_edges.append(edge)
			
			if len(good_edges) == n:
				all_edges_good = True
		    
			edges = good_edges
			
		return good_edges
	

class PolyConstructor:

	def __init__(self, initial_verts):
		self.points = initial_verts
		self.complete = False
		self.merged = False
		
	def __str__(self):
		pstr = 'PolyConstructor:'
		for point in self.points:
			pstr += ' (' +str(point) + ') '
		return pstr
		
	def try_merge(self, other):
		if self.complete:
			return False

		did_merge = self._try_merge(other.points)

		if not did_merge:
			did_merge = self._try_merge(other.flip().points)

		if did_merge:
			other.merged = True
			

		return did_merge
	
	def check_for_loop(self):
	
		loop = None
		n = len(self.points)
		start = self.points[0]
		end = self.points[-1]
		for i in range(n):
			if (i != 0) and PolyConstructor.pequal(start, self.points[i]):

				loop = Poly(self.points[:i])
				self.points = self.points[i:]
				
				break
			if (i != n-1) and PolyConstructor.pequal(end, self.points[i]):

				loop = Poly(self.points[i+1:])	
				self.points = self.points[:i+1]
				break
				
		return loop

	def flip(self):
		self.points = np.flip(self.points, axis=0)
		return self
		


	
	def _try_merge(self, other_pts):
		did_merge = False
		if PolyConstructor.pequal( self.points[0], other_pts[-1]):
			#prepend
			self.points = np.concatenate( (other_pts[:-1], self.points), axis = 0)
			did_merge = True

		if PolyConstructor.pequal(self.points[-1], other_pts[0]):
			if did_merge:
				self.complete = True
				#and just trim off the last one, which is redundant, since we just added the new points
				self.points = self.points[:-1]				
			else:
				self.points = np.concatenate((self.points[:-1], other_pts), axis=0)
				did_merge = True

		return did_merge
		
	def pequal(p0, p1):
		return (p0[0] == p1[0]) and (p0[1] == p1[1])

	
		
class Poly:

	def __init__(self, points, holes = None):
		self.setPoints(points)
		if holes is None:
			holes = []
		self.holes = holes
		self.cnvx_partition = None
		self.is_interior = False
	
	def __str__(self):
		pstr = 'Poly:'
		for point in self.points:
			pstr += ' (' +str(point) + ') '
		pstr += 'Is Interior: ' + str(self.is_interior)
		return pstr

	#build the A matrix and b vector s.t. Ax-b < 0vec iff point in polygon		
	def to_linear_constraints(self):
		#only works if we're working with a convex polygon
		if not self.is_convex():
			return None
			
		#number of constraints equal to number edges
		edges = self.edges()
		ne = len(edges)
		A = np.zeros((ne, 2))
		b = np.zeros((ne, 1))
		
		pt = self.get_interior_point()
		for i in range(ne):
			Ai, bi = edges[i].as_constraint()
			if Ai @ pt - bi >0:
				Ai *= -1
				bi *= -1
			A[i] = Ai
			b[i] = bi
			
		return (A, b)
		
	def contains_point(self, pt):
		in_poly =  point_in_poly(self.points, pt)
		in_hole = False
		if in_poly:
			for hole in self.holes:
				in_hole = hole.contains_point(pt)
				if in_hole:
					break
		return (in_poly and not in_hole)
		
		
	def area(self):
		area = poly_area(self.points)
		for hole in self.holes:
			area -= hole.area()
		return area
		
	def is_convex(self):
		return len(self.get_reflex_vertices()) == 0
		
	def get_reflex_vertices(self):
		self.make_cc()
		n = len(self.points)
		reflex=[]
		for i in range(n):
			a = self.points[i]
			b = self.points[(i+1)%n]
			c = self.points[(i+2)%n]
			if not is_ccw(a,b,c):
				reflex.append(self.points[(i+1)%n])
				
		return reflex
		
	def add_hole(self, other):
		if self == other:
			print("WARNING: tried to add object to itself as a hole")
			print(other)
		else:
			self.holes.append(other)
			other.is_interior = True
		
	def partition(self, algo=0):
		
		#simplify holes before partitioning
		for hole in self.holes:
			hole.make_cc()
			hole.reduce()
		
		if algo == 0:
			points_list = cgal_partition.optmal_convex_parition(self.points)
			self.cnvx_partition = [Poly(points) for points in points_list]
		elif algo == 1:
			points_list = cgal_partition.greene_approx_convex_partition(self.points)		
			self.cnvx_partition = [Poly(points) for points in points_list]
		elif algo == 2:
			points_list = cgal_partition.hertel_mehlhorn(self.points)
			self.cnvx_partition = [Poly(points) for points in points_list]
		elif algo == 4:
			self.cnvx_partition = self.hertel_mehlhorn()
		elif algo == 5:
			self.cnvx_partition = self.triangulate()
			
		#also partition the holes
		for hole in self.holes:
			hole.partition(algo)
		
	def triangulate(self):
		#create the data structure required by the triangle library
		n = len(self.points)
		points = self.points
		edges = [[i, (i+1)%n] for i in range(n)]
		hole_pts = []
		for hole in self.holes:
			offset = len(points)
			hn = len(hole.points)
			edges += [[i+offset, ((i+1)%n) + offset] for i in range(hn)]
			points = np.concatenate((points, hole.points), axis=0)
			hole_pts.append(hole.get_interior_point())
			
		tri_dict = {'vertices': points, 'segments':edges}
		if len(self.holes) > 0:
			tri_dict['holes'] = hole_pts 

		#triangulate with the triangle library
		assert len(points) > 2, 'Trying to triangulate fewer than 2 points: ' + str(self)
		t = triangle.triangulate(tri_dict)
		t_verts = t['vertices'].tolist()
		tris = t['triangles'].tolist()
		partition = [Poly(np.array([t_verts[vi] for vi in tr])) for tr in tris]

		#remove triangles corresponding to holes/areas outside of the polygon
		filtered = []
		for poly in partition:
			pt = poly.get_interior_point()
			if self.contains_point(pt):
				filtered.append(poly)
		return filtered
			
	def hertel_mehlhorn(self):
		#first, triangulate
		partition = self.triangulate()
		
		#now for the actual HM
		
		#Setup data structures need for HM
		edges = []
		edge_2_poly={}
		fixed_edges = self.edges()
		for hole in self.holes:
			fixed_edges += hole.edges()
		
		for i in range(len(partition)):
			poly = partition[i]
			for edge in poly.edges():
				if edges not in fixed_edges:
					if edge in edges:
						edge_2_poly[edges.index(edge)].append(i)
					else:
						idx = len(edges)
						edges.append(edge)
						edge_2_poly[idx] = [i]
						
		#Now cycle through edges and remove what we can
		removing_edges = True
		mask = [i for i in range(len(edges))]
		while removing_edges:
			new_mask = []
			for m in mask:
				did_merge = False
				edge = edges[m]
				polys = [partition[i] for i in edge_2_poly[m]]
				if len(polys) ==2 : 
					poly0 = polys[0]
					poly1 = polys[1]
					e_0 = poly0.edges()
					e_1 = poly1.edges()
					idx0 = e_0.index(edge)
					idx1 = e_1.index(edge)
					plyc0  = PolyConstructor(np.roll(poly0.points, len(poly0.points)-idx0-1, axis=0))
					plyc1  = PolyConstructor(np.roll(poly1.points, len(poly1.points)-idx1-1, axis=0))
					
					assert plyc0.try_merge(plyc1), 'Merge should always be succesful as polys share an edge'
					new_poly = Poly(plyc0.points)
					if new_poly.is_convex():
						did_merge = True
						poly0.setPoints(new_poly.points)
						#Have edges that linked to poly1 link to updated poly0
						partition[edge_2_poly[m][1]] = poly0
					if not did_merge:
						new_mask.append(m)
						
			if len(new_mask) == len(mask):
				removing_edges = False
			else:
				mask = new_mask
		#remove duplicate objects
		partition = list(set(partition))
		
		return partition
		
	def setPoints(self, pts):
		assert len(pts) > 2, 'Tried to give polygon fewer than 3 vertices: ' + str(pts)
		self.points = pts
		
	def edges(self):
		eds = []
		n = len(self.points)
		for i in range(n):
			eds.append( Poly.Edge(self.points[i], self.points[(i+1)%n]) ) 
		return eds	


	def reduce(self, eps=0.1):
		#rdp might actually handle pruning colinear points for us
		# but it seems to give us a marginal speedup, so leaving in for now
		self.prune_colinear()
		#n_orig = len(self.points)
		reduced= rdp(self.points, epsilon = eps)
		#print('reduced from %d to %d vertices'%(n_orig, len(reduced)))
		if len(reduced) > 2:
			self.setPoints(reduced)
		else:
			print('Info: Reduced Poly only contains two points. Reduction ignored\n' + str(self))
		
	def prune_colinear(self):
		pruning = True

		i = 0
		n = len(self.points)
		while pruning and n >3:#once we're down to three points, no need to continue reducing.

			p0 = self.points[(i-1)%n]
			p1 = self.points[(i)%n]
			p2 = self.points[(i+1)%n]

			#check to see if these points are on the same line
			if colinear(p0, p1, p2):
				#remove this one dude
				self.setPoints(np.delete( self.points, i%n, axis=0))
			else:
				i += 1
				pruning = (i < n)
				
			n = len(self.points)
	
			
	def get_interior_point(self):
		p0 = self.points[0]
		p1 = self.points[1]
		p2 = self.points[2]
		a1 = p0 - p1
		a2 = p2 - p1
		theta = (np.arctan2(a1[1], a1[0]) + np.arctan2(a2[1], a2[0]))/2
		v = np.array([np.cos(theta), np.sin(theta)])
		r = 0.001
		
		pi = p1 + r*v
		
		if not self.contains_point(pi):
			pi = p1 - r*v
			if not self.contains_point(pi):
				print('WARNING: get_interior_point did not get an interior point')
		return pi


	def plot(self, show_partition = True):
		if show_partition and self.cnvx_partition is not None:
			for cnvx in self.cnvx_partition:
				cnvx.plot(show_partition)
		else:
			plot_points = self.points
			p0 = self.points[0]
			plot_points = np.concatenate((plot_points, [p0]), axis=0)
			plt.plot(plot_points[:,0], plot_points[:,1], 'k-')
		
	
	def is_cc(self):
		return self._signed_area() < 0

	def make_cc(self):
		if not self.is_cc():
			self.setPoints(np.flip(self.points, axis=0))

	def _signed_area(self):
		area = _signed_area(self.points)
		return area
		
	#Auxillary data structure for Hertel-Mehlhorn algorithm
	class Edge:
	
		def __init__(self, p0, p1):
			#enfore point order
			p0o = p0
			p1o = p1
			if p0[0]> p1[0] or (p0[0] == p1[0] and p0[1] > p1[1]):
				p0o = p1
				p1o = p0
				
			self.points=np.array([p0o, p1o])
			
		def __str__(self):
			return 'Poly.Edge: ' + str(self.points)
			
		def __eq__(self, obj):
			is_equal = False
			if isinstance(obj, Poly.Edge):
				dif = self.points - obj.points
				is_equal = np.sum(abs(dif)) == 0					
			return is_equal
			
		def as_constraint(self):
			p0 = self.points[0]
			d = self.points[1] - p0

			#x coords are identical, we're dealing with a vertical line
			if d[0] == 0:
				A = np.array([1,0])
				b = p0[0]
			elif d[1] == 0:
				A = np.array([0,1])
				b = p0[1]
			else:
				#ok this is just a normal guy
				A = np.array([-1*d[1]/d[0], 1])
				b = ( A @ p0)
			return A, b
				


#Exposing several geometry related functions which may be useful in general
#May make more sense to move these to a separate file


def colinear(p0, p1, p2):
	a1 = p1-p0
	a2 = p2-p0
	#because of numerical issues, we need to check to see if this is close
	eps = 0.01
	return abs(np.arctan2(a1[1], a1[0]) - np.arctan2(a2[1], a2[0])) < eps


def point_in_poly(verts, pt):
	#just put another point way out there	
	end_pt = np.max(verts, axis=0) + 10 
	n = len(verts)
	
	n_inter = 0
	
	for i in range(n):
		a1 = verts[i]
		b1 = verts[(i+1)%n]
		n_inter += segments_intersect(a1, b1, pt, end_pt)

	#odd number of intersections indicates point was in interior
	return (n_inter%2 == 1)

def segments_intersect(a1, b1, a2, b2):
	return is_ccw(a1, b1, a2) != is_ccw(a1, b1, b2) and is_ccw(a2, b2, a1) != is_ccw(a2, b2, b1)
		
def is_ccw(p0, p1, p2):
	return (p1[0] - p0[0])*(p2[1] - p1[1]) - (p2[0] - p1[0])*(p1[1] - p0[1]) > 0	


def poly_area(points):
    return abs(_signed_area(points))
    
def _signed_area(points):
	n = len(points)
	area = 0
	for i in range(n):
		pi = points[i]
		pj = points[(i+1)%n]
		area += (pj[0] - pi[0])*(pj[1] + pi[1])
	area = area/2
	
	return area		
		
def to_CGAL_points(points):
    CGAL_pts = []
    for pt in points:
        CGAL_pt = Point_2(pt[0], pt[1])
        CGAL_pts.append(CGAL_pt)
        
    return CGAL_pts
