from itertools import permutations
import numpy as np

#The MIQCP Solver
import gurobipy as gp
from gurobipy import GRB


def _unique_cycles(n):
	permutes = permutations(range(n))
	pruned_permutes = []
	for permute in permutes:
		if permute[0] == 0:
			pruned_permutes.append(permute)
	return pruned_permutes


def _solve_one_instance(P, regions):
	n = len(P)
	#To use Gurboi, the objective must be linear or quadratic	
	m = gp.Model('min_pairwise_distance')
	#Let's keep this nice and quiet
	m.Params.OutputFlag = 0
	#setup dummy variables required for linear objective
	s = []
	for i in range(n):
		s.append(m.addVar())

	m.setObjective(np.sum(s))
	xs = [m.addVar() for i in range(n)]
	ys = [m.addVar() for i in range(n)]

	#add constraints to couple new dummy vars to relay postitions
	for i in range(n):
		edge_s = P[i]
		edge_e = P[(i+1)%n]
		#add the quadratic constraint that ||xi - xj||^2 = sij^2
		M = np.array([[1,0,0, 0, 0],[0,-1,1, 0, 0],[0,1,-1, 0, 0], [0, 0, 0, -1, 1], [0, 0, 0, 1, -1]])
		xc=[s[i], xs[edge_s], xs[edge_e], ys[edge_s], ys[edge_e]]
		m.addMQConstr(M,None, GRB.EQUAL, 0, xc, xc)

	#constrain relay points to lie within their regions
	eta = []
	for i in range(n):
		reg = regions[i]
		As = np.zeros((0,2))
		bs = []
		Cs = []

		n_vars = 2
		for poly in reg.polygons:
			for cnvx in poly.cnvx_partition:
				n_vars += 1
				Aik,bik = cnvx.to_linear_constraints()
				As = np.concatenate((As,Aik), axis = 0)
				C = 1*10000*np.ones(len(bik)) #actually need to find a value of this constant
				Cs.append(C)
				bs += (bik[:,0]+C).tolist()

		eta_i = [m.addVar(vtype=GRB.BINARY) for k in range(n_vars-2)]
		eta.append(eta_i)
		m.addConstr(np.sum(eta_i) == 1)
		eta.append(eta_i)
		n_constraints = len(bs)
		A = np.zeros((n_constraints, n_vars))
		A[:,:2] = As
		idx = 0
		for j in range(n_vars - 2):
			C = Cs[j]
			lc = len(C)
			A[idx:idx+lc, 2+j] = C
			idx += lc

		b = np.array(bs)

		LMC = m.addMConstr(A, [xs[i], ys[i] , *eta_i], GRB.LESS_EQUAL, b)

	#and now solve
	m.params.NonConvex = 2
	m.optimize()

	#TODO - feasability checking/handling
	#assert m.status == 2, '

	#and extract the optimal values
	x = []
	for i in range(n):
		x.append([xs[i].x, ys[i].x])
		
	return m.getObjective().getValue(), np.array(x)
	

def TSPN_BF(regions):
	n = len(regions)

	#for each unique cycle (exclude congruent, directionally symmetric cycles)
	#solve the average weighted distance problem using Gurobi

	bsf = np.inf
	argmin = None
	Ps = _unique_cycles(n)
	for P in Ps:
		print('Working on Permutation ' + str(P))
		val, x = _solve_one_instance(P, regions)
		#check if this is the best we've seen
		if val <= bsf:
			print('Optimal Solution Improved to %.2f'%(val))
			bsf = val
			argmin = (x, P)
	return bsf, argmin
