import numpy as np

#The MIQCP Solver
import gurobipy as gp
from gurobipy import GRB

def min_cycle(regions, P):
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
	xs = [m.addVar(lb=reg.xmin, ub=reg.xmax) for reg in regions]
	ys = [m.addVar(lb=reg.ymin, ub=reg.ymax) for reg in regions]

	#add constraints to couple new dummy vars to relay postitions
	for i in range(n):
		edge_s = P[i]
		edge_e = P[(i+1)%n]
		#add the quadratic constraint that ||xi - xj||^2 = sij^2
		M = np.array([[1,0,0, 0, 0],[0,-1,1, 0, 0],[0,1,-1, 0, 0], [0, 0, 0, -1, 1], [0, 0, 0, 1, -1]])
		xc=[s[i], xs[edge_s], xs[edge_e], ys[edge_s], ys[edge_e]]
		m.addMQConstr(M,None, GRB.EQUAL, 0, xc, xc)

	_add_regional_constraints(m, regions, n, xs, ys)

	return _extract_solution(m, n, xs, ys)

def min_PWD(regions, pi, formulation=0, verbose = False):
	if formulation == 0:
		return _min_PWD_MIQCP(regions, pi, verbose)
	else:
		return _min_PWD_PWLC(regions, pi, verbose)
"""
Formulate the problem using piecewise linear constraints (PWLC)
"""
def _min_PWD_PWLC(regions, pi, verbose = False):
	#use a series of piecewise linear constraints rather than binary variables
	
	#To use Gurboi, the objective must be linear or quadratic	
	m = gp.Model('min_PWD_PWLC_gb')
	#Let's keep this nice and quiet
	if not verbose:
		m.Params.OutputFlag = 0
	#setup dummy variables required for linear objective
	n_regions = len(regions)
	
	xs, ys = _add_x_vars(m, regions, n_regions, pi)

	#constrain relay points to lie within their regions
	#using a series of linear and piecewise linear equality constraints
	for i in range(n_regions):
		reg = regions[i]
		theta_tildes = []
		k=0
		for poly in reg.polygons:
			for cnvx in poly.cnvx_partition:
				k+=1
				Aik,bik = cnvx.to_linear_constraints()
				#recall Aik [x,y]^T -b<0 if [x,y] in polygon
				#=> -Aik[x,y]^T +b > 0 if [x,y] in polygon
				#=>z = -Aik[x,y]^T +b => z+Aik[x,y]=b
				n_edges = len(bik)
				Aik_aug = np.zeros((n_edges,2+n_edges))
				Aik_aug[:,:2] = Aik
				Aik_aug[:,2:] = np.eye(n_edges)

				zs = [m.addVar(lb=float('-inf')) for j in range(n_edges)]
				m.addMConstr(Aik_aug, [xs[i], ys[i], *zs ], GRB.EQUAL, bik)
				
				ztildes = [m.addVar(lb=float('-inf'), name='z tilde %d.%d.%d'%(i,k,j)) for j in range(n_edges)]
				for j in range(n_edges):
					#first of the piecewise constraints
					m.addGenConstrPWL(zs[j], ztildes[j], [-1,0,1],[-1000, 1, 1])
				
				#dummy variable for checking how many constraints are satisfied
				theta = m.addVar(lb=float('-inf'), name='theta %d.%d'%(i,k))
				a = np.ones((1, n_edges+1))
				a[0,-1] = -1
				m.addMConstr(a, [*ztildes, theta], GRB.EQUAL, np.array([0]))
				
				# saturate so that theta_tilde is 1 if all constraints satisfied,
				# ~0 otherwise
				theta_tilde = m.addVar(name='theta tilde %d.%d'%(i,k))
				theta_tildes.append(theta_tilde)
				m.addGenConstrPWL(theta, theta_tilde, [0,n_edges-0.001,n_edges,n_edges+1],
							[0,0, 1, 1])
					
		#sum of theta tildes should be greater than 1,
		# i.e., all constraints for at least one cnvx region should be satisfied
		a = np.ones((1,len(theta_tildes)))
		m.addMConstr(a, theta_tildes, GRB.GREATER_EQUAL, np.array([1]))
				
	    
	return _extract_solution(m, n_regions, xs, ys)

"""
Formulate the prolbem as MIQCP
"""
def _min_PWD_MIQCP(regions, pi, verbose = False):
	#To use Gurboi, the objective must be linear or quadratic	
	m = gp.Model('min_PWD_MIQCP_gb')
	#Let's keep this nice and quiet
	if not verbose:
		m.Params.OutputFlag = 0
	#setup dummy variables required for linear objective
	n_regions = len(regions)
	
	xs, ys = _add_x_vars(m, regions, n_regions, pi)

	_add_regional_constraints(model, regions, n_regions, xs, ys)
		    
	return _extract_solution(m, n_regions, xs, ys)

#Shared setup

def _extract_solution(model, n_regions, xs, ys):
	#and now solve
	model.params.NonConvex = 2
	model.optimize()
	
	#TODO - feasability checking/handling
	#assert m.status == 2, '
	
	#and extract the optimal values
	x = []
	for i in range(n_regions):
		x.append([xs[i].x, ys[i].x])

	#return the points as well as the optimal value
	return np.array(x), model.getObjective().getValue()

def _add_regional_constraints(model, regions, n_regions, xs, ys):
	#constrain relay points to lie within their regions
	eta = []
	for i in range(n_regions):
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
	    
		eta_i = [model.addVar(vtype=GRB.BINARY) for k in range(n_vars-2)]
		eta.append(eta_i)
		model.addConstr(np.sum(eta_i) == 1)
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
	    
		LMC = model.addMConstr(A, [xs[i], ys[i] , *eta_i], GRB.LESS_EQUAL, b)

def _add_x_vars(model, regions, n_regions, pi):
	s = []
	s_obj = []
	for i in range(n_regions):
		for j in range(i+1, n_regions):
			s.append(model.addVar())
			s_obj.append(pi[i]*pi[j]*s[-1])	

	model.setObjective(np.sum(s_obj))
	xs = [model.addVar(lb=reg.xmin, ub=reg.xmax) for reg in regions]
	ys = [model.addVar(lb=reg.ymin, ub=reg.ymax) for reg in regions]
	
	#add constraints to couple new dummy vars to relay postitions
	ij = 0
	for i in range(n_regions):
		for j in range(i+1, n_regions): 
			#add the quadratic constraint that ||xi - xj||^2 = sij^2
			M = np.array([[1,0,0, 0, 0],[0,-1,1, 0, 0],[0,1,-1, 0, 0], [0, 0, 0, -1, 1], [0, 0, 0, 1, -1]])
			xc=[s[ij], xs[i], xs[j], ys[i], ys[j]]
			model.addMQConstr(M,None, GRB.EQUAL, 0, xc, xc)
			ij += 1

	return xs, ys