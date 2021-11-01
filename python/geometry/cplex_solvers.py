import numpy as np


#import CPLEX libraries
from docplex.mp.model import Model
from docplex.mp.constr import QuadraticConstraint, ComparisonType

import sys  
from pathlib import Path
path = Path(__file__)
top_dir = path.parent.parent
sys.path.insert(0, str(top_dir.absolute())+"/utils")
import fssmc as fmc

def min_PWD(regions, W, policy_type = 0):
	n_regions = len(regions)
	
	if policy_type == 0:
		pi = W
		P = np.tile(pi,(n_regions,1))
	elif policy_type == 1:
		P = W
		pi = fmc.stationary(P)
		
	
	mdl = Model(name='min_PWD_cplx')
	s=[]
	k=[]
	for i in range(n_regions):
		for j in range(i+1, n_regions):
			s.append(mdl.continuous_var())
			k.append(pi[i]*P[i,j])

	mdl.minimize(mdl.sum([k[i]*s[i] for i in range(len(s))]))

	xs, ys = _add_xs_and_ys(n_regions, regions, mdl)

	#add constraints to couple new dummy vars to relay postitions
	qcs=[]
	for i in range(n_regions):
		for j in range(i+1, n_regions):
			qcs.append( QuadraticConstraint(mdl, mdl.sum_squares( [(xs[i]-xs[j]), (ys[i]-ys[j])] ),
											ComparisonType.LE, mdl.sum_squares([s[len(qcs)]]) ))      
	mdl.add_quadratic_constraints(qcs)

	_add_regional_constraints(n_regions, regions, mdl, xs, ys)

	return _solve(mdl, xs, ys, n_regions)

def min_cycle(regions, order, verbose = False):
	mdl = Model(name='min_PWD_cplx')
	n_regions = len(regions)

	s = mdl.continuous_var_list(n_regions)

	mdl.minimize(mdl.sum([s[i] for i in range(len(s))]))

	xs, ys = _add_xs_and_ys(n_regions, regions, mdl)

	#add constraints to couple new dummy vars to relay postitions
	qcs = []
	for i in range(n_regions):
		j = (i+1)%n_regions
		qcs.append( QuadraticConstraint(mdl, mdl.sum_squares( [(xs[order[i]]-xs[order[j]]), (ys[order[i]]-ys[order[j]])] ),
										ComparisonType.LE, mdl.sum_squares([s[i]]) ))       
	mdl.add_quadratic_constraints(qcs)

	_add_regional_constraints(n_regions, regions, mdl, xs, ys)

	return _solve(mdl, xs, ys, n_regions)


def _solve(mdl, xs, ys, n_regions):
	sol = mdl.solve()
	argmin = None
	min_val = float('inf')
	if sol is not None:
		#extract x and y values
		argmin = []
		for i in range(n_regions):
			argmin.append([xs[i].solution_value, ys[i].solution_value])
		argmin = np.array(argmin)
		#extract min val
		min_val = sol.get_objective_value()
	return argmin, min_val

def _add_xs_and_ys(n_regions, regions, mdl):
	#set up variables needed for objective
	xs = mdl.continuous_var_list(n_regions, lb=[reg.xmin for reg in regions], ub=[reg.xmax for reg in regions])
	ys = mdl.continuous_var_list(n_regions, lb=[reg.ymin for reg in regions], ub=[reg.ymax for reg in regions])
	return xs, ys


def _add_regional_constraints(n_regions, regions, mdl, xs, ys):
	C = 10000 #actually need to find a value of this constant, for now just make really big, adjust if needed
	#Constrain relay points to lie within their regions
	for i in range(n_regions):

		reg = regions[i]    
		#constrain to be within the bounding box, which should speed things up
		mdl.add_constraint_(mdl.ge_constraint(xs[i], reg.xmin))
		mdl.add_constraint_(mdl.le_constraint(xs[i], reg.xmax))
		mdl.add_constraint_(mdl.ge_constraint(ys[i], reg.ymin))
		mdl.add_constraint_(mdl.le_constraint(ys[i], reg.ymax))
		mk=0
		eta_i = []
		for poly in reg.polygons:
			for cnvx in poly.cnvx_partition:
				eta_ik = mdl.binary_var()
				eta_i.append(eta_ik)
				mk+=1
				Aik,bik = cnvx.to_linear_constraints()
			
				bik[:,0] += C

				for j in range(len(bik)):
					#add the linear constraint
					mdl.add_constraint_(mdl.le_constraint(Aik[j,0]*xs[i] + Aik[j,1]*ys[i] +C*eta_ik ,bik[j,0]))
	
		#we must be in exactly one of the subregions
		mdl.add_constraint(mdl.eq_constraint(mdl.sum(eta_i), 1))