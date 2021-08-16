import pyomo.environ as pymo
from pyomo.opt import SolverFactory, TerminationCondition
import numpy as np

def min_PWD(regions, pi, formulation=0, verbose = False):
	if formulation == 0:
		return _min_PWD_MICP(regions, pi, verbose)
	else:
		return _min_PWD_MIQCP(regions, pi, verbose)

def min_cycle(regions, order, verbose = False):
	model = pymo.ConcreteModel()

	model.n = len(regions)

	model.order = order
	
	_add_x_vars(model, regions)

	#Setup s dummy variables
	model.s = pymo.Var([i for i in range(model.n)], bounds=(0, None))

	def AvgSwtichingTime(model):
		expr=0
		for i in range(model.n):
			expr+=model.s[i]
		return expr

	model.o = pymo.Objective(rule=AvgSwtichingTime)

	#Constraint x and y to s
	def C_StoXY(model):
		l=0
		order = model.order
		x = model.x
		y = model.y
		n = model.n
		s = model.s
		for i in range(model.n):
			yield (x[order[i]]-x[order[(i+1)%n]])**2 + (y[order[i]]-y[order[(i+1)%n]])**2 == s[i]**2

	model.Dummy = pymo.ConstraintList(rule=C_StoXY)

	_add_region_constraints(model, regions)

	if verbose:
		model.pprint()

	return _solve(model, verbose)

def min_PWD_over_P(pi, D, verbose = False):
	model = pymo.ConcreteModel()

	n = len(pi)
	model.n = n
	model.pi = pi
	model.D = D
	model.P = pymo.Var([ i*n +j for i in range(n) for j in range(n) if i!=j], bounds=(0,1))

	model.irreducible_aperiodic = pymo.ConstraintList()
	
	for i in range(n):
		for j in range(n):
			if i!=j:
				model.irreducible_aperiodic.add(0.0001<=model.P[i*n +j])
				
	def detail_balance(model):
		pi = model.pi
		for i in range(model.n):
			for j in range(model.n):
				if i!=j:
					yield pi[i]*model.P[i*n +j] - pi[j]*model.P[j*n +i] == 0
	model.detail_balance = pymo.ConstraintList(rule=detail_balance)

	def valid_dist(model):
		for i in range(model.n):
			expr = 0
			for j in range(model.n):
				if i!=j:
					expr+=model.P[model.n*i + j]
			yield expr ==1
	model.valid_dist = pymo.ConstraintList(rule=valid_dist)

	def AvgSwtichingTime(model):
		return sum([model.pi[i]*model.P[i*n +j]*model.D[i,j] for i in range(model.n) for j in range(model.n) if i!=j])
	model.o = pymo.Objective(rule=AvgSwtichingTime)

	if verbose:
		model.pprint()

	with SolverFactory('SHOT') as opt:
		results = opt.solve(model, load_solutions = False, tee=verbose)
		if results.solver.termination_condition != TerminationCondition.optimal:
			print('Optimal solution was not found')
			argmin = None
			val = float('inf')
		else:
			model.solutions.load_from(results)
			val = model.o()
			argmin = np.zeros((n,n))
			for i in range(model.n):
				for j in range(model.n):
					if i!=j:
						argmin[i,j]=model.P[i*n +j]()

	return argmin, val




#Private functions

def _min_PWD_MICP(regions, pi, verbose=False):
	model = pymo.ConcreteModel()

	model.n = len(regions)

	model.pi = pi
	
	_add_x_vars(model, regions)

	def AvgSwtichingTime(model):
		terms = [ model.pi[i]*model.pi[j]*(((model.x[i]-model.x[j])**2 + (model.y[i]-model.y[j])**2)**0.5) for i in range(model.n) for j in range(i+1, model.n)]
		return sum(terms)

	model.o = pymo.Objective(rule=AvgSwtichingTime)


	_add_region_constraints(model, regions)

	if verbose:
		model.pprint()

	return _solve(model, verbose)

def _min_PWD_MIQCP(regions, pi, verbose=False):
	model = pymo.ConcreteModel()

	model.n = len(regions)

	model.pi = pi
	
	_add_x_vars(model, regions)

	#Setup s dummy variables
	model.s = pymo.Var([i for i in range(int(model.n*(model.n+1)/2))], bounds=(0, None))

	def AvgSwtichingTime(model):
		l = 0
		expr=0
		for i in range(model.n):
			for j in range(i+1, model.n):
				expr+=model.pi[i]*model.pi[j]*model.s[l]
				l+=1
		return expr

	model.o = pymo.Objective(rule=AvgSwtichingTime)

	#constrain x and y to s
	def C_StoXY(model):
		l=0
		for i in range(model.n):
			for j in range(i+1, model.n):
				yield (model.x[i]-model.x[j])**2 + (model.y[i]-model.y[j])**2 == model.s[l]**2
				l+=1

	model.Dummy = pymo.ConstraintList(rule=C_StoXY)

	_add_region_constraints(model, regions)

	if verbose:
		model.pprint()

	return _solve(model, verbose)


def _solve(model, verbose = False):
	with SolverFactory('SHOT') as opt:
		results = opt.solve(model, load_solutions = False, tee=verbose)
		if results.solver.termination_condition != TerminationCondition.optimal:
			print('Optimal solution was not found')
			argmin = None
			val = float('inf')
		else:
			model.solutions.load_from(results)
			val = model.o()
			argmin = []
			for i in range(model.n):
				argmin.append([model.x[i](), model.y[i]()] )

	return np.array(argmin), val


def _add_x_vars(model, regions):
	#TODO - possibly add bounds
	#add x and y variables
	def xb(model, i, regions):
		return (regions[i].xmin, regions[i].xmax)
	model.x = pymo.Var([i for i in range(model.n)], bounds=lambda model, i: xb(model, i, regions))

	def yb(model, i, regions):
		return (regions[i].ymin, regions[i].ymax)
	model.y = pymo.Var([i for i in range(model.n)], bounds=lambda model, i: yb(model, i, regions))

def _add_region_constraints(model, regions):
	#constrain (x,y) to be in relay regions
	mks = [reg.mk for reg in regions]
	eta_idxs = np.cumsum(mks)
	model.eta = pymo.Var([i for i in range(eta_idxs[-1])], domain=pymo.Binary)


	#Constrain xi to be in at least one of the subgregions of relay region i
	model.InRegions=pymo.ConstraintList()
	#at least one must be set
	start = 0
	end = 0
	for i in range(model.n):
		start = end
		end = eta_idxs[i]
		expr=0
		for j in range(start, end):
			expr += model.eta[j]
		model.InRegions.add(expr >= 1)

	#and finally, our hefty linear constraints

	model.CnvxPoly = pymo.ConstraintList()
	C = 100000 #actually need to find a value of this constant
	l = 0
	for i in range(model.n):
		reg =regions[i]
		for poly in reg.polygons:
			for cnvx in poly.cnvx_partition:
				Aik, bik = cnvx.to_linear_constraints()
				bik[:,0] += C

				for j in range(len(bik)):
					model.CnvxPoly.add(Aik[j,0]*model.x[i] + Aik[j,1]*model.y[i] + C*model.eta[l] <= bik[j,0])
				l+=1