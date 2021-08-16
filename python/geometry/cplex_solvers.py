import numpy as np


#import CPLEX libraries
from docplex.mp.model import Model
from docplex.mp.constr import QuadraticConstraint, ComparisonType

def min_PWD(regions, pi):
		""""
		IBM's CPlex studio cannot handle quadratic equality constraints nor anything beyond quadratic terms in the objective,
		so this method will always result in an error. See 
		https://www.ibm.com/docs/en/icos/12.7.1.0?topic=smippqt-miqcp-mixed-integer-programs-quadratic-terms-in-constraints for details.
		"""
	
			
		mdl = Model(name='min_PWD_cplx')
		n_regions = len(regions)
		s=[]
		k=[]
		for i in range(n_regions):
			for j in range(i+1, n_regions):
				s.append(mdl.continuous_var())
				k.append(pi[i]*pi[j])

		mdl.minimize(mdl.sum([k[i]*s[i] for i in range(len(s))]))

		#set up variables needed for objective
		xs = mdl.continuous_var_list(n_regions, lb=[reg.xmin for reg in regions], ub=[reg.xmax for reg in regions])
		ys = mdl.continuous_var_list(n_regions, lb=[reg.ymin for reg in regions], ub=[reg.ymax for reg in regions])

		#add constraints to couple new dummy vars to relay postitions
		qcs=[]
		for i in range(n_regions):
			for j in range(i+1, n_regions):
				qcs.append( QuadraticConstraint(mdl, mdl.sum_squares( [(xs[i]-xs[j]), (ys[i]-ys[j])] ), ComparisonType.EQ, mdl.sum_squares([s[len(qcs)]]) ))
		#Not sure this is necessary        
		mdl.add_quadratic_constraints(qcs)


		C = 10000 #actually need to find a value of this constant
		#Constrain relay points to lie within their regions
		for i in range(n_regions):
			reg = regions[i]    
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


		sol = mdl.solve()
		argmin = None
		min_val = float('inf')
		if sol is not None:
			#extract x and y values
			argmin = []
			for i in range(n_regions):
				argmin.append([xs[i].solution_value(), ys[i].solution_value()])
			argmin = np.array(argmin)
			#extract min val
			min_val = sol.get_objective_value()

		return argmin, min_val