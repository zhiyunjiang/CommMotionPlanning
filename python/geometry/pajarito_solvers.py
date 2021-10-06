import julia 
import os

#pajarito isn't compatible with the latest version of JuMP, which seems to be causing problems. Ignoring this fow now

J_dir = os.path.dirname(os.path.abspath(__file__))

def min_PWD(regions, pi):

	j = julia.Julia(debug = True, compiled_modules = False)
	j_min_pwd = j.include(J_dir + '/min_pwd.jl')
	
	data = []

	n = len(pi)
	w = []
	for i in range(n):
		for j in range(i+1, n):
			w.append(pi[i]*pi[j])

	bounds = []
	poly_bounds = []
	m = 0
	for reg in regions:
		bounds.append([reg.xmin, reg.ymin,
		 				reg.xmax, reg.ymax])
		mk=0
		poly_bounds_i = []
		for poly in reg.polygons:
			for cnvx in poly.cnvx_partition:
				mk+=1
				Aik,bik = cnvx.to_linear_constraints()
			
				#augment Aik with row of Cs
				Aik_tilde = np.concatenate( (Aik,C*np.ones((len(bik),1))) )
				#augment bik with C
				bik[:,0] += C
				bounds_dict = {'A': Aik_tilde, 'b': b_ik}
				poly_bounds_i.append(bounds_dict)

		m += mk
		poly_bounds.append(poly_bounds_i)

	argmin, min_val = j_min_pwd(n, m, w, bounds, poly_bounds)

	return argmin, min_val