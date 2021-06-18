import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from cgal_partition_2 cimport PyList, optmal_convex_parition_from_list, greene_approx_convex_partition_from_list, hertel_mehlhorn_from_list


def optmal_convex_parition(vertices):
	cdef PyList c_verts = py_to_cpp_verts(vertices)
	#call the magic function
	c_part_verts = optmal_convex_parition_from_list(c_verts)
	return cpp_partition_to_py_verts(c_part_verts)
	
def greene_approx_convex_partition(vertices):
	cdef PyList c_verts = py_to_cpp_verts(vertices)
	#call the magic function
	c_part_verts = greene_approx_convex_partition_from_list(c_verts)
	return cpp_partition_to_py_verts(c_part_verts)
	
def hertel_mehlhorn(vertices):
	cdef PyList c_verts = py_to_cpp_verts(vertices)
	#call the magic function
	c_part_verts = hertel_mehlhorn_from_list(c_verts)
	return cpp_partition_to_py_verts(c_part_verts)

cdef PyList py_to_cpp_verts(verts):
	cdef PyList cverts;
	for vert in verts:
		cverts.push_back(vert)
		
	return cverts
	
cdef cpp_partition_to_py_verts(vector[PyList] cverts): 
	partition_polys = []
	for pverts in cverts:
		n_verts = pverts.size()
		np_verts = np.zeros((n_verts,2))
		for i in range(n_verts):
			vert = pverts[i]
			np_verts[i, 0] = vert.first
			np_verts[i, 1] = vert.second
			
		partition_polys.append(np_verts)
	return partition_polys


