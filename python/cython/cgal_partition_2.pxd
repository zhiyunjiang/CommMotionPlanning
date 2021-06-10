from libcpp.vector cimport vector
from libcpp.pair cimport pair

cdef extern from "cgal_partition_2_wrapper.h":
	ctypedef vector[pair[float, float]] PyList

	vector[PyList] optmal_convex_parition_from_list(PyList pt_pts)
	
	vector[PyList] greene_approx_convex_partition_from_list(PyList pt_pts)

	vector[PyList] hertel_mehlhorn_from_list(PyList pt_pts)
