#include <utility>
#include <vector>

typedef std::vector<std::pair<float, float>>	PyList;

//Just want access to these three CGAL partitioning functions
std::vector<PyList> optmal_convex_parition_from_list(PyList pt_pts);

std::vector<PyList> greene_approx_convex_partition_from_list(PyList pt_pts);

std::vector<PyList> hertel_mehlhorn_from_list(PyList pt_pts);

