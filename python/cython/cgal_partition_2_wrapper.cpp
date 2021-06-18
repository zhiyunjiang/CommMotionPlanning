#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <cassert>
#include <list>
#include <utility>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Partition_traits_2<K>                         Traits;
typedef Traits::Polygon_2                                   Polygon_2;
typedef Traits::Point_2                                     Point_2;
typedef std::list<Polygon_2>                                Polygon_list;
typedef std::vector<std::pair<float, float>>		      PyList;

void polygon_from_list(Polygon_2& polygon, PyList& py_pts)
{
	for (size_t i=0; i<py_pts.size(); i++){
		std::pair<float, float>  pt = py_pts[i];
		polygon.push_back(Point_2(pt.first, pt.second));
	}
}

void list_from_polygons(Polygon_list& partition_list, std::vector<PyList>& vertices)
{
	for (Polygon_2 poly : partition_list){

		//iterate over the polygon's vertices
		PyList pVerts;
		for (Point_2 p : poly.container()) {
			//format as pair
			std::pair<float, float> py_pt (p.x(), p.y());
			//push onto vector
			pVerts.push_back(py_pt);
		}
		vertices.push_back(pVerts);
	}
}

std::vector<PyList> optmal_convex_parition_from_list(PyList py_pts) {
	//make the polygon from the list
	Polygon_2             polygon;
	Polygon_list          partition_polys;
	polygon_from_list(polygon, py_pts);
	
	CGAL::optimal_convex_partition_2(polygon.vertices_begin(),
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));
	
	assert(CGAL::partition_is_valid_2(polygon.vertices_begin(),
                                     polygon.vertices_end(),
                                     partition_polys.begin(),
                                     partition_polys.end()));
                                     
	//now convert partition_polys into a list of lists
	//each list representing a polygon
	std::vector<PyList> vertices;
	list_from_polygons(partition_polys, vertices);
	
	return vertices;
}

std::vector<PyList> greene_approx_convex_partition_from_list(PyList py_pts){

	Polygon_2             polygon;
	Polygon_list          partition_polys;
	polygon_from_list(polygon, py_pts);
	
	CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(),
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));
	
	assert(CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(),
                                     polygon.vertices_end(),
                                     partition_polys.begin(),
                                     partition_polys.end()));
                                     
	//now convert partition_polys into a list of lists
	//each list representing a polygon
	std::vector<PyList> vertices;
	list_from_polygons(partition_polys, vertices);
	
	return vertices;
}

//approx_convex_partition_2
std::vector<PyList> hertel_mehlhorn_from_list(PyList py_pts){
	Polygon_2             polygon;
	Polygon_list          partition_polys;
	polygon_from_list(polygon, py_pts);
	
	CGAL::approx_convex_partition_2(polygon.vertices_begin(),
                                    polygon.vertices_end(),
                                    std::back_inserter(partition_polys));
	
	assert(CGAL::partition_is_valid_2(polygon.vertices_begin(),
                                     polygon.vertices_end(),
                                     partition_polys.begin(),
                                     partition_polys.end()));
                                     
	//now convert partition_polys into a list of lists
	//each list representing a polygon
	std::vector<PyList> vertices;
	list_from_polygons(partition_polys, vertices);
	
	return vertices;
}
