#Calculate baseline policy by solving the traveling salesman policy.
import numpy as np

#solve using Google's OR tools
from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp

def calcGTSPPolicy(pt_lists):
	"""
	Following the discussion in "The Traveling Salesman Problem and Its Variations", (Gutin and Punnen, 2007), we reformaulte our GTSP as a standard TSP by setting the distance between points in the same set as large negative values.
	"""
	data = _buildProblemData(pt_lists)
	
	manager = pywrapcp.RoutingIndexManager(len(data['distance_matrix']),
					data['num_vehicles'], data['depot'])
	routing = pywrapcp.RoutingModel(manager)


	def _distanceCallback(from_index, to_index):
		"""Returns the distance between the two nodes."""
		# Convert from routing variable Index to distance matrix NodeIndex.
		from_node = manager.IndexToNode(from_index)
		to_node = manager.IndexToNode(to_index)
		dist = data['distance_matrix'][from_node][to_node]
		return dist
		
		
	transit_callback_index = routing.RegisterTransitCallback(_distanceCallback)
	
	routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)
	
	search_parameters = pywrapcp.DefaultRoutingSearchParameters()
	search_parameters.first_solution_strategy = (routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)
	
	solution = routing.SolveWithParameters(search_parameters)
	if solution:
		print_solution(manager, routing, solution)
		
		
	routes = get_routes(solution, routing, manager)
	# Display the routes.
	for i, route in enumerate(routes):
		print('Route', i, route)
	


    
def print_solution(manager, routing, solution):
	"""Prints solution on console."""
	print('Objective: {} miles'.format(solution.ObjectiveValue()))
	index = routing.Start(0)
	plan_output = 'Route for vehicle 0:\n'
	route_distance = 0
	while not routing.IsEnd(index):
		plan_output += ' {} ->'.format(manager.IndexToNode(index))
		previous_index = index
		index = solution.Value(routing.NextVar(index))
		route_distance += routing.GetArcCostForVehicle(previous_index, index, 0)
	plan_output += ' {}\n'.format(manager.IndexToNode(index))
	print(plan_output)
	plan_output += 'Route distance: {}miles\n'.format(route_distance)
    
    
def get_routes(solution, routing, manager):
	"""Get vehicle routes from a solution and store them in an array."""
	# Get vehicle routes and store them in a two dimensional array whose
	# i,j entry is the jth location visited by vehicle i along its route.
	routes = []
	for route_nbr in range(routing.vehicles()):
		index = routing.Start(route_nbr)
		route = [manager.IndexToNode(index)]
	while not routing.IsEnd(index):
		index = solution.Value(routing.NextVar(index))
		route.append(manager.IndexToNode(index))
		routes.append(route)
	return routes
    

def _buildProblemData(pt_lists):
	M = 50000
	partition_sizes = [len(pts) for pts in pt_lists] 
	bookmarks = np.cumsum(partition_sizes)
	n = np.sum(bookmarks[-1]) #size of distance matrix, 0 indexed
	#actually need to "shift" bookmarks back
	bookmarks = np.insert(bookmarks, 0, 0)

	D = np.ones((n,n))
	all_pts = np.empty((0,2))
	for pts in pt_lists:
		all_pts = np.concatenate((all_pts, pts))

	for i in range(n):
		id0 = np.where(bookmarks <= i)[0][-1]
		id1 = np.where(bookmarks > i)[0][0]
		ni = bookmarks[id1] - bookmarks[id0]
		ip1 = bookmarks[id0] + ( (i-bookmarks[id0]+1) % ni)
		for j in range(n):
			D[i,j] = np.linalg.norm(all_pts[ip1] - all_pts[j], axis=0)
		#and set the "next" node in set to -M
		D[i, ip1] = 0

	print(D)	
	data = {}
	data['distance_matrix'] = D
	data['num_vehicles'] = 1
	data['depot'] = 0	

	return data
	
	
	
#GUROBI imports
import gurobipy as gp
from gurobipy import GRB

def calcGTSPPolicyWMathProg(pt_lists):
	partition_sizes = [len(pts) for pts in pt_lists] 
	bookmarks = np.cumsum(partition_sizes)
	bookmarks = np.insert(bookmarks, 0, 0)
	n = np.sum(bookmarks[-1]) #size of distance matrix, 0 indexed

	all_pts = np.empty((0,2))
	for pts in pt_lists:
		all_pts = np.concatenate((all_pts, pts))


	#Build out the optimization problem in Gurobi
	m = gp.Model('GTSP')
	#Let's keep this nice and quiet
	m.Params.OutputFlag = 0
	D = np.zeros((n,n))
	wxs=[]
	xs=[]#binary variables, 1 - edge inlucded in optimal cycle, 0 - otherwise
	Mxs = np.empty((n,n))
	ys = [m.addVar(vtype=GRB.BINARY) for i in range(n)]#variables indicating whether each node is included in optimal cycle
	
	for i in range(n):
		for j in range(i,n):
			D[i,j] = np.linalg.norm(all_pts[i] - all_pts[j], axis=0)
			D[j,i] = D[i,j]
			#I could not include edge variables for points in same partition, but
			#leaving in for now to simplify coding
			xs.append(m.addVar(vtype=GRB.BINARY))
			wxs.append(D[i,j]*xs[-1])
			Mxs[i,j] = xs[-1]
			Mxs[j,i] = xs[-1]
			#constraints related to edges/nodes
			
	for i in range(n):
		#add constraint
		one_in = np.concatenate((Mxs[:i],Mxs[i+1:]))
		m.addLConstr(np.sum(one_in), GRB.EQUAL, 2*ys[i])
		
	for i in range(len(pt_lists)):
		partition_ys = ys[bookmarks[i]:bookmarks[i+1]]
		m.addLConstr(np.sum(partition_ys)), GRB.EQUAL, 1)
		
	#TODO add last constraint
	
	
	m.setObjective(np.sum(wxs))

	#and then solve
	

	











	
