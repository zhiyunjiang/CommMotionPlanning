import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

#The MIQCP Solver
import gurobipy as gp
from gurobipy import GRB

#Import a few utility functions...
import sys 
from pathlib import Path
path = Path(__file__)
top_dir = path.parent.parent
sys.path.insert(0, str(top_dir.absolute())+"\\utils")
sys.path.insert(0, str(top_dir.absolute())+"\\geometry")
sys.path.insert(0, str(top_dir.absolute())+"\\polling_systems")
sys.path.insert(0, str(top_dir.absolute())+"\\comm_channel")
#So we can import my local libs
import CommChannel as CC
import PollingSystem as PS
import MarkovianRP as MRP
from CoordTransforms import toRawFromGrid
from CoordTransforms import toGridFromRaw
import pointcloud as PC

#Delay-Tolerant Relay System
class DTR:

	#TODO eventually, Beta should be calculated from comm specifications
	def __init__(self, channels, ls, beta, th = -92, p_th=0.7):
		assert len(channels)%2 == 0, "Expected even number of channels"
		assert len(channels) == 2*len(ls), "Mismatch: number of regions and number of queue parameters"

		self.n = len(ls)
		self.channels = channels
		self.p_th = p_th
		self.gamma_th = th
		self.region = channels[0].region
		self.res = channels[0].res
		self.findJointConnectivityField()
		self.Xis = [_indices_to_pts(cf, self.region, self.res) for cf in self.cfields]
		self.cregions = [ PC.PointCloud(Xi['points']) for Xi in self.Xis]
		for pc in self.cregions:
    			pc.partition(algo=4, alpha = 0.5/self.res)

		self.ps = PS.PollingSystem(ls, beta)

		
	def findJointConnectivityField(self):
		#old way to do it, which just needs both to be above the threshold
		#self.cfs = [cc.getConnectionField(th) for cc in channels]
		#self.cfields = [(i+1)*(self.cfs[2*i]*self.cfs[2*i+1]) for i in range(len(ls))]
		
		pfs = [cc.getPConField(self.gamma_th) for cc in self.channels]
		self.cfields = [1*((pfs[2*i]*pfs[2*i+1]) > self.p_th) for i in range(self.n)]

	def optimize(self, do_plot=True, x_opt_method = 0, verbose = False):
		converged = False
		#initialize policy and polling locations
		pi = MRP.build_ed_policy(self.ps.Ls, self.ps.beta)
		#randomly sample points from Xis
		X = self._pick_random_X()
		Xs = [X]
		S = XtoS(np.copy(X))
		it_count = 0
		max_it = 500
		Wprev = np.inf
		while (not converged) and (it_count < max_it):
			#first, find optiaml pi
			res = self.ps.calc_optiaml_rp(S)
			pi = res.x
			W = res.fun
			if it_count%50 == 0 and verbose:
				print("Optimized Policy Waiting Time: %.4f"%(W))
			#Now update the coordinates
			if x_opt_method == 0:
				X = self._simple_gradient(X, S, pi, _lr(it_count))
			elif x_opt_method == 1:
				X = self._PSO(X, S, pi)
			elif x_opt_method == 2:
				base_temp = 100
				X = self._Metropolis(X,S,pi, temp=100/(1 + it_count//25))
			elif x_opt_method ==3:
				X, _ = self._find_min_PWDalt(pi)
			Xs.append(np.copy(X))
			rp = MRP.RandomRP(pi)
			S =  XtoS(X)
			W = self.ps.calc_avg_wait(rp, S)
			if verbose:
				print('Transition probabilities: ', pi)
				print('Points: ', X)
				print("Optimized Location Waiting Time: %.4f"%(W))
			if it_count%50 == 0 and verbose:
				print("Optimized Location Waiting Time: %.4f"%(W))
			it_count += 1

			eps = 0.01
			if abs(W - Wprev) <= eps:
				converged = True
			Wprev = W
		if do_plot:
			self.plot_optimization(X)
		if verbose:
			print("Optimized Waiting Time: %.4f"%(Wprev))
		return W, pi, X
		
	def naive_policy(self):
		weights = np.ones(self.n)
		X,_ = self._find_min_PWD(weights)
		pi = self.ps.Ls/np.sum(self.ps.Ls)
		rp = MRP.RandomRP(pi)
		W = self.ps.calc_avg_wait(rp, XtoS(X))
		return W, pi, X

	def plot_optimization(self, x, save = False, ls=12):
		fig = plt.figure(figsize=(12,12))
		#plot the connectivity fields
		color_list=['r', 'g', 'b', 'c']
		for i in range(self.n):
			Xi = self.Xis[i]
			pts = Xi['points']
			plt.plot(pts[:,0],  pts[:,1], '.' + color_list[i], label='Relay Region %d'%(i+1))
			plt.plot(x[i,0], x[i,1], '+k', markersize=20)
			
		plt.xlim(self.region[1], self.region[0])
		plt.ylim(self.region[3], self.region[2])
		plt.xlabel('x (m)')
		plt.ylabel('y (m)')
		plt.legend(prop={'size':12})
		plt.show()
		if save:
			plt.savefig('sim_pairs_%d_pth_%.2f_gammath_%d'%(self.n, self.p_th, self.gamma_th),format='png')

	#Private Functions
	
	def _find_min_PWD(self, pi):
		#To use Gurboi, the objective must be linear or quadratic	
		m = gp.Model('min_pairwise_distance')
		#Let's keep this nice and quiet
		#m.Params.OutputFlag = 0
		#setup dummy variables required for linear objective
		n_regions = len(self.cregions)
		s = []
		s_obj = []
		for i in range(n_regions):
			for j in range(i+1, n_regions):
				s.append(m.addVar())
				s_obj.append(pi[i]*pi[j]*s[-1])	

		m.setObjective(np.sum(s_obj))
		xs = [m.addVar() for i in range(n_regions)]
		ys = [m.addVar() for i in range(n_regions)]
		
		#add constraints to couple new dummy vars to relay postitions
		ij = 0
		for i in range(n_regions):
			for j in range(i+1, n_regions): 
				#add the quadratic constraint that ||xi - xj||^2 = sij^2
				M = np.array([[1,0,0, 0, 0],[0,-1,1, 0, 0],[0,1,-1, 0, 0], [0, 0, 0, -1, 1], [0, 0, 0, 1, -1]])
				xc=[s[ij], xs[i], xs[j], ys[i], ys[j]]
				m.addMQConstr(M,None, GRB.EQUAL, 0, xc, xc)
				ij += 1
	
		#constrain relay points to lie within their regions
		eta = []
		for i in range(n_regions):
			reg = self.cregions[i]
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
		    
			eta_i = [m.addVar(vtype=GRB.BINARY) for k in range(n_vars-2)]
			eta.append(eta_i)
			m.addConstr(np.sum(eta_i) == 1)
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
		    
			LMC = m.addMConstr(A, [xs[i], ys[i] , *eta_i], GRB.LESS_EQUAL, b)
		    
		#and now solve
		m.params.NonConvex = 2
		m.optimize()
		
		#TODO - feasability checking/handling
		#assert m.status == 2, '
		
		#and extract the optimal values
		x = []
		for i in range(n_regions):
			x.append([xs[i].x, ys[i].x])

		#return the points as well as the optimal value
		return np.array(x), m.getObjective().getValue()
		
	def _find_min_PWDalt(self, pi):
		#use a series of piecewise linear constraints rather than binary variables
		
		
		#To use Gurboi, the objective must be linear or quadratic	
		m = gp.Model('min_pairwise_distance')
		#Let's keep this nice and quiet
		#m.Params.OutputFlag = 0
		#setup dummy variables required for linear objective
		n_regions = len(self.cregions)
		s = []
		s_obj = []
		for i in range(n_regions):
			for j in range(i+1, n_regions):
				s.append(m.addVar())
				s_obj.append(pi[i]*pi[j]*s[-1])	

		m.setObjective(np.sum(s_obj))
		xs = [m.addVar() for i in range(n_regions)]
		ys = [m.addVar() for i in range(n_regions)]
		
		#add constraints to couple new dummy vars to relay postitions
		ij = 0
		for i in range(n_regions):
			for j in range(i+1, n_regions): 
				#add the quadratic constraint that ||xi - xj||^2 = sij^2
				M = np.array([[1,0,0, 0, 0],[0,-1,1, 0, 0],[0,1,-1, 0, 0], [0, 0, 0, -1, 1], [0, 0, 0, 1, -1]])
				xc=[s[ij], xs[i], xs[j], ys[i], ys[j]]
				m.addMQConstr(M,None, GRB.EQUAL, 0, xc, xc)
				ij += 1
	
		#constrain relay points to lie within their regions
		for i in range(n_regions):
			reg = self.cregions[i]
			theta_tildes = []
			for poly in reg.polygons:
				for cnvx in poly.cnvx_partition:
					Aik,bik = cnvx.to_linear_constraints()
					#recall Aik [x,y]^T -b<0 if [x,y] in polygon
					n_edges = len(bik)
					zs = [m.addVar(lb=float('-inf')) for i in range(n_edges)]
					Aik_aug = np.zeros((n_edges,2+n_edges))
					Aik_aug[:,:2] = Aik
					Aik_aug[:,2:] = np.eye(n_edges)
					m.addMConstr(Aik_aug, [xs[i], ys[i], *zs ], GRB.EQUAL, bik)
					
					ztildes = [m.addVar(lb=float('-inf')) for i in range(n_edges)]
					for j in range(n_edges):
						#first of the piecewise constraints
						m.addGenConstrPWL(zs[j], ztildes[j], [-1,0,1],[-1000, 1, 1])
					
					#dummy variable for checking how many constraints are satisfied
					theta = m.addVar(lb=float('-inf'))
					a = np.ones((1, n_edges+1))
					a[-1] = -1
					m.addMConstr(a, [*ztildes, theta], GRB.EQUAL, np.array([0]))
					
					# saturate so that theta_tilde is 1 if all constraints satisfied,
					# ~0 otherwise
					theta_tilde = m.addVar()
					theta_tildes.append(theta_tilde)
					m.addGenConstrPWL(theta, theta_tilde, [0,n_edges-0.001,n_edges,n_edges+1],
								[0,0, 1, 1])
						
			#sum of theta tildes should be greater than 1,
			# i.e., at least all constraints for one cnvx region should be satisfied
			a = np.ones((1,len(theta_tildes)))
			m.addMConstr(a, theta_tildes, GRB.GREATER_EQUAL, np.array([1]))
					
		    
		#and now solve
		m.params.NonConvex = 2
		m.optimize()
		
		#TODO - feasability checking/handling
		#assert m.status == 2, '
		
		#and extract the optimal values
		x = []
		for i in range(n_regions):
			x.append([xs[i].x, ys[i].x])

		#return the points as well as the optimal value
		return np.array(x), m.getObjective().getValue()

	def _simple_gradient(self, X, S, pi, lr):
		Xnext = np.copy(X)
		for i in range(self.n):
			grad = self._grad_xi(i, X, S, pi)
			Xnext[i] = X[i] - lr*grad
			idx = toGridFromRaw(self.region, self.res, Xnext[i])
			#make sure we're still in the region
			if not self.cfields[i][idx[0][0], idx[1][0]]:
				Xnext[i] = X[i]

		return Xnext

	def _Metropolis(self, X, S, pi, max_it = 100, temp=100):
		rp = MRP.RandomRP(pi)
		Xnext = np.copy(X)
		shape = self.cfields[0].shape
		neighborhood = [[1,1],[1,0],[1,-1],[0,1],[0,0],[0,-1],[-1,1],[-1,0],[-1,-1]]
		for it in range(max_it):
			for i in range(self.n):
				idx = toGridFromRaw(self.region, self.res, Xnext[i]).reshape(2)
				scores = np.Inf*np.ones(9)
				for j in range(9):
					neighbor = idx + neighborhood[j]
					#index in region and is connected
					if (neighbor >= 0).all() and (neighbor < shape).all() and self.cfields[i][neighbor[0], neighbor[1]] :
						Xtemp = np.copy(Xnext)
						xitemp = toRawFromGrid(self.region, self.res, neighbor)
						Xtemp[i] = xitemp
						S =  XtoS(Xtemp)
						scores[j] = self.ps.calc_avg_wait(rp, S)
				#move to neighbors with some probability
				#normalize
				t = temp/(1 + (it)//2)
				scores = np.exp(scores)**(-1/t)
				scores = scores/np.sum(scores)
				neighbor_idx = np.random.choice(9,p=scores)
				neighbor = idx + neighborhood[neighbor_idx]
				Xnext[i] = toRawFromGrid(self.region, self.res, neighbor)
		#print(Xnext)
		return Xnext

	def _PSO(self, X, S, pi):
		return X

	def _pick_random_X(self):
		X = np.zeros((len(self.Xis), 2))
		for i  in range(len(self.Xis)):
			pts = self.Xis[i]['points']
			xi_idx = np.random.choice(len(pts))
			X[i] = pts[xi_idx]
		return X

	def _grad_xi(self, i, X, S, pi):
		dists = np.maximum(S[i,:], 1e-8).reshape((self.n, 1)) #avoid divide by 0, will zero out elsewhere
		partials = np.array([ self._partial_sij(i, j, S, pi) for j in range(self.n)]).reshape((1,self.n))
		weights = (1/dists)*(X[i] - X)
		return  partials @ weights

	def _partial_sij(self,i, j, S, pi):
		s = np.transpose(pi) @ S @ pi
		s_2 = np.transpose(pi) @ S**2 @ pi
		rhok = np.reshape(self.ps.Ls * self.ps.beta, (1, self.n))

		term1 = ( 2*(rhok - rhok**2) @ (1/pi)*(1-self.ps._Rho()) 
		 			+ self.ps._Rho()*(2*S[i,j]*s - s_2 )/s**2 )
		return term1*pi[i]*pi[j] - (pi[i]*rhok[0][j] + pi[j]*rhok[0][i])

def _lr(i):
	if i < 100:
		return 1
	elif i < 200:
		return 1e-1
	else:
		return 1e-3 

def XtoS(X, v = 1):
    n = len(X)
    S = np.zeros((n,n))
    for i in range(n):
        for j in range(i+1,n):
            dist = np.linalg.norm(X[i] - X[j])
            S[i,j] = dist/v
            S[j,i] = dist/v
    return S

def _indices_to_pts(r, region, res):
    #First get the non-zero indices
    idcs = np.array(np.where(r>0)).T
    pts = idcs - [region[1], region[3]]
    pts = pts/res
    return {'points': pts, 'indices': idcs}
