import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

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
import shot_solvers as SHOT
import gurobi_solvers as GB
import cplex_solvers as CPLX

#Delay-Tolerant Relay System
class DTR:

	def __init__(self, channels, region, ls, beta, th = -92, p_th=0.7):
		assert len(channels)%2 == 0, "Expected even number of channels"
		assert len(channels) == 2*len(ls), "Mismatch: number of regions and number of queue parameters"

		self.n = len(ls)
		self.channels = channels
		self.p_th = p_th
		self.gamma_th = th
		self.region = region
		pfs = [cc.getPConField(self.gamma_th) for cc in self.channels]
		self.cfields = [1*((pfs[2*i]*pfs[2*i+1]) > self.p_th) for i in range(self.n)]
		self.Xis = [_indices_to_pts(self.cfields[i], channels[i*2].region, channels[i*2].res) for i in range(self.n)]
		
		self.cregions = [ PC.PointCloud(Xi['points']) for Xi in self.Xis]
		for i in range(self.n):
		
			self.cregions[i].partition(algo=4, alpha = 0.5/self.channels[i].res)

		self.ps = PS.PollingSystem(ls, beta)

	def shiftRegion(self, rid, offset):
		new_region = np.array(self.channels[rid*2].region) + np.array([offset[0], offset[0], offset[1], offset[1]])
		self.channels[rid*2].region = new_region
		self.channels[rid*2+1].region = new_region
		self.Xis[rid] = _indices_to_pts(self.cfields[rid], self.channels[rid*2].region, self.channels[rid*2].res) 
		
		self.cregions[rid] = PC.PointCloud(self.Xis[rid]['points']) 
		self.cregions[rid].partition(algo=4, alpha = 0.5/self.channels[rid].res)

	def optimize(self, do_plot=True, x_opt_method = 0, verbose = False, v=1, X0 = None):
		converged = False
		#initialize policy and polling locations
		pi = MRP.build_ed_policy(self.ps.Ls, self.ps.beta).pi
		if X0 is None:
			#randomly sample points from Xis
			X = self._pick_random_X()
		else:
			X = X0
		argmin = (X, pi)
		min_W = float('inf')
		Xs = [X]
		S = XtoS(np.copy(X), v)
		it_count = 0
		it_wo_improvement = 5
		max_it = 50
		max_it_wo_improvement = 10
		Wprev = np.inf

		while (not converged) and (it_count < max_it) and (it_wo_improvement < max_it_wo_improvement):
			#first, find optiaml pi
			pi, W = self.ps.calc_optimal_rp(S)
			#Now update the coordinates
			if x_opt_method == 0:
				X = self._simple_gradient(X, S, pi, _lr(it_count))
			elif x_opt_method == 1:
				X = self._PSO(X, S, pi)
			elif x_opt_method == 2:
				base_temp = 100
				X = self._Metropolis(X,S,pi, temp=100/(1 + it_count//25))
			elif x_opt_method ==3:
				#X, _ = SHOT.min_PWD(self.cregions, pi, 1, verbose)
				# X, _ = GB.min_PWD(self.cregions, pi, 1, verbose)
				X, _ = CPLX.min_PWD(self.cregions, pi)
			Xs.append(np.copy(X))
			rp = MRP.RandomRP(pi)
			S =  XtoS(X,v)
			W = self.ps.calc_avg_wait(S, rp)
			if verbose:
				print('Transition probabilities: ', pi)
				print('Points: ', X)
				print("Optimized Waiting Time: %.4f"%(W))
			it_count += 1
			it_wo_improvement += 1
			if  W < min_W:
				min_W = W
				argmin= (X, pi)
				it_wo_improvement = 0
			eps = 0.001
			if abs(W - Wprev) <= eps:
				converged = True
			Wprev = W
		if do_plot:
			self.plot_optimization(X)
		if verbose:
			print("Optimized Waiting Time: %.4f"%(Wprev))
		return min_W, argmin[1], argmin[0]

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
					

	def _simple_gradient(self, X, S, pi, lr):
		Xnext = np.copy(X)
		for i in range(self.n):
			grad = self._grad_xi(i, X, S, pi)
			Xnext[i] = X[i] - lr*grad
			idx = toGridFromRaw(self.channels[i].region, self.channels[i].res, Xnext[i])
			#make sure we're still in the region
			if not self.cfields[i][idx[0][0], idx[1][0]]:
				Xnext[i] = X[i]

		return Xnext

	def _Metropolis(self, X, S, pi, max_it = 100, temp=100):
		rp = MRP.RandomRP(pi)
		Xnext = np.copy(X)
		neighborhood = np.array([[1,1],[1,0],[1,-1],[0,1],[0,0],[0,-1],[-1,1],[-1,0],[-1,-1]])
		for it in range(max_it):
			for i in range(self.n):
				idx = toGridFromRaw(self.channels[i*2].region, self.channels[i*2].res, Xnext[i]).reshape(2)
				scores = np.Inf*np.ones(9)
				shape = self.cfields[i].shape
				for j in range(9):
					neighbor = idx + neighborhood[j]
					#index in region and is connected
					if (neighbor >= 0).all() and (neighbor < shape).all() and self.cfields[i][neighbor[0], neighbor[1]] :
						Xtemp = np.copy(Xnext)
						xitemp = toRawFromGrid(self.channels[i].region, self.channels[i].res, neighbor)
						Xtemp[i] = xitemp
						S =  XtoS(Xtemp)
						scores[j] = self.ps.calc_avg_wait(S, rp)
				#move to neighbors with some probability
				#normalize
				t = temp/(1 + (it)//2)
				scores = np.exp(scores)**(-1/t)
				scores = scores/np.sum(scores)
				neighbor_idx = np.random.choice(9,p=scores)
				neighbor = idx + neighborhood[neighbor_idx]
				Xnext[i] = toRawFromGrid(self.channels[i].region, self.channels[i].res, neighbor)
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

		term1 = ( 2*(rhok - rhok**2) @ (1/pi)*(1-self.ps.RhoSys()) 
		 			+ self.ps.RhoSys()*(2*S[i,j]*s - s_2 )/s**2 )
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
    pts = toRawFromGrid(region, res, idcs)
    return {'points': pts, 'indices': idcs}
