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

#Delay-Tolerant Relay System
class DTR:

	#TODO eventually, Beta should be calculated from comm specifications
	def __init__(self, channels, ls, beta, th = -92):
		assert len(channels)%2 == 0, "Expected even number of channels"
		assert len(channels) == 2*len(ls), "Mismatch: number of regions and number of queue parameters"

		self.n = len(ls)
		self.region = channels[0].region
		self.res = channels[0].res
		self.cfs = [cc.getConnectionField(th) for cc in channels]
		self.cfields = [(i+1)*(self.cfs[2*i]*self.cfs[2*i+1]) for i in range(len(ls))]
		self.Xis = [_indices_to_pts(cf, self.region, self.res) for cf in self.cfields]
		self.cregions = [ PC.PointCloud(Xi['points']) for Xi in self.Xis]
		for pc in self.cregions:
    			pc.partition(algo=4)
		self.channels = channels
		self.ps = PS.PollingSystem(ls, beta)

	def optimize(self, do_plot=True, x_opt_method = 0):
		converged = False
		#initialize policy and polling locations
		pi = MRP.build_ed_policy(self.ps.Ls, self.ps.beta)
		#randomly sample points from Xis
		X = self._pick_random_X()
		Xs = [X]
		S = XtoS(np.copy(X))
		it_count = 0
		max_it = 500
		Wprev = 0
		while (not converged) and (it_count < max_it):
			#first, find optiaml pi
			res = self.ps.calc_optiaml_rp(S)
			pi = res.x
			W = res.fun
			if it_count%50 == 0:
				print("Optimized Policy Waiting Time: %.4f"%(W))
			#Now update the coordinates
			if x_opt_method == 0:
				X = self._simple_gradient(X, S, pi, _lr(it_count))
			elif x_opt_method == 1:
				X = self._PSO(X, S, pi)
			elif x_opt_method == 2:
				base_temp = 100
				X = self._Metropolis(X,S,pi, temp=100/(1 + it_count//25))
			Xs.append(np.copy(X))
			rp = MRP.RandomRP(pi)
			S =  XtoS(X)
			W = self.ps.calc_avg_wait(rp, S)
			if it_count%50 == 0:
				print("Optimized Location Waiting Time: %.4f"%(W))
			it_count += 1

			if (W - Wprev) == 0:
				converged = True
			Wprev = W
		if do_plot:
			self.plot_optimization(Xs)

		print("Optimized Waiting Time: %.4f"%(Wprev))
		return pi, X

	def plot_optimization(self, Xs):
		#plot the connectivity fields
		color_array=['white', 'red', 'green', 'cyan']
		cmap = ListedColormap(color_array)
		ticks =[0.75/2,1.5-0.75/2, 1.5+0.75/2, 3-0.75/2]
		tick_labels = ['No Connectivity','X1', 'X2', 'X3']
		CC.plotField(self.region, np.sum(self.cfields, axis=0), 'Regions of Joint Connectivity', cmap=cmap,
        			ticks=ticks, tick_labels = tick_labels, do_show=False)

		#now print the trajectories and final positions
		Xs = np.array(Xs)
		for i in range(self.n):
			xis = Xs[:,i,:]
			plt.plot(xis[:,0], xis[:,1], label='t%d'%(i+1))
			plt.scatter(xis[-1,0], xis[-1,1], marker="+", c='k', s=35)

		plt.legend()
		plt.show()

	#Private Functions

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
            S[i,j] = dist
            S[j,i] = dist
    return S

def _indices_to_pts(r, region, res):
    #First get the non-zero indices
    idcs = np.array(np.where(r>0)).T
    pts = idcs - [region[1], region[3]]
    pts = pts/res
    return {'points': pts, 'indices': idcs}
