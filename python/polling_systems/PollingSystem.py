import numpy as np
from numpy.random import default_rng
from MarkovianRP import RandomRP
from StaticRP import CyclicRP, TableRP
from Queue import Queue
import warnings
import scipy.optimize
from scipy.optimize import minimize
import matplotlib.pyplot as plt

class PollingSystem:
	"""
	Currently assumes:
	 - exhaustive service discipline at all queues
	 - uniform, deterministic service times
	 - deterministic switching times
	"""
	def __init__(self, Ls, beta):
		self.n = len(Ls) #Number of queues in the system
		self.Ls = np.array(Ls) #Lambdas, i.e. arrival rates at each system
		self.beta = beta #service time
		if self.RhoSys() >= 1:
			warnings.warn("System with traffic %f will be unstable."%(self.RhoSys()))
			
	def LSys(self):
		"""
		Arrival rate for entire system
		"""
		return np.sum(self.Ls)

	def RhoSys(self):
		"""
		Traffic for entire system
		"""
		return self.beta * self.LSys()


	def calc_optimal_rp(self, S):
		#start out assuming equal distances
		x0 = self._opt_pi_equal_distances()
		sys_wait = lambda x: self._calc_avg_wait_random(S, x)
		n = len(x0)
		lb = 0.0005
		bounds = [(lb, 1-(n-1)*lb) for i in range(self.n)]
		A = np.ones(self.n)
		sum_to_one = scipy.optimize.LinearConstraint(A, 1, 1)
		res = minimize(sys_wait , x0, method='SLSQP', bounds = bounds, constraints = sum_to_one, tol=1e-6)
		argmin = res.x
		min_val = res.fun
		return argmin, min_val

	def _opt_pi_equal_distances(self):
		#assumes all exhaustive services
		rhos = self.Ls*self.beta
		pi = np.sqrt(rhos*(1-rhos))
		return pi/sum(pi)#normalize


	def calc_avg_wait(self, S, rp):
		"""
		Currently only handles random routing policy
		"""
		#TODO - handle Table policies once I can access
		#"Polling with a General-Service Order Table", (Baker, Rubin 1987)
		if isinstance(rp, RandomRP):
			return self._calc_avg_wait_random(S, rp.pi)
		elif isinstance(rp, CyclicRP):
			return self._calc_avg_wait_cyclic(S)
		else:
			print("Theoretical Waiting Time calculation only implemented for clyclic and random policies.\n")
			return -1
		
	def simulate(self, rp, S, tmax, q = 0):

		#initialization
		t = 0.0
		xt=[]
		wt=[]
		polling_instants = []
		is_traveling = False
		queues = [Queue() for i in range(self.n)]
		stage = 0
		total_travel_time = 0

		Tcounts = np.zeros((self.n, self.n))
		Tsums = np.zeros((self.n, self.n))
		Tlast_visit = np.zeros(self.n)#technically off, but shouldn't make a big deal
		while t <=tmax:
			x = [len(queue.waiting) for queue in queues]
			xt.append(np.concatenate( (np.array([t, q]), x) ))
			sys_avg_wait = self._calc_sim_stats(queues)
			wt.append([t,sys_avg_wait])

			if is_traveling:#we were traveling, now we're polling
				polling_instants.append([t,q])
			#check the number of requests currently at queue q
			reqs = x[q]
			if reqs > 0:
				#service those
				server_time = self.beta
				queues[q].start_service(t)
				is_traveling = False
			else:
				Tcounts[:,q]+=1
				Tsums[:,q]+= t - Tlast_visit
				Tlast_visit[q] = t

				#decide where we're going next
				q_prev = q
				q = rp.next(q)
				server_time = S[q_prev,q]
				total_travel_time += server_time
				is_traveling = True

			#account for anything that happens during the server time
			t_next = t + server_time

			#register arrivals that showed up during server_time
			xt = self._register_arrivals(server_time, t, queues, q, is_traveling, xt)

			t = t_next
			if is_traveling:
				stage += 1

		sys_avg_wait = self._calc_sim_stats(queues)
		wt.append([tmax, sys_avg_wait])
		S_bar = 0
		if stage>0:
			S_bar = total_travel_time/stage

		return xt, wt, queues, total_travel_time, S_bar, polling_instants, Tsums/Tcounts


	def waitMG1(self):
		term1 = 1/self.RhoSys()
		term2 = ( self.beta*self.RhoSys()**2 )/ ( 2*(1 - self.RhoSys()) )
		return term1 * term2


	def plotWvsPi(self, S):
		"""
		Currently only works for systems with 3 nodes
		"""
		if self.n == 3:
			gran = 100
			pi1 = np.linspace(0,1, gran)
			pi2 = np.linspace(0,1, gran)
			Z = np.zeros([gran, gran])
			for i in range(gran):
				for j in range(gran):
					Z[i,j] = self._calc_avg_wait_random( np.array([pi1[i], pi2[j], 1 - pi1[i] - pi2[j]]), S)
			#now plot this SOB
			plt.imshow(Z, extent=[0,1,0,1], cmap='hot', interpolation='nearest')
			plt.show()
		

	"""
	Private Methods
	"""
	def _calc_sim_stats(self, queues):
		sum_of_all_wait_times = 0
		n_serviced = 0
		for q in queues:
			sum_of_all_wait_times += sum(q.wait_times)
			n_serviced += len(q.wait_times)
		if n_serviced > 0:
			sys_avg_wait = sum_of_all_wait_times/n_serviced
		else:
			sys_avg_wait = 0

		return sys_avg_wait

	#Useful system properties

	def _S_avg(self, S, pi):
		return pi.T @ S @ pi

	def _t_avg(self, S, pi):
		return self._S_avg(S, pi)/(1-self.RhoSys())

	def _S_2(self, S, pi):
		return pi.T @ (S**2) @ pi		

	def _D_2(self, S, pi):
		#the mony ball
		return sum([pi[i]*self._var_D_i(S, pi, i) for i in range(self.n)]) + (self._t_avg(S, pi)*self.LSys())**2

	def _var_D_i(self, S, pi, i):
		tbar = self._t_avg(S, pi)
		rho_i = self._Rho_i(i)

		dbari = tbar*self.Ls[i]/pi[i]

		var_i = dbari *(1+ 2*(rho_i)/(1-rho_i) )

		return var_i

	def _Rho_i(self, i):
		return self.beta*self.Ls[i]

	def _Li_mc_avg_at_i(self, S, pi, i):
		return (self._t_avg(S, pi)*(1-self._Rho_i(i))*self.Ls[i])/pi[i] # true given at i


	def _Tij_avg(self, S, pi, i, j): #time since last visit to i given we're at j
		sj_bar = S[j,:]@pi
		rhos = self.beta * self.Ls
		tbar = self._t_avg(S, pi)
		mask = np.array([j for j in range(self.n) if j != i])
		s_no_i = (S[:, mask] @pi[mask])/(1.0-pi[i])


		return ((sj_bar + tbar*rhos[j]/pi[j]) + (1/pi[i])*tbar*(self.RhoSys() - rhos[i]) + 
		( (1-pi[i])/pi[i])*sum([pi[k]*s_no_i[k] for k in range(self.n) if k!= i])
		+ (1-pi[i])*s_no_i[i] )

	def _Tij_avg_alt(self, S, pi, i, j):
		sj_bar = S[j,:]@pi
		rhos = self.beta * self.Ls
		tbar = self._t_avg(S, pi)
		return (sj_bar + tbar*rhos[j]/pi[j]) + (1-pi[i])*tbar/pi[i]

	def _Li_mc_avg(self, S, pi, i):
		rhos = self.beta*self.Ls
		tbar = self._t_avg(S, pi)
		mask = np.array([j for j in range(self.n) if j != i])

		S_no_i = (S[:, mask] @pi[mask])/(1.0-pi[i])#normalize

		return tbar*self.Ls[i]*(1-rhos[i]) + (1-pi[i])*self.Ls[i]*(
			S_no_i[i] + (1/pi[i])*(sum([pi[k]*S_no_i[k] for k in range(self.n) if k != i ])
												 + (self.RhoSys() - rhos[i])*tbar)
			)

	def _LSys_mc_avg(self, S, pi):
		return sum([self._Li_mc_avg(S,pi, i) for i in range(self.n)])

	def _calc_avg_wait_random(self, S, pi):
		P = np.tile(np.reshape(pi,self.n ), (self.n, 1))
		T = np.array([ [ self._Tij_avg(S, pi, i, j) for j in range(self.n)] for i in range(self.n)])

		return self._calc_avg_wait_markovian(P, S, T, pi)

	def _calc_avg_wait_markovian(self, P, S, T, pi = None):
		"""
		Based on "Waiting Times in Polling Systems with Markovian Server Routing" (Boxma & Westerate, 1989)
		See also "Efficient Visit Orders for Polling Systems" (Boxma, Levy, & Westerate, 1993)
		"""
		if pi is None:
			v, M = np.linalg.eig(P.T)
			pi = M[:,0]/sum(M[:,0])
		sk = np.diag(P@S)
		sbar = pi.T @ sk
		s_2 = pi.T @ np.diag( P@(S**2))
		rhok = self.Ls * self.beta

		term1 = 1/self.RhoSys()
		term2 = ( self.beta*self.RhoSys()**2 )/ ( 2*(1 - self.RhoSys()) )
		
		term3 = (1/sbar)*(sum([pi[i]*sum([P[i,j]*S[i,j]*sum([rhok[k]*T[k,i] for k in range(self.n) if k!=i])  for j in range(self.n) ]) for i in range(self.n)]))

		term4 =  (self.RhoSys() * s_2) /(2*sbar )
		
		return np.reshape(term1 * ( term2 + term3 + term4 ), 1)[0]
		
	def _calc_avg_wait_random_uni_si(self, S, pi):

		"""
		Holds when s_ij = s_ik for all j, k (more generally, when first and second moment of s_ij's are 
		equal for any fixed i)

		Based on "Waiting Times in Polling Systems with Markovian Server Routing" (Boxma & Westerate, 1989)
		See also "Efficient Visit Orders for Polling Systems" (Boxma, Levy, & Westerate, 1993)
				 "The Analysis of Random Polling Systems" (Kleinrock & Levy, 1988)
		"""
		sk = S @ pi
		sk_2 = S**2 @ pi
		rhok = np.reshape(self.Ls * self.beta, (1, self.n))

		term1 = 1/self.RhoSys()
		term2 = ( self.beta*self.RhoSys()**2 )/ ( 2*(1 - self.RhoSys()) )
		term3 = (np.transpose(pi) @ sk)/(1-self.RhoSys()) * ( (rhok - rhok**2) @ (1/pi) )
		term4 = rhok @ sk
		term5 =  (self.RhoSys() * np.transpose(pi) @ sk_2) /(2* ( np.transpose(pi) @ sk) )
		
		return np.reshape(term1 * ( term2 + term3 - term4 + term5 ), 1)[0]

	def _cycle_times(self, S, pi):
		return self._t_avg(S, pi)*(1/pi)

	
	def _calc_avg_wait_cyclic(self, S):
		"""
		Based on "Workloads and Waiting Times in Single-Server Systems with Multiple Customer Classes" (Boxma, 1989)
		"""
		K_n = np.reshape(self.Ls/self.LSys(), (self.n, 1))
		K_nm = (K_n @ K_n.T ) - np.diag(K_n)
		s = np.sum(S)
		s_2 = s**2
		rhok = np.reshape(self.Ls * self.beta, (1, self.n))
		term1 = 1/self.RhoSys()
		term2 = ( self.beta*self.RhoSys()**2 )/ ( 2*(1 - self.RhoSys()) )
		term3 = (self.LSys()/( 2*(1 - self.RhoSys())))*(self.beta**2)*np.sum(K_nm)
		term4 = self.RhoSys()*s_2/(2*s)
		term5 = (s/ (2*(1 - self.RhoSys()) ))*(self.RhoSys()**2 - np.sum(rhok**2))
		return term1*(term2 + term3 + term4 + term5)
		
	def _calc_avg_wait_table(self):
		pass
	
	def _register_arrivals(self, server_time, t_start, queues, q_current, is_traveling, xt):
		rng = default_rng()
		if is_traveling:
			q_current = -1

		for q in range(self.n):
			t = 0
			while t <= server_time:
				#draw the next arrival from the exponential distribution
				t += rng.exponential(scale=1/self.Ls[q])
				if t<=server_time:
					queues[q].add(t_start+t)
		x = [len(queue.waiting) for queue in queues]
		xt.append(np.concatenate( ([t_start+server_time, q_current], x) ))
		return xt