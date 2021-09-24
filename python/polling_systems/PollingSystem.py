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
		#start out assuming frequencies directly proportional to arrival rates
		x0 = self.Ls/self.LSys()
		sys_wait = lambda x: self._calc_avg_wait_random(x, S)
		bounds = [(0.001,1) for i in range(self.n)]
		A = np.ones(self.n)
		constraint = scipy.optimize.LinearConstraint(A, 1, 1)
		return minimize(sys_wait , x0, method='SLSQP', bounds = bounds, constraints = constraint, tol=1e-6)


	def calc_avg_wait(self, rp, S):
		"""
		Currently only handles random routing policy
		"""
		#TODO - handle Table policies once I can access
		#"Polling with a General-Service Order Table", (Baker, Rubin 1987)
		if isinstance(rp, RandomRP):
			return self._calc_avg_wait_random(rp.pi, S)
		elif isinstance(rp, CyclicRP):
			return self._calc_avg_wait_cyclic(S)
		else:
			print("Theoretical Waiting Time calculation only implemented for clyclic and random policies.\n")
			return -1
		
	def simulate(self, rp, S, tmax, q = 0):

		#generate arrival times at each queue
		arrival_times = self._generate_arrivals(tmax)
		t = 0.0
		xt=[]
		wt=[]
		is_traveling = False
		queues = [Queue() for i in range(self.n)]
		stage = 0
		total_travel_time = 0
		while t <=tmax:
			x = [len(queue.waiting) for queue in queues]
			xt.append(np.concatenate( (np.array([t, q]), x) ))
			sys_avg_wait = self._calc_sim_stats(queues)
			wt.append([t,sys_avg_wait])
			#check the number of requests currently at queue q
			reqs = x[q]
			if reqs > 0:
				#service those
				server_time = self.beta
				queues[q].start_service(t)
				is_traveling = False
			else:
				#
				q_prev = q
				while q == q_prev:#ignore empty cycles
					q = rp.next(q)
					stage +=1
				server_time = S[q_prev,q]
				total_travel_time += server_time
				is_traveling = True
			#account for anything that happens during the server time
			t_next = t + server_time
			if len(arrival_times) > 0:
				xt, arrival_times = self._register_arrivals(arrival_times, t_next, queues, q, is_traveling, xt)

			#register arrivals that showed up during server_time
			#xt = self._register_arrivals_2(server_time, t, queues, q, is_traveling, xt)

			t = t_next

		sys_avg_wait = self._calc_sim_stats(queues)
		wt.append([tmax, sys_avg_wait])
		S_tilde = 0
		if stage>0:
			S_tilde = total_travel_time/stage
		return xt, wt, queues, S_tilde


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
		
	def _calc_avg_wait_random(self, pi, S):

		"""
		Based on "Efficient Visit Orders for Polling Systems" (Boxma, Levy, & Westerate, 1993)
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

	def _calc_avg_wait_markovian(self, P, S):
		pass
		
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
	
	def _register_arrivals_2(self, server_time, t_start, queues, q_current, is_traveling, xt):
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

	def _register_arrivals(self, arrival_times, t_next, queues, q, is_traveling, xt):
		t_arrival = arrival_times[0][0]
		while t_arrival <= t_next:
			arrival = arrival_times[0]

			#add to the correct queue
			queues[arrival[1]].add(arrival[0])

			q_at_arv = q
			if is_traveling:
				q_at_arv = -1

			x = [len(queue.waiting) for queue in queues]
			xt.append(np.concatenate( ([arrival[0], q_at_arv], x) ))
			
			if len(arrival_times)>=2:
				arrival_times = arrival_times[1:]
				t_arrival = arrival_times[0][0]
			else:
				#we just serviced the only arrival
				arrival_times = []
				break
		return xt, arrival_times

	def _generate_arrivals(self, tmax):

		arrival_times = []
		rng = default_rng()
		for q in range(self.n):
			t = 0
			while t <= tmax:
				#draw the next arrival from the exponential distribution
				t += rng.exponential(scale=1/self.Ls[q])
				if t<=tmax:
					arrival_times.append((t,q))

		arrival_times = np.array(arrival_times, dtype=[('time_inst',float), ('queue', int)])
		return np.sort(arrival_times, order='time_inst')
				
