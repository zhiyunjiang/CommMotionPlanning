"""
CYCLICAL AND ROUTING TABLE POLLING SYSTEM ROUTING POLICIES
"""
import numpy as np
from collections import Counter


class StaticRP:

	def __init__(self, sequence):
		self.n = len(sequence)
		self.seq = sequence
		self._cidx = -1
	
	#q parameter kept for compatibility with non-deterministic routing policies
	def next(self, q=None):
		self._cidx = (self._cidx + 1)%self.n
		return self.seq[self._cidx]
		
	def reset(self):
		self._cidx = -1

	def __eq__(self, obj):
		if ( not isinstance(obj, StaticRP) ) or (self.n != obj.n):
			return False

		ref = self.seq[0]
		#find all ocurrences of ref in other, see if any match up
		idxs = [i for i,x in enumerate(obj.seq) if x == ref]
		match = False
		for idx in idxs:

			match = True
			for i in range(self.n):
				if self.seq[i] != obj.seq[(i+idx)%self.n]:
					match = False
					break
			if match:
				break
		return match

class CyclicRP(StaticRP):
	pass
	
class TableRP(StaticRP):
	pass


def SRPInOrder(n):
	sequence = [i for i in range(n)]
	return CyclicRP(sequence)
	
def SRPFromPis(pis, eps=0.01, use_obs = False):
	sequence = _sequence_from_pi(pis, eps, use_obs)

	return TableRP(sequence)


def MinSSRPFromPis(pis, S, eps=0.01):
	sequence = _sequence_from_pi(pis, eps, use_pi_obs = True)
	optimizer = SeqDPOptimizer(sequence, S)
	sequence, _ = optimizer.optimize()
	return TableRP(sequence)

def _sequence_from_pi(pis, eps, use_pi_obs = False):
	pis = np.array(pis)
	length = len(pis)
	if use_pi_obs:
		n = length
		P = np.zeros((n,n))
		for i in range(n):
		    for j in range(n):
		        if i != j:
		            P[i,j] = pis[j]
		    P[i,:] /= sum(P[i,:])
		v, M = np.linalg.eig(P.T)

		pis = M[:,0]/sum(M[:,0]) 
		
	expected_counts = length*np.array(pis)
	l_max = 1000
	while not _expected_counts_OK(expected_counts, eps) and length <= l_max:
		length += 1
		expected_counts = length*np.array(pis)
		
	counts = np.rint(expected_counts)
	length = int(np.sum(counts))
	#now apply the Golden ratio
	phi_inv = 0.5*(np.sqrt(5) - 1)
	phi_mods = np.array([((i+1)*phi_inv)%1 for i in range(length)])
	mods_sorted = np.sort(phi_mods)
	sequence = [-1 for i in range(length)]
	
	running_count = 0
	for i in range(len(counts)):
		count = int(counts[i])
		for j in range(count):
			idx = np.where(mods_sorted == phi_mods[running_count + j] )
			sequence[idx[0][0]] = i
		running_count += count

	return sequence

	
def _expected_counts_OK(ecnts, eps):
	OK = True
	for cnt in ecnts:
		if np.rint(cnt) == 0 or ( (cnt%1) > eps and (cnt%1) < 1-eps):
			OK = False
			break
	return OK

class SeqDPOptimizer:

	def __init__(self, sequence, S):
		cntr = Counter(sequence)
		self.opts = cntr.keys()
		self.n =  len(self.opts)
		self.m = len(sequence)
		self.spacings = {}
		for opt in self.opts:
			idxs = [i for i,x in enumerate(sequence) if x == opt]
			spacing = [idxs[j+1]-idxs[j] for j in range(len(idxs) - 1)]
			spacing.append(idxs[0]+self.m-idxs[-1])
			self.spacings[opt]=spacing

		self.optimal_cost2go = {}
		self.S = S

	def optimize(self):
		state = [-1 for i in range(self.m)]
		#print('All opts:')
		#print(self.opts)
		min_wait, argmin = self._optimize_recur(state, self.opts)
		return self._argmin_to_sequence(argmin), min_wait

	def _argmin_to_sequence(self, argmin):
		state = [-1 for i in range(self.m)]
		# print(argmin)
		for arg in argmin:
			v,c,state = self._tx(state, arg[0], arg[1])
		return state

	def _optimize_recur(self, state, opts):
		# print('Optimizing with State: ')
		# print(state)
		state_key = SeqDPOptimizer._state_to_str(state)
		if state_key in self.optimal_cost2go:
			optimal = self.optimal_cost2go[state_key]
			return optimal[0], optimal[1]

		min_2go = float('inf')
		argmin = [(-1, -1)]
		for opt in opts:
			for i in range(len(self.spacings[opt])):
				is_valid, stage_cost, nxt_state = self._tx(state, opt, i)
				if is_valid:
					nxt_opts = [o for o in opts if o != opt]

					#handle base case
					cost_2_go = 0
					argmin_2go = []
					if len(nxt_opts) != 0:
						# print('had available actions')
						cost_2_go, argmin_2go = self._optimize_recur(nxt_state, nxt_opts)
					cost_from_here = cost_2_go + stage_cost

					if cost_from_here < min_2go:
						# print('Updating min2go and argmin:')
						# print(min_2go)
						# print(argmin)
						min_2go = cost_from_here
						argmin = [(opt,i)] + argmin_2go

		self.optimal_cost2go[state_key] = (min_2go, argmin)
		return min_2go, argmin

	def _state_to_str(state):
		return ','.join(str(x) for x in state)

	def _tx(self, state, opt, i_spc):
		valid_transition = True
		
		nxt_state = [s for s in state]#create deep copy of the state
		
		cost = 0

		i_start = nxt_state.index(-1)
		nxt_state[i_start] = opt
		
		spacing = self.spacings[opt]
		n = len(spacing)
		for j in range(n):
			i_prev = (i_start -1) % self.m
			if nxt_state[i_prev] != -1:
				cost += self.S[opt, state[i_prev]]
			i_nxt = (i_start +1) % self.m
			if nxt_state[i_nxt] != -1:
				cost += self.S[opt, state[i_nxt]]
			nxt_state[i_start] = opt
			spc = spacing[(j+i_spc)%n]
			i_start = (i_start + spc) % self.m

			if ((nxt_state[i_start] != -1) and j <(n-1) ) or ( (j==n-1) and (nxt_state[i_start] != opt) ):
				valid_transition = False
				break

		return valid_transition, cost, nxt_state