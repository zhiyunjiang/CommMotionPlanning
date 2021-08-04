import numpy as np

class MarkovianRP:

	def __init__(self, P):
		self.P = P # transition matrix
		self.n = len(P)

	def next(self,q):
		probs = self.P[q,:]# get all columns associated with q
		return np.random.choice(self.n, p=probs)

	#STATIC FUNCTIONS#
	def _is_valid_dist(pi):
		return  ( np.sum(pi) == 1 and (pi >0).all() )

class RandomRP(MarkovianRP):

	def __init__(self, pi):
		#make sure it's in the right shape and that it's a numpy array
		self.pi = np.reshape( pi, (len(pi), 1) )

		#verify legitimate probaility distribution
		#assert MarkovianRP._is_valid_dist(self.pi), "Invalid probability distribution: " + str(self.pi)

		#create a P from the pi
		P = np.tile(pi, (len(pi), 1))
		super().__init__(P)

#build a policy which is strictly proportional to the 
def build_proportional_policy(ls):
	pi = ls/np.sum(ls)
	return RandomRP(pi)

def build_ed_policy(ls, beta):
	rhok = ls*beta
	pi = np.sqrt(rhok*(1-rhok))
	#normalize
	pi = pi/np.sum(pi)
	return RandomRP(pi)
