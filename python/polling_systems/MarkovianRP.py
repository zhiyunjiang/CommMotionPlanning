import numpy as np

class MarkovianRP:

	def __init__(self, P):
		self.P = P # transition matrix
		self.n = len(P)

	def next(self,q, not_self=False):
		probs = self.P[q,:]# get all columns associated with q
		n = self.n
		if not_self:
			n = [i for i in range(self.n) if i != q]
			probs = np.array(probs[n])
			probs = probs/sum(probs)

		return np.random.choice(n, p=probs)

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

#build a policy which is strictly proportional to lambdas
def build_proportional_policy(ls):
	pi = ls/np.sum(ls)
	return RandomRP(pi)

def build_ed_policy(ls, beta):
	rhok = ls*beta
	pi = np.sqrt(rhok*(1-rhok))
	#normalize
	pi = pi/np.sum(pi)
	return RandomRP(pi)

def rnd2P(pi, S, alpha=0.5, verbose=False):
	
	import picos
	n = len(pi)
	#calculate true probability transitions
	P = np.zeros((n,n))
	for i in range(n):
	    for j in range(n):
	        if i != j:
	            P[i,j] = pi[j]
	    P[i,:] /= sum(P[i,:])
	v, M = np.linalg.eig(P.T)

	pi_obs = M[:,0]/sum(M[:,0])
	    
	model = picos.Problem()
	P = picos.RealVariable('P', (n,n))
	X = picos.RealVariable('X', (n,n))
	W = picos.Constant('W', S)
	v1 = picos.Constant('1_n', np.ones(n))
	q = picos.Constant('q', np.sqrt(pi_obs))
	Pi = picos.Constant('pi', pi_obs)

	PI_sqrt = picos.Constant('PI^.5', np.diag(np.sqrt(pi_obs)))
	PI_negsqrt = picos.Constant('PI^-.5', np.diag(1/np.sqrt(pi_obs)) )
	T = np.eye(n) - PI_sqrt*P*PI_negsqrt + q*q.T
	M = picos.block([[T, "I"],["I", X]])
	model.set_objective('min', picos.min(alpha*picos.trace(X) +(1-alpha)*Pi.T*(P^W)*v1))

	#and now some constraints
	#LMI linking slack variable X to P
	model.add_constraint(M.hermitianized >> 0)
	#row stochastic
	model.add_list_of_constraints([ sum(P[i,:]) == 1 for i in range(n)])
	#reversible MC
	model.add_list_of_constraints([ pi_obs[i]*P[i,j]==pi_obs[j]*P[j,i] for i in range(n) for j in range(i+1, n)])
	#entires of P must be no greater than 1...
	model.add_list_of_constraints([ P[i,j]<=1 for i in range(n) for j in range(n) if i!=j])
	#... and no less than 0 to be legitimate probability measure
	model.add_list_of_constraints([ 0<=P[i,j] for i in range(n) for j in range(n) if i!=j])
	#Avoid self loops, as these just create empty cycles
	model.add_list_of_constraints([P[i,i]==0 for i in range(n)])

	if verbose:
		print(model)

	solution = model.solve()

	Pnp = np.array(P.value)
	#clean up any numerical issues
	Pnp[Pnp<1e-6]=0
	for i in range(n):
		Pnp[i,:]/=np.sum(Pnp[i,:])

	return MarkovianRP(Pnp)
