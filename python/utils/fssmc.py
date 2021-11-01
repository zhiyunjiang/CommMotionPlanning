#useful functions for finite state space markov chains
import numpy as np

def stationary(P):
	pi = None
	if is_valid_transition_matrix(P):
		es, M = np.linalg.eig(P.T)
		evec1 = M[:,np.isclose(es, 1)]
		evec1 = evec1[:,0]
		pi = evec1/sum(evec1)
		pi = pi.real
	return pi

def is_valid_transition_matrix(P):
	return (np.sum(P,axis=1) == 1).all() and (P >=0).all()
	