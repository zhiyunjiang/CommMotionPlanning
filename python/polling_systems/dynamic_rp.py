import numpy as np

def nextL(q_lengths, lambdas, beta, S, q, q_next):
	switch_time = S[q, q_next]
	lam_sys = np.sum(lambdas)
	return np.sum(q_lengths) + lam_sys*switch_time + (lam_sys - lambdas[q_next])*( (q_lengths[q_next] + lambdas[q_next]*switch_time)/(1 - lambdas[q_next]*beta)) - q_lengths[q_next]


def shortest_lengths(q, q_lengths, lambdas, beta, S):
	n = len(q_lengths)
	nextLs = np.array([nextL(q_lengths, lambdas, beta, S, q, i) for i in range(n)])
	nextLs[q] = float('inf')
	return np.argmin(nextLs)


class DPRP:

	def __init__(self, lambdas, beta, S, fmap, fully_observed = True):
		self.lambdas = lambdas
		self.beta = beta
		self.S = S
		self.fmap = fmap
		self.fully_observed = fully_observed

	def next(self, q, x):
		q_lengths = x
		if not self.fully_observed:
			q_lengths = x*self.lambdas

		return self.fmap(q, q_lengths, self.lambdas, self.beta, self.S)