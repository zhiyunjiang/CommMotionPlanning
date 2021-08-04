"""
CYCLICAL AND ROUTING TABLE POLLING SYSTEM ROUTING POLICIES
"""
import numpy as np


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
		
class CyclicRP(StaticRP):
	pass
	
class TableRP(StaticRP):
	pass


def SRPInOrder(n):
	sequence = [i for i in range(n)]
	return CyclicRP(sequence)
	
def SRPFromPis(pis, eps=0.01):
	pis = np.array(pis)
	length = len(pis)
	expected_counts = length*np.array(pis)
	while not _expected_counts_OK(expected_counts, eps):
		length += 1
		expected_counts = length*np.array(pis)
		
	counts = np.rint(expected_counts)
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
	return TableRP(sequence)
	
def _expected_counts_OK(ecnts, eps):
	OK = True
	for cnt in ecnts:
		if np.rint(cnt) ==0 or ( (cnt%1) > eps and (cnt%1) < 1-eps):
			OK = False
			break
	return OK
	

