from itertools import permutations
import numpy as np
import time


import sys 
from pathlib import Path
path = Path(__file__)
top_dir = path.parent.parent
sys.path.insert(0, str(top_dir.absolute())+"\\geometry")
import shot_solvers as SHOT
import gurobi_solvers as GB

def _unique_cycles(n):
	permutes = permutations(range(n))
	pruned_permutes = []
	for permute in permutes:
		if (permute[0] == 0) and not _flipped(permute) in pruned_permutes:
			pruned_permutes.append(permute)
	return pruned_permutes
	
def _flipped(perm):
	return (perm[0], *perm[len(perm):0:-1])


	
def lb(dm, p):
	n = len(p)
	bound = 0
	for i in range(n):
		bound += dm[p[i], p[(i+1)%n]]
	return bound
	

def TSPN_BF(regions):
	t = time.time()
	n = len(regions)

	#for each unique cycle (exclude congruent, directionally symmetric cycles)
	#solve the average weighted distance problem using Gurobi
	
	min_dists = np.zeros((n,n))
	for i in range(n):
		for j in range(i+1, n):
			d,_ = regions[i].distToPC(regions[j])
			min_dists[i,j] = d
			min_dists[j, i] = d
	bsf = np.inf
	argmin = None
	Ps = _unique_cycles(n)
	print('Total of %d Permutations to Try'%(len(Ps)))
	for P in Ps:
		print('Working on Permutation ' + str(P))
		if lb(min_dists, P) < bsf:
			# x, val = SHOT.min_cycle(regions, P)
			#will use the Gurboi solver, which seems to be more accurate
			x, val = GB.min_cycle(regions, P)
			#check if this is the best we've seen
			if val <= bsf:
				print('Optimal Solution Improved to %.2f'%(val))
				bsf = val
				argmin = (x, P)
		else:
			print("skiping based on lower bound")
	print("Elapsed Time: %.2f" %(time.time()-t))
	return bsf, argmin
	
	
#def TSPN_MIP
