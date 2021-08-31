import numpy as np
def split_seg(p1, p2, npu=20):
	deltas = p2-p1
	length = np.linalg.norm(deltas)
	n = int(npu*length)
	path = np.zeros((n,2))
	path[0] = p1
	for i in range(1,n-1):
		path[i] = path[i-1]+(1/n)*deltas
	path[-1] = p2
	return path

def split_seg_by_time(p1, p2, v0, v1, t, dt=0.1):
	n = max(int(t//dt),2)
	a = (v1-v0)/t
	path = np.zeros((n,2))
	v = v0
	path[0] = p1
	for i in range(1,n-1):
		v_next = v+a*dt
		path[i] = path[i-1] + (v+v_next)*dt/2
		v = v_next
	path[-1] = p2
	return path
