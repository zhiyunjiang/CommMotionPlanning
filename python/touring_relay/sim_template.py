import numpy as np
import matplotlib
#allow for latex markup in matplotlib figures
matplotlib.rcParams['text.usetex'] = False
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap


#Import a few utility functions...
import sys  
from pathlib import Path
sys.path.insert(0, "../comm_channel")
sys.path.insert(0, "../polling_systems")
sys.path.insert(0, "../geometry")
sys.path.insert(0, "../utils")

#So we can import my local libs
import CommChannel as CC
import qos
import pointcloud as PC
import PollingSystem as PS
import MarkovianRP as MRP
import StaticRP as SRP
import dtr
from CoordTransforms import toRawFromGrid
from CoordTransforms import toGridFromRaw
fs = 50
ms = 50
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"],
    'font.size': fs})
    

COLORS = ['xkcd:aqua', 'xkcd:coral', 'xkcd:wheat', 'xkcd:green', 'xkcd:orange', 'xkcd:azure', 'xkcd:yellow', 'xkcd:cyan']
def create_channels(cps, region, res, GAMMA_TH, sub_regions = None):
	
	if sub_regions is None:
		sub_regions = [region for i in range(len(cps)//2)]
	ccs = [CC.CommChannel(cps[i], sub_regions[(i)//2], res) for i in range(len(cps))]
	
	for cc in ccs:
		cc.generateSH();cc.generateMP(2)#using GP model for channel sim
		
	cfs = [cc.getConnectionField(GAMMA_TH) for cc in ccs]
	true_joint_con_fields = [1*(cfs[i*2]*cfs[(i*2)+1]) for i in range(len(ccs)//2)]
	true_joint_con_pts = [field_to_pts(true_joint_con_fields[i], sub_regions[i], res) for i in range(len(ccs)//2)]
	
	return ccs, cfs, true_joint_con_fields, true_joint_con_pts
	
def predict_channels(res, ccs, true_joint_con_fields, GAMMA_TH, p_th=0.7, jcm = 0):
	sres = res#//2
	pct_sample = 0.01#1% sampling

	pcs = []
	for cc in ccs:
		region = cc.region
		n_samples = int(pct_sample*(region[0] - region[1])*(region[2] - region[3])*sres**2)
		print('Drawing %d samples from Channel %d'%(n_samples, len(pcs)+1))
		xs, vals = cc.sampleChannel(n_samples)
		pcs.append(CC.PredictedChannel(cc.cp, region, sres, xs, vals))
		pcs[-1].setPth(p_th)
		print('Completed PredictedChannel %d'%(len(pcs)))
		
	#Generate probability of connectivity fields
	pfs = [pc.getPConField(GAMMA_TH) for pc in pcs]
	
	if jcm == 0:
		pred_joint_con_fields = [1*(pfs[i*2]*pfs[(i*2)+1]>p_th) for i in range(len(ccs)//2)]
	elif jcm == 1:
		pred_joint_con_fields = [1*((pfs[i*2]>p_th)*(pfs[(i*2)+1]>p_th)) for i in range(len(ccs)//2)]
	pred_joint_con_pts = [ field_to_pts(pred_joint_con_fields[i], ccs[i*2].region, sres) for i in range(len(ccs)//2)]
	
	#calculate the probabilty of the point being in the predicted regions
	p_pred_connected = []
	for i in range(len(ccs)//2):
		#convert pred points back to grid points, but in res resolution
		idxs = toGridFromRaw(ccs[(i*2)].region, res, pred_joint_con_pts[i])
		p_pred_connected.append(prob_pred_in_true(true_joint_con_fields[i], idxs))
		
	return pcs, pfs, pred_joint_con_pts, p_pred_connected
	
def plotDecompFig(n, tjcps, pfs, qBase, region, pt_cld):
	fig = plt.figure(figsize=(22,18))
	#colors = ['xkcd:aqua', 'xkcd:coral', 'xkcd:wheat', 'xkcd:green']
	Z = pfs[0] * pfs[1]
	nx, ny = Z.shape
	X = np.linspace(region[1], region[0], nx)
	Y = np.linspace(region[3], region[2], ny)
	pts = tjcps[0]
	color=COLORS[1]
	plt.plot(pts[:,0], pts[:,1], 's', color=color , markeredgewidth=0.0, markersize=5, alpha=0.35, zorder=-10)
	#dummy series for label
	plt.plot([-1000],  [-1000], 's', color=color, markeredgewidth=0.0, markersize=20, alpha=0.35, label='True Relay Region')

	C = plt.contourf(X, Y, Z.T,[0.5, 0.6, 0.7, 0.8, 0.9, 1 ], alpha = 0.5, colors = ['xkcd:aqua', 'xkcd:wheat', 'xkcd:green', 'xkcd:yellow', 'xkcd:azure', ])
	# C = plt.contour(X, Y, Z.T,[0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 ])
	cb = fig.colorbar(C)
	cb.set_label('Predicted Probability of Connecctivity', labelpad=20)
	#cb.legend()

	for poly in pt_cld.polygons:
		poly.plot()
	#dummy series for label
	plt.plot( [-1000, -1000],  [-1000, -900], '-', color='k', markersize=20, label='Convex Partition ($p_{{th}}=0.7$)')

	#plot base stations
	plt.scatter([qBase[0][0]], [qBase[0][1]], color=color, marker='v', s=400, edgecolor='k', label='Source')
	plt.scatter([qBase[1][0]], [qBase[1][1]], color=color, marker='^', s=400, edgecolor='k', label='Destination')	    

	plt.xlim(region[1],region[0])
	plt.ylim(region[2],region[3])
	plt.xlabel('x (m)')
	plt.ylabel('y (m)')
	plt.legend(loc='upper center', bbox_to_anchor=[0.5,-.075], ncol=2)
	

def plot_bs(qBase, ax):
	#plot base stations
	n = len(qBase)//2
	for i in range(n):
		ax.scatter([qBase[2*i][0]], [qBase[2*i][1]],
		color=COLORS[i], marker='v', s=200, edgecolor='k', zorder=200)
		ax.scatter([qBase[2*i+1][0]], [qBase[2*i+1][1]],
		color=COLORS[i], marker='^', s=200, edgecolor='k', zorder=200)

	#dummy series for legend formatting
	ax.scatter([-100], [-100], marker='v', s=10*ms, color='w', edgecolor='k', label='Source')
	ax.scatter([-100], [-100], marker='^', s=10*ms, color='w', edgecolor='k', label='Destination')

def plot_relay_regions(n, tjcps, pjcps, els, pi, ax):
	for i in range(n):
		#plot the true field
		pts = tjcps[i]
		ax.plot(pts[:,0],  pts[:,1], '.', color=COLORS[i])

		#plot the predicted field
		pts = pjcps[i]
		ax.plot(pts[:,0],  pts[:,1], '.', color='k', alpha=0.1)
		#dummy series for better legend formatting
		s_label = 'Relay Region %d'%(i+1)
		if els is not None:
			s_label += ' ($\\lambda_%d = %.2f$, $\\tilde{\\pi}_%d = %.2f$)'%(i+1, els[i], i+1, pi[i])

		ax.plot([-100], [-100], '.', color=COLORS[i], markersize=ms, label=s_label)

	#dummy series for label
	ax.plot([-1000],  [-1000], '.', color='k', markersize=ms, alpha=0.25, label='Predicted Relay Regions')

def set_lims(region, ax):
	ax.set_xlim(region[1],region[0])
	ax.set_ylim(region[2],region[3])

def label_axes(ax):
	ax.set_xlabel('x (m)')
	ax.set_ylabel('y (m)')

def plotCFwithOverlay(n, tjcps, pjcps, qBase, region, els = None, pi = None, ax = None):
	if (ax is None):
		fig = plt.figure(figsize=(15,15))
		ax = fig.add_axes([0,0,1,1])
	
	plot_relay_regions(n, tjcps, pjcps, els, pi, ax)

	plot_bs(qBase, ax)

	set_lims(region, ax)

	label_axes(ax)
	
	

def calc_AORP(dt_sys, vel):
	W, pi, X = dt_sys.optimize(x_opt_method=3, do_plot = False, verbose=False, v=vel)
	AORP = {'WT': W, 'X': X, 'pi': pi}
	return AORP

def plot_regional_decomposition(dt_sys, tjcps, pjcps, qBase, region, els=None, pi = None, ax = plt):
	plotCFwithOverlay(dt_sys.n, tjcps, pjcps, qBase, region, els, pi, ax)
	
	for creg in dt_sys.cregions:
	    	creg.plot_polys(ax = ax)

def plot_AORP_W_TSPN(dt_sys, AORP, TSPNP, tjcps, pjcps, qBase, region, els, pi):
	fig, (ax1, ax2) = plt.subplots(1,2, figsize=(32, 16))

	#Plot the AORP
	plot_regional_decomposition(dt_sys, tjcps, pjcps, qBase, region, els, pi, ax = ax2)
	X=AORP['X']
	_plot_relay_points(X, ax2)
	pi = AORP['pi']
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

	lws = [pi_obs[i]*P[i,j] for i in range(n) for j in range(i+1, n)]

	base_width = 8/(2*np.max(lws))
	min_width=0.2
	for i in range(n):
		for j in range(i+1, n):
			ax2.plot(X[[i,j], 0], X[[i,j], 1], 'k', linewidth = max(2*base_width*pi_obs[i]*P[i,j], min_width) )

	#dummy series for all edges
	ax2.plot([-1000, -900], [-1000, -1000], 'k', label='Routes', markersize = ms)
	ax2.invert_yaxis()
	ax2.set_title('AORP')


	#plot the TSPNP
	plot_regional_decomposition(dt_sys, tjcps, pjcps, qBase, region, els, pi_obs, ax = ax1)
	X=TSPNP['X']
	order = list(TSPNP['SEQ'])
	order.append(0)#complete the loop
	ax1.plot(X[order,0], X[order,1], 'k', linewidth = base_width/n, zorder=100, label='Route')
	_plot_relay_points(X, ax1)
	ax1.invert_yaxis()
	ax1.set_title('TSPNP')

	handles, labels = ax1.get_legend_handles_labels()
	fig.legend(handles, labels, fontsize=fs, loc='upper center', bbox_to_anchor=[0.5,0.025], ncol=2)
	#fig.subplots_adjust(wspace=0.12)

def pi_to_P(pi):
	n = len(pi)

	#calculate true probability transitions
	P = np.zeros((n,n))
	for i in range(n):
	    for j in range(n):
	        if i != j:
	            P[i,j] = pi[j]
	    P[i,:] /= sum(P[i,:])
	v, M = np.linalg.eig(P.T)

	pi_tilde = M[:,0]/sum(M[:,0]) 
	return P, pi_tilde

def plot_AORP_routes(AORP, ax):
	X=AORP['X']
	_plot_relay_points(X, ax)
	pi = AORP['pi']
	n = len(pi)

	P, pi_obs = pi_to_P(pi) 

	lws = [pi_obs[i]*P[i,j] for i in range(n) for j in range(i+1, n)]

	base_width = 8/(2*np.max(lws))
	min_width=0.2
	for i in range(n):
		for j in range(i+1, n):
			ax.plot(X[[i,j], 0], X[[i,j], 1], 'k', linewidth = max(2*base_width*pi_obs[i]*P[i,j], min_width) )

	#dummy series for all edges
	ax.plot([-100, -90], [-100, -100], 'k', label='Routes')

	return lws, base_width

def plot_AORP(dt_sys, AORP, tjcps, pjcps, qBase, region, els, pi, ax = plt):
	plot_regional_decomposition(dt_sys, tjcps, pjcps, qBase, region, els, pi, ax)

	lws, base_width = plot_AORP_routes(AORP, ax)

	#_format_legend()
	ax.invert_yaxis()
	return lws, base_width

def plot_TSPNP(dt_sys, TSPNP, tjcps, pjcps, qBase, region, base_width = 8):
	plot_regional_decomposition(dt_sys, tjcps, pjcps, qBase, region)

	
	X=TSPNP['X']
	order = list(TSPNP['SEQ'])
	order.append(0)#complete the loop
	plt.plot(X[order,0], X[order,1], 'k', linewidth = base_width, zorder=100, label='Route')
	_plot_relay_points(X)
	_format_legend()
	plt.gca().invert_yaxis()


def plot_els(els):
	x = np.arrange(len(els))
	plt.bar(x, els)
def plot_pis(pis):
	x = np.arrange(len(pi))
	plt.bar(x, pis)

def plot_els_w_pis(els, pis):

	plt.rcParams.update({'font.size': 22})
	x = np.arange(len(els))
	width = 0.35
	fig, ax = plt.subplots()
	
	s1 = ax.bar(x-width/2, np.around(els, decimals=2), width, label='Norm. Arriv. Rt.')
	s2 = ax.bar(x+width/2, np.around(pis, decimals=2), width, label='Visit Freq.')
	
	ax.set_xticks(x)
	labels = ['Q '+str(i+1) for i in range(len(pis))]
	ax.set_xticklabels(labels)
	ax.legend()
	
	ax.bar_label(s1, padding = 3)
	ax.bar_label(s2, padding = 3)
	ax.figure.set_size_inches(10,8)
	plt.ylim(0,1.1)
	#fig.tight_layout()
	
def run_sims(ps, AORP, TSPNP, hrs, mins, seconds, motion_power, tx_power, v=1):
	from importlib import reload 
	reload(SRP)
	minutes = hrs*60+mins
	seconds = minutes*60 + seconds
	#look at policy-invariant values
	print("Theotretical MB serviced: " + str(ps.LSys()*seconds))
	AP = ps.RhoSys()*tx_power + (1-ps.RhoSys())*motion_power
	print("Theoretical Energy Consumption (J): " + str(AP*seconds))

	print('\tTh. WT\tWT\tE (J)\tMBS\tMBR')	
	#and now run the sims
	S =  dtr.XtoS(AORP['X'],v)
	pi = AORP['pi']
	#now run a bunch of simulations and look at the averages

	#first run for AORP
	aorp = MRP.RandomRP(pi)
	AORP_res, AORP_xt = runsimsforpolicy(ps, aorp, S, motion_power, tx_power, seconds)
	print('AORP\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f'%(AORP['WT'], AORP_res['WT'], AORP_res['E'], AORP_res['MBS'], AORP_res['MBR']))

	# #also look at what happens if we try the Markoviian routing policy
	# mrp = MRP.rnd2P(pi, S, alpha=0.95)
	# MRP_res, MRP_xt = runsimsforpolicy(ps, mrp, S, motion_power, tx_power, seconds)
	# print('MRP\t---\t%.2f\t%.2f\t%.2f\t%.2f'%(MRP_res['WT'], MRP_res['E'], MRP_res['MBS'], MRP_res['MBR'] ))

	rtable = SRP.SRPFromPis(pi, eps=0.01)
	rtable_res, rtable_xt = runsimsforpolicy(ps, rtable, S, motion_power, tx_power, seconds)
	print('Tab\t---\t%.2f\t%.2f\t%.2f\t%.2f'%(rtable_res['WT'], rtable_res['E'], rtable_res['MBS'], rtable_res['MBR']))

	# rtable = SRP.SRPFromPis(pi, eps=0.01, use_obs = True)
	# rtable_res, rtable_xt = runsimsforpolicy(ps, rtable, S, motion_power, tx_power, seconds)
	# print('Tab-O\t---\t%.2f\t%.2f\t%.2f\t%.2f'%(rtable_res['WT'], rtable_res['E'], rtable_res['MBS'], rtable_res['MBR']))

	# ortable = SRP.MinSSRPFromPis(pi, S, eps=0.1)
	# if rtable != ortable:
	# 	ortable_res, ortable_xt = runsimsforpolicy(ps, ortable, S, motion_power, tx_power, seconds)
	# 	print('OTab\t---\t%.2f\t%.2f\t%.2f\t%.2f'%(ortable_res['WT'], ortable_res['E'], ortable_res['MBS'], ortable_res['MBR']))	
	# #print(ortable.seq)
	#print(rtable.seq)
	#Run baseline (TSPN) policy
	S_TSPN = dtr.XtoS(TSPNP['X'],v)
	tspnp = SRP.StaticRP(list(TSPNP['SEQ']))
	tspn_res, tspn_xt = runsimsforpolicy(ps, tspnp, S_TSPN, motion_power, tx_power, seconds)
	print('TSPN\t---\t%.2f\t%.2f\t%.2f\t%.2f'%(tspn_res['WT'], tspn_res['E'], tspn_res['MBS'], tspn_res['MBR']))
	
	return AORP_res, AORP_xt, rtable_res, rtable_xt, tspn_res, tspn_xt
	
def plotLast10Min(xt):
	xt = np.array(xt)
	n = xt.shape[1]
	fig = plt.figure(figsize=[20,10])
	for i in range(2,n):
		plt.plot(xt[:,0]/60, xt[:,i], label = 'Queue %d'%(i-1))

	plt.xlim((xt[-1,0]/60)-10,xt[-1,0]/60)
	plt.ylim(0,400)
	plt.legend()
	plt.xlabel('Minutes of Operation')
	plt.ylabel('Queue Length (MB)')
	

def save_plt(name):
	#ax = plt.gca()
	#box = ax.get_position()
	#ax.set_position([0, 0, box.width*1.1, box.height*1.1])
	plt.savefig(name, format='png', bbox_inches='tight')
#######################################################
#######################################################
# Helper functions
#######################################################
#######################################################
def runsimsforpolicy(ps, rp, S, motion_power, tx_power, seconds, n_trials = 20):
	FAWT = 0
	TMBS = 0 #total mb serviced
	TMBR = 0 #total mb remaining
	TE = 0
	for i in range(n_trials):
		xt, wt, queues, total_travel_time = ps.simulate(rp, S, seconds)
		FAWT += wt[-1][1]
		MB_serviced = sum([len(q.wait_times) for q in queues])
		TMBS += MB_serviced
		TMBR += sum([len(q.waiting) for q in queues])
		TE += motion_power*total_travel_time + tx_power*ps.beta*MB_serviced

	return {'WT': FAWT/n_trials, 'MBS': TMBS/n_trials, 'MBR': TMBR/n_trials, 'E': TE/n_trials}, xt
	
def prob_pred_in_true(tjcf, idxs):
	n_in = 0
	if not len(idxs):
		return 0
		
	for pt in idxs:
		if tjcf[pt[0], pt[1]]:
			n_in += 1

	return n_in/len(idxs)

	
def field_to_pts(field, region, res):
	idcs = np.array(np.where(field>0)).T 
	return toRawFromGrid(region, res, idcs)

def _plot_relay_points(X, ax = plt):
	ax.plot(X[:,0], X[:,1], '*', markersize=25, markerfacecolor='palegreen', markeredgecolor='k', zorder = 101)
	ax.plot([1000], [100], '*', markersize=35, markerfacecolor='palegreen', markeredgecolor='k', label='Relay Positions')
	# ax.xlabel('x (m)')
	# ax.ylabel('y (m)')

def _format_legend(ax = plt):
	return ax.legend(fontsize=fs, loc='upper center', bbox_to_anchor=[0.5,-.05], ncol=3)