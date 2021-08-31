"""
###################################################
% Predicted Channel
###################################################
Parameters for the probabilistic modeling of channel gain.
Characterizes the path loss, shadowing, and multipath componenets 
"""

import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import warnings
import random

#Import a few utility functions
import sys  
from pathlib import Path
path = Path(__file__)
top_dir = path.parent.parent
sys.path.insert(0, str(top_dir.absolute())+"/utils")
sys.path.insert(0, str(top_dir.absolute())+"/geometry")
from CoordTransforms import toRawFromGrid
from CoordTransforms import toGridFromRaw
import ChannelSim as ChSim


class ChannelParams:

	def __init__(self, qBase, nPL, kPL, sigmaSH, decorrSH, decorrMP, lam, kRic, corrMP, psdAtFC, sigmaMP):
		self.qBase = np.array(qBase) #base station location in [x,y], generally assumed to be meters
		self.nPL = nPL #path loss exponent
		self.kPL = kPL #path loss constant
		self.sigmaSH = sigmaSH #shadowing standard deviation
		self.decorrSH = decorrSH #shadowing decorrelation distance. Referred to beta in some papers
		self.decorrMP = decorrMP #multipath decorrelation distance.
		self.lam = lam #wavelength
		self.kRic = kRic #  if modeling multipath as rician, kRic parameterizes the distribution
		self.corrMP = corrMP # if corrMP == True, treating multipath affects as correlated. Else treat as uncorrelated
		self.psdAtFC = psdAtFC 
		self.sigmaMP = sigmaMP #if modeling shadowing as a zero-mean, log normal distribution (uncorrelated)


class CommChannel:

	def __init__(self, cp, region, res, gamma_PL = None):
		self.cp = cp
		self.region = region
		self.res = res

		# the rectangular region specifying the environment
		x_max = self.region[0]
		x_min = self.region[1]
		y_max = self.region[2]
		y_min = self.region[3]

		if x_max <= x_min or y_max <= y_min:
			warnings.warn('Region\'s max values must be greater than its min values') 

		self.gx, self.gy = np.meshgrid(np.linspace(x_min,x_max, int( np.rint( (x_max-x_min)*res) ) ), \
		np.linspace(y_min,y_max, int( np.rint((y_max-y_min)*res) )))

		self.yImax, self.xImax = self.gy.shape[0], self.gy.shape[1] 

		if gamma_PL is None:
			self.generatePL()
		else:
			self.gammaPLdB = gamma_PL


		#initialize shadowing and multipath to 0
		self.gammaSHdB = 0;
		self.gammaMPdB = 0;

		K_ric = self.cp.kRic;
		# TODO - implement makedist
		#self.ricDist = makedist('Rician', 's', sqrt(K_ric/(K_ric + 1)), 'sigma',1/sqrt(2*(1+K_ric)))


	def generatePL(self):
		"""
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% generatePL
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Generate the PL component of the channel based on the channel
		% params
		% Input:
		% self - reference to the CommChannel object
		%
		% Output:
		% gamma_PL_dB - matrix of PL powers (in dB) over the region
		"""
		x_b, y_b = self.cp.qBase
		gd = np.sqrt((self.gx - x_b)**2 + (self.gy - y_b)**2)
		gd = np.where(gd < 1/self.res, 1/self.res, gd)
		self.gammaPLdB = self.cp.kPL - 10*self.cp.nPL*np.log10(gd)
		self.gammaPLdB = self.gammaPLdB.T #corerct for meshgrids backward indexing

	def generateSH(self, nSin=5000):
		"""
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% generateSH
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Generate a realization of the channel's shadowing component over
		% the entire region
		% Input:
		% self - reference to the CommChannel object
		%
		% Output:
		% gamma_SH_dB - matrix of SH powers (in dB) over the region
		"""
		print("Generating shadowing...")
		self.gammaSHdB = ChSim.generate_shadowing(self.cp.sigmaSH**2, self.cp.decorrSH,
				        nSin, self.cp.psdAtFC, self.gx, self.gy, self.res)
		self.gammaSHdB = self.gammaSHdB.T
		print("Shadowing generation complete.")

	def generateMP(self, mode=1):
		""" 
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% generateMP
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Generate a realization of the channel's multipath fading 
		% component over the entire region
		% Input:
		% self - reference to the CommChannel object
		% mp_model - flag indicating which multipath model to use. 1
		%               indicates Rician (default). 2 indicates log normal.
		%               If 2, MP is assumed to be uncorrelated.
		%
		% Output:
		% gamma_MP_dB - matrix of MP gain (in dB) over the region
		"""
		print("Generating MP...")
		if mode == 1: #%Rician MP
			gamma_MP_lin = ChSim.generate_multipath(self.cp.lam, self.gx, self.gy,
						        self.res, self.cp.kRic, self.cp.corrMP)
			self.gammaMPdB = 10*np.log10(gamma_MP_lin)
			self.gammaMPdB = self.gammaMPdB.T
		elif mode == 2:
			#log normal model - MP are assumed to be uncorrelated
			sz = self.gammaPLdB.shape
			self.gammaMPdB = self.cp.sigmaMP*np.random.standard_normal(sz)
		else:
			warnings.warn("Invlalid mode passed to CommChannel.generateMP")
		print("MP generation complete.")

	def getGammaTOTdB(self):
		return self.gammaPLdB + self.gammaSHdB + self.gammaMPdB

	def getGammaTOTdBAtPoint(self, pt):
		grid_pt = toGridFromRaw(self.region, self.res, pt)
		return self._getGammaTOTdBAtGridPoint(grid_pt)

	def sampleChannel(self, n, sres = None):
		"""
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% sampleChannel
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Sample the channel
		% Input:
		% self - reference to the CommChannel object
		% n - number of samples
		%
		% Output:
		% sample_pos = n X 2 matrix with the x,y location of the samples in
		%               continuous coordinates
		% sample_vals = Channel gain in dB at the points sampled
		% sres = resolution for sampling grid. Useful to set when the channel's resolution is very fine
		"""
		if sres is None:
			sres = self.res
		nx = ((self.region[0] - self.region[1])*sres)//1
		ny = ((self.region[2] - self.region[3])*sres)//1
		#avoid duplicates by sampling w/o replacement
		samples = random.sample(range(nx*ny), n)
		grid_pts = np.array([list(divmod(sample, ny)) for sample in samples])
		sample_pos = toRawFromGrid(self.region, sres, grid_pts)
		sample_vals = self.getGammaTOTdBAtPoint(sample_pos)
		return sample_pos, sample_vals

	"""
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% getReqTXPowerW
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Calculates the required TX power at a point in the workspace
	% for a given BER, spectral efficiency, and noise power.
	%
	% Input:
	% self - reference to the CommChannel object
	% pt - point in question
	% qos - QoSParams object containiong BER, spectral efficiency, and
	%       noise power
	%
	% Output:
	% required transmit power at the given point for the given QoS requirement
	"""
	def getReqTXPowerAtPoint(self, pt, qos):
		#TODO - handle points that are outside of the region
		channel_gain = self.getGammaTOTdBAtPoint(pt)
		return qos.reqTXPower(channel_gain)

	def getConnectionField(self, gamma_th):
		"""
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% getConnectionField
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Computes binary field of connectivity for a given channel power
		% threshold.
		%
		% Input:
		% self - reference to the CommChannel object
		% gamma_th - minimum channel power required for communication (dB)
		%
		% Output:
		% connection_field - matrix with entires equal to 1 if reliable
		%                       communication possible at that point, 0
		%                       otherwise.
		"""
		gamma = self.getGammaTOTdB()
		return (gamma >= gamma_th)

	def getConnectedPoints(self, gamma_th):
		r = self.getConnectionField(gamma_th)
		idcs = np.array(np.where(r>0)).T
		return toRawFromGrid(self.region, self.res, idcs)

	"""
	Plotting Functions
	"""
	def plot_channel(self, chnlCmpFlg = 7, twod=False):
		fnt_siz = 20

		vals, title = self._getChannelVals(chnlCmpFlg)
		if twod:
			fig, ax = plt.subplots(figsize=(15,10))
			region = self.region
			pos = ax.imshow(vals.T, origin = 'lower', vmax=self.cp.kPL, 
				extent=[region[1],region[0],region[3],region[2]])
			fig.colorbar(pos, ax=ax)
		else:
			fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(15,10))
			ax.plot_surface(self.gx, self.gy, vals.T, linewidth=0)#flip back to meshgrid ordering
			ax.set_zlabel(title, fontsize= fnt_siz,  fontweight = 'bold')
		plt.xlabel('x (m)', fontsize= fnt_siz,  fontweight = 'bold')
		plt.ylabel('y (m)', fontsize= fnt_siz,  fontweight = 'bold')
			

	def plotConnectivityField(self, gamma_th):
		conn_field = self.getConnectionField(gamma_th)
		title = 'Connected Regions for Gamma_th = %.2f dB'% (gamma_th)
		plotConnectionField(self.region, conn_field, title)

	"""
	Private Helper Functions
	"""
	def _getGammaTOTdBAtGridPoint(self, grid_pts):
		gamma_TOT_dB = self.getGammaTOTdB()
		xis = tuple(grid_pts[:,0].tolist())
		yis = tuple(grid_pts[:,1].tolist())
		return gamma_TOT_dB[(xis, yis)]

	def _getChannelVals(self, chnlCmpFlg):
		#3-bit binary number
		# 1st bit - include PL
		# 2nd bit - include SH
		# 3rd bit - include MP
		vals = 0
		title = 'Channel Gain ('
		if (chnlCmpFlg % 2) == 1: #first bit is on
			vals += self.gammaPLdB
			title += 'PL +'
		if ( (chnlCmpFlg//2)%2) == 1: #second bit is on
			vals += self.gammaSHdB
			title += ' SH +'
		if ( (chnlCmpFlg//4) %2) == 1: #second bit is on
			vals += self.gammaMPdB
			title += ' MP '
		title =title[:-1]
		title += ") (dB)"
		return vals, title

class PredictedChannel(object):

	def __init__(self, cp, region, res, obs_pos, obs_vals, use_estimates = False):
		self.cp = cp

		self.region = region
		self.res = res

		# the rectangular region specifying the environment
		x_max = self.region[0]
		x_min = self.region[1]
		y_max = self.region[2]
		y_min = self.region[3]

		if x_max <= x_min or y_max <= y_min:
			warnings.warn('Region\'s max values must be greater than its min values') 

		self.gx, self.gy = np.meshgrid(np.linspace(x_min,x_max, int( np.rint( (x_max-x_min)*res) ) ), \
		np.linspace(y_min,y_max, int( np.rint((y_max-y_min)*res) )))


		self.obsPos = obs_pos
		self.obsVals = obs_vals
		self.useEstimates = use_estimates
		self.nObs = len(obs_vals)
		#TODO - include code for parameter estiamtion based on observations
		#       Currently will just use the true values contained in cp
		self._obsPL = self.estimatePLdBAtPoint(self.obsPos)
		self._obsDif = self.obsVals - self._obsPL
		self._obsCov = self._calcObsCov()
		self.UInv = np.linalg.inv(self._obsCov + self.rho()*np.eye(self.nObs))
		self._p_th = 0.85
		self._setStats()

	@property
	def p_th(self):
		return self._p_th

	def setPth(self, p):
		assert 0<=p<=1, 'Probability threshold must be in [0,1]'
		self._p_th = p

	def estimatePLdBAtPoint(self, pts):
		"""
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% estimatePLAtPoint
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Estimate the path loss gain at a point
		% Input:
		% self - reference to the PredictedChannel object
		% pts - either a single point [x,y] or matrix of points 
		%       [x1,y1;...xn,yn]. Assumes points are in raw coordinates
		% Output:
		% Returns the path loss gain in dB
		"""
		bs = self.cp.qBase
		pts = np.array(pts)

		return self.kPL() - 10*self.nPL()*np.log10(np.sqrt(np.sum((bs - pts)**2, pts.ndim-1)))
		
	def getConnectedPoints(self, gamma_th):
		r = self.getConnectionField(gamma_th)
		idcs = np.array(np.where(r)).T
		return toRawFromGrid(self.region, self.res, idcs)

	def getConnectionField(self, gamma_th):
		"""
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% getConnectionField
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Fetches a matrix representing a binary field with 1 indicating
		% communication is possible and 0 indicating connection is not
		% possible. This is based on a minimum required channel gain, which
		% in turn can be calculated from other requirements (e.g. SNR, BER,
		% etc...)
		% Input:
		% self - reference to the PredictedChannel object
		% gamma_th - the minimum required channel gain for communication
		%
		% Output:
		% conn_field - matrix of the same dimension as the grid-scaled
		%               workspace. Represents the binary connection field.
		"""
		x_vals = self.gx[0,:]
		y_vals = self.gy[:,0]
		grd_sz = self.gx.shape
		conn_field = np.zeros(np.flip(grd_sz))
		for i in range(grd_sz[1]):
			for j in range(grd_sz[0]):
				conn_field[i,j] = self.posteriorPConn([x_vals[i], y_vals[j]], gamma_th) >= self.p_th
		return conn_field

	def posteriorPConn(self, pt, gamma_th):
		"""
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% posteriorPConn
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Calculates the probability that the channel gain is above a giveb
		% channel gain threshold in dB. Prediction based on Gaussian process
		% Input:
		% self - reference to the PredictedChannel object
		% pt - a single point [x,y]
		% gamma_th - the minimum required channel gain for communication
		%
		% Output:
		% Returns the probability that the channel gain at pt is above
		%           gamma_th
		"""
		grid_pt = toGridFromRaw(self.region, self.res, pt)
		mean = self.means[grid_pt[:,0], grid_pt[:,1]]
		variance = self.vars[grid_pt[:,0], grid_pt[:,1]]

		return 1 - scipy.stats.norm.cdf(gamma_th, mean, np.sqrt(variance))
                
	def getPConField(self, gamma_th):
		grd_sz = self.gx.shape
		pField = np.zeros(np.flip(grd_sz))
		for i in range(grd_sz[1]):
			for j in range(grd_sz[0]):
				pField[i,j] = 1 - scipy.stats.norm.cdf(gamma_th, self.means[i,j], np.sqrt(self.vars[i,j]))

		return pField


	def kPL(self):
		if self.useEstimates:
			kPL = self.cp.kPL
			warnings.warn("Estimation not yet implemented")
		else:
			kPL = self.cp.kPL
		return kPL

	def nPL(self):
		if self.useEstimates:
			nPL = self.cp.nPL
			warnings.warn("Estimation not yet implemented")
		else:
			nPL = self.cp.nPL
		return nPL

	def rho(self):
		if self.useEstimates:
			rho = self.cp.sigmaMP**2
			warnings.warn("Estimation not yet implemented")
		else:
			rho = self.cp.sigmaMP**2
		return rho

	def beta(self):
		if self.useEstimates:
			beta = self.cp.decorrSH
			warnings.warn("Estimation not yet implemented")
		else:
			beta = self.cp.decorrSH
		return beta

	def alpha(self):
		if self.useEstimates:
			rho = self.cp.sigmaSH**2
			warnings.warn("Estimation not yet implemented")
		else:
			alpha = self.cp.sigmaSH**2
		return alpha

	"""
	Plotting Functions
	"""

	def plotConnectivityField(self, gamma_th):
		"""
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% plotConnected2D
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Plots the binary connectivity field over the workspace given a
		% minimum required channel gain and threshold probability.
		% Input:
		% self - reference to the PredictedChannel object
		% p_th - probability threshold. Since we're working with a
		%           predicted channel, this is the probability that the
		%           channel is good enough for communication
		% gamma_th - the minimum required channel gain for communication
		"""
		conn_field = self.getConnectionField(gamma_th)
		title = 'Connected Regions for Gamma_th = %.2f dB, p_th >= %.2f'% (gamma_th, self.p_th)
		plotConnectionField(self.region, conn_field, title)
	"""
	Private Functions
	"""

	def _calcObsCov(self):            
		diffs = np.zeros((self.nObs,self.nObs));
		diffs = np.array([ np.sqrt(np.sum((self.obsPos - pt)**2,1)) for pt in self.obsPos])
		cov = self.cp.sigmaSH**2*np.exp(-1*diffs/self.beta())
		return cov

	def _setStats(self):
		grd_sz = self.gx.shape
		x_vals = self.gx[0,:]
		y_vals = self.gy[:,0]
		self.means = np.zeros(np.flip(grd_sz))
		self.vars = np.zeros(np.flip(grd_sz))
		for i in range(grd_sz[1]):
			for j in range(grd_sz[0]):
				pt = [ x_vals[i], y_vals[j] ]
				self.means[i,j] = self._posteriorExpecteddB(pt)
				self.vars[i,j] = self._calcVariance(pt)		

	def _posteriorExpecteddB(self, pt):
		phi = self._varPosWithObs(pt)
		K = phi.T @ self.UInv
		mean = self._calcMean(K, pt)
		return mean

	def _calcMean(self, K, pt):
		return self.estimatePLdBAtPoint(pt) + K.T @ self._obsDif

	def _calcVariance(self, pt):
		k, phi = self._K(pt);
		var = self.alpha() + self.rho()  - k @ phi
		return var

	def _K(self, pt):
		phi = self._varPosWithObs(pt)
		k = phi.T @ self.UInv
		return k, phi

	def _varPosWithObs(self, pt):
		pt = np.array(pt)
		diffs = np.sqrt(np.sum((self.obsPos - pt)**2, 1))
		cov = self.cp.sigmaSH**2*np.exp(-1*diffs/self.beta());
		return cov


def plotConnectionField(region, field, title):
	plotField(region, field, title, is_binary = True)

def plotField(region, field, title, is_binary=False, cmap='cividis', tick_labels = None,
	ticks = None, vmin = None, vmax = None, do_show=True):
	field = field.astype(int)

	if is_binary:
		cmap = ListedColormap(['grey', 'green'])
		vmax, vmin = 1,0
	else:
		if vmax is None:
			vmax = np.max(field)
		if vmin is None:
			vmin = np.min(field)

	fig, ax = plt.subplots()
	fig.set_size_inches(12, 10)
	ws = ax.imshow(np.transpose(field), cmap=cmap, extent=[region[1],region[0],region[3],region[2]],
	vmax=vmax, vmin=vmin, interpolation='nearest', origin='lower', aspect='auto')
	cbar = fig.colorbar(ws)
	if is_binary:
		cbar.set_ticks([0.25, 0.75])
		cbar.set_ticklabels(["Disconnected", "Connected"])
	else:
		if ticks is not None:
			cbar.set_ticks(ticks)
		if tick_labels is not None:
			cbar.set_ticklabels(tick_labels)
	plt.title(title)

	if do_show:
		plt.show()