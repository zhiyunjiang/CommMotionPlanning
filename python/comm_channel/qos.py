import numpy as np


class QoSReq:

	def __init__(self, ber, r, rx_noise):
		assert (0<= ber <=1), 'QosReq: BER must be in range [0,1]'
		assert r > 0, 'QoSReq: spectral efficiency r must be positive'
		assert rx_noise.mW >= 0, 'QoSReq: Receiver noise must be non-negative'
		
		self.ber = ber
		self.r = r
		self.rx_noise = rx_noise
		self.k = -1.5/np.log(5*self.ber)
		
		
	def reqTXPower(self, channel_gain_dB):
		CNR_lin = self._toLinCNR(channel_gain_dB)
		return self._calcPower(1/CNR_lin)
		
	def thresholdChannelGain(self, txpwr):
		CNR_lin = self._calcLinCNRThreshold(txpwr)
		return self._calcChannelPowerdBmTh(CNR_lin)
		
	def expReqTXPower(self, exp_ch_gain_dB, exp_ch_gain_var):
		exp_CNR_lin = self._toLinCNR(exp_ch_gain_dB)
		
		exp_CNR_lin_inverse = np.exp((np.log(10)/10)**2 * exp_ch_gain_var/2)/exp_CNR_lin
		
		return self._calcPower(exp_CNR_lin_inverse)
		
	def _toLinCNR(self, channel_power_dBm):
		return 10**(channel_power_dBm/10)/self.rx_Noise.mW
		
	def _calcPower(self, CNR_lin_inverse):
		power_mW = ((2**(self.r) -1)/self.k)*CNR_lin_inverse
		return TXPwr(Pwr.mW2dBm(power_mW))
		
	def _calcLinCNRThreshold(self, txpwr):
		CNR_inv = txpwr.W*1000/((2**self.r - 1)/self.k)
		return 1/CNR_inv
		
	def _calcChannelPowerdBmTh(self, CNR_lin_th):
		return 10*np.log10(self.rx_noise.mW*CNR_lin_th)
		
class Pwr(object):

	def __init__(self, dBm):
		self._dBm = dBm
		self._dBW = Pwr.dBm2dBW(dBm)
		self._mW = Pwr.dBm2mW(dBm)
		self._W = Pwr.dBm2W(dBm)
		
	def __str__(self):
		return '%f dBm\t%f dBW\t%f mW\t%f W'%(self._dBm, self._dBW, self._mW, self._W)
		
	@property
	def dBm(self):
		return self._dBm
	@property
	def dBW(self):
		return self._dBW
	@property
	def mW(self):
		return self._mW
	@property
	def W(self):
		return self._W
	
	#some handy static functions
	def dBm2dBW(dBm):
		return dBm - 30
		
	def dBm2mW(dBm):
		return 10**(dBm/10)
		
	def dBm2W(dBm):
		return Pwr.dBm2mW(dBm)/1000
		
	def mW2dBm(mW):
		return 10*np.log10(mW)
	
		
		
