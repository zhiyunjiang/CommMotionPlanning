"""
Channel Simulation Functions
"""

##########################################################################################
# Channel Simulator
# usage
#-------------------------
# [gamma_TOT_dB, gamma_PL_SH_dB, gamma_PL_dB, gamma_SH_dB, gamma_MP_LIN, g_x, g_y] 
# = channel_simulator(region, q_b, K_PL, n_PL, alpha, beta, N_sin, PSD_at_f_c, lambda, K_ric, res, corr_mp)
# 
# inputs
#-------------------------
# region     : a vector containing the boundary of the target rectangular region 
# q_b        : position of the base station (remote station or transmitter)
# K_PL, n_PL : path-loss parameters (gamma_PL_dB = K_PL - 10*n_PL * log10(d)) 
# alpha      : power of the shadowing (in dB)
# beta       : decorrelation distance (in meter)
# N_sin      : # of sinusoids to use 
# PSD_at_f_c : the amplitude difference (in dB) between PSD of shadowing at
#              cutoff frequency and at frequency 0
#              A more detailed description of this variable can be found in
#              the Documentation.
# lambda     : the wavelength of transmission (in meter)
# K_ric      : Rician K factor (K_ric = 0 results in Rayleigh distribution.)
# res        : the resolution of the grid (in samples/m) = 1/(dis between 2 sample along x or y axis)
# corr_mp    : corr_mp = 1 -> correlated multipath and corr_mp = 0 -> uncorrelated multipath
# 
# outputs
#-------------------------
# gamma_TOT_dB   : total channel (path loss + shadowing + multipath) (dB)
# gamma_PL_SH_dB : path loss + shadowing (dB)
# gamma_PL_dB    : path-loss only (dB) 
# gamma_SH_dB    : shadowing only (dB)
# gamma_MP_LIN   : multipath fading component in linear domain 
# g_x, g_y       : the corresponding 2D grid 
#
# Original codes courtesy of Alireza Ghaffarkhah and Alejandro Gonzalez-Ruiz
# Updated by Herbert Cai (April 2016)
#
# If you have any questions or comments regarding the codes,
# please contact Winston Hurst at winstonhurst@ece.ucsb.edu
#
# Dept. of Electrical and Computer Engineering
# University of California, Santa Barbara
###########################################################################
import numpy as np
import scipy.special

def channel_simulator(region, q_b, K_PL, n_PL, alpha, beta, N_sin, PSD_at_f_c, lam, K_ric, res, corr_mp):


	# the rectangular region specifying the environment
	x_max = region[0]
	x_min = region[1]
	y_max = region[2]
	y_min = region[3]

	# the position of the base station
	x_b = q_b[0]
	y_b = q_b[1]
	                         
	#######################################
	# generating the grid
	#######################################
	g_x, g_y = np.meshgrid(np.linspace(x_min,x_max, int( np.rint( (x_max-x_min)*res) ) ), \
		np.linspace(y_min,y_max, int( np.rint((y_max-y_min)*res) )))
	M, N = g_x.shape[0], g_x.shape[1]

	#######################################
	# path loss
	#######################################
	print('\ngrid size = %d pixels\n'%(M*N*res))
	print('generating path loss...\n')

	g_d = np.sqrt((g_x - x_b)**2 + (g_y - y_b)**2)

	# prevent the path loss to become very large if the samples are very close
	# to the base station
	g_d = np.where(g_d < 1/res, 1/res, g_d)

	# generating the path loss
	gamma_PL_dB = K_PL - 10*n_PL*np.log10(g_d)


	#######################################
	# shadowing
	#######################################
	print('generating shadowing...\n')
	gamma_SH_dB = generate_shadowing(alpha, beta, N_sin, PSD_at_f_c, g_x, g_y, res)


	#######################################
	# multipath fading
	#######################################
	print('generating multipath fading...\n')
	gamma_MP_LIN = generate_multipath(lam, g_x, g_y, res, K_ric, corr_mp)


	#######################################
	# overall channel
	#######################################

	gamma_PL_SH_dB = gamma_PL_dB + gamma_SH_dB
	gamma_PL_SH_LIN = np.power(10, (gamma_PL_SH_dB / 10) )
	gamma_TOT_LIN = gamma_PL_SH_LIN * gamma_MP_LIN
	gamma_TOT_dB = 10*np.log10(gamma_TOT_LIN)

	return gamma_TOT_dB, gamma_PL_SH_dB, gamma_PL_dB, gamma_SH_dB, gamma_MP_LIN, g_x, g_y




###########################################################################
# Shadow fading generator
#
# usage 
#-----------------------
# gamma_sh_dB = generate_shadowing(alpha, beta, N, PSD_at_f_c, g_x, g_y, res)
# 
# inputs 
#-----------------------
# alpha      : power of the shadowing (in dB)
# beta       : decorrelation distance (in meter)
# N          : # of sinusoids to use 
# PSD_at_f_c : the amplitude difference (in dB) between PSD of shadowing at
#              cutoff frequency and at frequency 0
#              A more detailed description of this variable can be found in
#              the Documentation.
# g_x, g_y   : 2D grid 
# res        : resolution in (samples / m)
#
# outputs
#-----------------------
# gamma_sh_dB : shadowing component (dB)
#
# References used in the comments =========================================
# This code follows the method described by:
# [1] X. Cai and G. B. Giannakis, “A two-dimensional channel
# simulation model for shadowing processes,IEEE Transactions
# on Vehicular Technology, vol. 52, no. 6, pp. 1558-1567, 2003.
#
# "JOR":(by Gonzales-Ruiz et al.)
# A Comprehensive Overview and Characterization of Wireless Channels for 
# Networked Robotic and Control Systems 
# This is our paper published in the Journal of Robotics, which has the
# mathematical details of the simulation (Section 6).
# =========================================================================
#
# Original codes courtesy of Alireza Ghaffarkhah and Alejandro Gonzalez-Ruiz
# Updated by Herbert Cai (April 2016)
# Translated to python by Winston Hurst (April 2021)
#
# If you have any questions or comments regarding the codes,
# please contact Winston Hurst at winstonhurst@ece.ucsb.edu
#
# Dept. of Electrical and Computer Engineering
# University of California, Santa Barbara
###########################################################################

def generate_shadowing(alpha, beta, N, PSD_at_f_c, g_x, g_y, res):


	#######################################
	# generating the shadowing component
	#######################################

	# standard deviation of shadowing
	sigma = np.sqrt(alpha);        

	a = 1/beta;

	# M is the number of radial frequencies, which is related to the number of sinusoids through: N = 2*M^2;
	M = int(np.round(np.sqrt(N/2)))
	                            
	# For mathematical details see [1] or the JOR paper (paper info. at the beginning of this file)

	PSD_at_f_c_lin = 10**(-PSD_at_f_c/10);                                  #see definition of epsilon on page 16 of JOR
	f_c = np.sqrt( 1/(4*np.pi**2)*((a**3/(PSD_at_f_c_lin))**(2/3)-a**2) )   #f_c is the cuttoff frequency, see definition of f_{r,c} on pager 16 of JOR
	P = 1 - a/(np.sqrt(a**2+4*np.pi**2*f_c**2))                              #P is the mean power within the cutoff frequency, we calculated this expression from Eq. 15 in [1]
	cn = np.sqrt(2/N)                                                       #coefficient (c_j) in Eq. 28 of JOR, the value coms from page 4 of [1]


	frm = np.zeros(M+1)  

	# compute the frequency frm recursively (Eq. 17 of [1])
	for m in range(1, M+1):
	    frm[m] = 1/(2*np.pi)*np.sqrt((P/(M*a)-1/np.sqrt(a**2 + 4*np.pi**2*frm[m-1]**2))**(-2)-a**2);  


	delta_phi = 2*np.pi/(4*M);         #see explanation below Eq. 29 of JOR

	# 2M angles uniformly distributed in ( -pi/2 , pi/2 )
	phi = np.linspace(-np.pi/2, np.pi/2-delta_phi/2 , 2*M);               

	fxn = np.zeros( (M, 2*M) )        #sampling frequencies along x and y
	fyn = np.zeros( (M, 2*M) )

	for m in range(0,M):
	    fxn[m,:] = (frm[m] + frm[m+1]) * np.cos(phi) /2
	    fyn[m,:] = (frm[m] + frm[m+1]) * np.sin(phi) /2

	theta = np.random.random_sample(fxn.shape)*2*np.pi          #the phase should be uniform between 0 and 2*pi

	sn = np.zeros(g_x.shape)

	for i in range(0, fxn.shape[0]-1):
	    for j in range(0, fxn.shape[1]-1):
	        sn =  sn + cn* np.cos(2*np.pi*(fxn[i,j]*g_x + fyn[i,j] *g_y) + theta[i,j])    #this is generating the actual shadowing 
	                                                                                   #by summing the sinusoids with appropriate frequencies and phase (see Eq. 28 of JOR)

	gamma_sh_dB = sigma*(sn-np.mean(sn))/np.std(sn);            #giving the shadowing the appropriate power

	return gamma_sh_dB



###########################################################################
# multipath fading generator
# This is generated using the filtering approach to generate correlated
# inphase and quadrature components assuming uniform angle of arrival, 
# which results in Jake's Bessel function for correlation.
#
# usage 
#-----------------------
# gamma_mp_lin = generate_multipath(lambda, g_x, g_y, res, K_ric, corr_mp)
#
# inputs 
#-----------------------
# lambda    : the wavelength of transmission (in meter)
# g_x, g_y  : 2D grid 
# res       : resolution in (samples / m)
# K_ric     : Rician K factor (k = 0 results in Rayleigh distribution)
# corr_mp   : corr_mp = 1 -> correlated multipath and corr_mp = 0 -> uncorrelated multipath
#
# outputs
#-----------------------
# gamma_mp_lin : multipath fading component in linear domain
#
# References used in the comments =========================================
# 1. "Goldsmith":
# Wireless Communications by A. Goldsmith
# This is the textbook that we will often refer to.
# =========================================================================
#
# Original codes courtesy of Alireza Ghaffarkhah and Alejandro Gonzalez-Ruiz
# Updated by Herbert Cai (April 2016)
#
# If you have any questions or comments regarding the codes,
# please contact Herbert Cai at hcai@ece.ucsb.edu
#
# Dept. of Electrical and Computer Engineering
# University of California, Santa Barbara
###########################################################################

def generate_multipath(lam, g_x, g_y, res, K_ric, corr_mp):

	M, N = g_x.shape;

	# We generate this because for multipath fading we need relative distance,
	# not distance between Tx and Rx.
	I, J = np.meshgrid(np.linspace(1,N+1, N), np.linspace(1,M+1, M))        

	I = I - np.round((N+1)/2);
	J = J - np.round((M+1)/2);

	d = np.sqrt(I**2 + J**2)/res;

	# This is the impulse response of the filter that is applied to the uncorrelated 
	# Gaussians (used in corr_gaussian_2D to generate correlated Gaussian).
	# generate a Bessel function (as shown in Fig. 3.5 of Goldsmith) in 2D
	h = scipy.special.jv(0, 2*np.pi*d/lam);          
	# inphase and quadrature components
	if corr_mp:
	    # correlated 
	    a = corr_gaussian_2D(0, 1, M, N, h);
	    b = corr_gaussian_2D(0, 1, M, N, h);
	else:
	    # uncorrelated 
	    a = np.random.normal(0, 1, (M, N))
	    b = np.random.normal(0, 1, (M, N))

	# Generating the correlated Rician
	# The standard way to generate a Rician is
	# to have R=sqrt(X^2 + Y^2) where X~(v,s^2) and Y~(0,s^2), where v and s
	# can be found as a function of K_ric in page 79 of Goldsmith.
	# Note that the average of resulting multipath fading power is one, as it should be.
	v = np.sqrt(K_ric/(K_ric + 1))
	s = 1/np.sqrt(2*(1+K_ric));

	a = a*s + v;
	b = b*s;

	gamma_mp_lin = a**2 + b**2
	#gamma_mp_dB = 10*log10(gamma_mp_lin);
	return gamma_mp_lin



###########################################################################
# Channel Simulator
# This function generates correlated Gaussian samples
#
# usage 
#-----------------------
# g = corr_gaussian_2D(mu, sigma, M, N, h)
#
# inputs
#-------------------------
# mu, sigma : the mean and std (standard deviation) of the Gaussian 
# M, N      : the desired size of the output matrix
# h         : the impulse response of the filter (used to generate gaussian with arbitrary correlation)
#
# outputs
#-------------------------
# g : the resulting Gaussian matrix with desired correlation
#
# References used in the comments =========================================
# 1. "JOR":(by Gonzales-Ruiz et al.)
# A Comprehensive Overview and Characterization of Wireless Channels for Networked Robotic and Control Systems 
# This is our paper published in the Journal of Robotics, which has the
# mathematical details of the simulation (Section 6).
# =========================================================================
#
# Original codes courtesy of Alireza Ghaffarkhah and Alejandro Gonzalez-Ruiz
# Updated by Herbert Cai (April 2016)
# Tranlsated to python by Winston Hurst (April 2021)
#
# If you have any questions or comments regarding the codes,
# please contact Winston Hurst at winstonhurst@ece.ucsb.edu
#
# Dept. of Electrical and Computer Engineering
# University of California, Santa Barbara
###########################################################################

def corr_gaussian_2D(mu, sigma, M, N, h):

	a = np.random.normal(0, 1, (M, N))            #generate an uncorrelated input

	# make it correlated by filtering, we multiply the samples in the frequency domain
	# with the appropriate filter and transform the result back to space domain
	# see Eq. 27 of JOR for details

	absH = np.sqrt(abs(np.fft.fft2(h)))
	A_f = absH * np.fft.fft2(a)                             
	a_f = np.fft.ifft2(A_f)
	g = a_f*sigma/np.std(a_f) + mu     #it gives the result the desired mean and std

	return g.real #in the processes of fft, ifft, a small complex component will appear. Ignore it.