import math
import numpy as np
from scipy.fftpack import fft,ifft

def eq_1p56(pt,freq,g,sigma,te,b,nf,loss,distance):
    c = 3.0e8
    wavelength = c/freq
    p_peak = 10*math.log10(pt)
    wavelength_sqdb = 10*math.log10(math.pow(wavelength,2))
    sigma_db = 10*math.log10(sigma)
    fourPIcub = 10*math.log10(math.pow(4.0*math.pi,3))
    k_db = 10*math.log10(1.38e-23)
    te_db = 10*math.log10(te)
    b_db = 10*math.log10(b)
    distance_pwr4_db = 10*np.log10(np.power(distance,4))
    
    num = p_peak+g+g+wavelength_sqdb+sigma_db
    den = fourPIcub+k_db+te_db+b_db+nf+loss+distance_pwr4_db
    snr = num-den
    return snr
    
def eq_1p60(tau,R,Sigma):
    Rref = 86e3
    tau_ref = 0.1e-6
    SNRref = 20
    snrref = math.pow(10,SNRref/10)
    Sigmaref = 0.1
    Lossp = 2
    lossp = math.pow(10,Lossp/10)
    
    rangeratio = np.power((Rref/R),4)
    snr = snrref*(tau/tau_ref)*(1/lossp)*(Sigma/Sigmaref)*rangeratio
    snr = np.log10(snr)
    return snr
    
def power_aperture(snr,tsc,sigma,distance,te,nf,loss,az_angle,el_angle):
    Tsc = 10*np.log10(tsc)
    Sigma = 10*np.log10(sigma)
    fourPI = 10*np.log10(4.0*math.pi)
    k_db = 10*math.log10(1.38e-23)
    Te = 10*math.log10(te)
    distance_pwr4_db = 10*np.log10(np.power(distance,4))
    omega = az_angle*el_angle/math.pow(57.296,2)
    Omega = 10*math.log10(omega)
    PAP = snr+fourPI+k_db+Te+nf+loss+distance_pwr4_db+Omega-Sigma-Tsc
    return PAP
    
def pulse_integration(pt,freq,g,sigma,te,b,nf,loss,distance,npulse,ci_nci):
    snr1 = eq_1p56(pt,freq,g,sigma,te,b,nf,loss,distance)
    if (ci_nci==1): #coherent integration
        snrout = snr1+10*np.log10(npulse)
    else:  #non-coherent integration
        if (ci_nci==2):
            snr_nci = np.power(10,snr1/10)
            val1 = np.power(snr_nci,2)/(4*npulse*npulse)
            val2 = snr_nci/npulse
            val3 = snr_nci/(2*npulse)
            SNR_1 = val3+np.sqrt(val1+val2)
            LNCI = (1+SNR_1)/SNR_1
            snrout = snr1+10*np.log10(npulse)-10*np.log10(LNCI)
    return snrout
    
def que_func(x):
# This function computes the value of the Q-function
# listed in Eq.(2.16), It uses the approximation in Eqs.(2.17) and (2.18)
    if x>=0:
    	denom = 0.661*x+0.339*np.sqrt(np.power(x,2)+5.51)
    	expo = np.exp(-np.power(x,2)/2.0)
    	fofx = 1.0-(1.0/np.sqrt(2.0*np.pi))*(1.0/denom)*expo
    else:
    	denom = 0.661*x+0.339*np.sqrt(np.power(x,2)+5.51)
    	expo = np.exp(-np.power(x,2)/2.0)
    	val = 1.0-(1.0/np.sqrt(2.0*np.pi))*(1.0/denom)*expo
    	fofx = 1.0-val
    return fofx
    
def marcumsq(a,b):
# This function uses Parl's method to compute PD
    max_test_value = 5000
    if (a<b):
      alphan0 = 1.0
      dn = a/b
    else:
      alphan0 = 0
      dn = b/a
    alphan_1 = 0
    betan0 = 0.5
    betan_1 = 0
    D1 = dn
    n = 0
    ratio = 2.0/(a*b)
    r1 = 0
    betan = 0
    alphan = 0
    while (betan<1000):
      n = n+1
      alphan = dn+ratio*n*alphan0+alphan
      betan = 1.0+ratio*n*betan0+betan
      alphan_1 = alphan0
      alphan0 = alphan
      betan_1 = betan0
      betan0 = betan
      dn = dn*D1
    PD = (alphan0/(2*betan0))*np.exp(-np.power((a-b),2)/2)
    if (a>=b):
    	PD = 1-PD
    return PD
    
def improv_fac(npulse,pfa,pd):
# this function computes the non-coherent integration improvement
# factor using the empirical formula defined in Eq.(2.49)
    fact1 = 1.0+np.log10(1.0/pfa)/46.6
    fact2 = 6.79*(1.0+0.235*pd)
    fact3 = 1.0-0.14*np.log10(npulse)+0.0183*np.power(np.log10(npulse),2)
    impr_of_npulse = fact1*fact2*fact3*np.log10(npulse)
    return impr_of_npulse
    
def incomplete_gamma(vt,npulse):
# this function implements Eq.(2.67) to compute the Incomplete Gamma
# this function needs "factor" function to run
    eps = 1.000000001
    if (npulse==1):
    	value1 = vt*np.exp(-vt)
    	value = 1.0-np.exp(-vt)
    	return value
    sumold = 1.0
    sumnew = 1.0
    calc1 = 1.0
    calc2 = npulse
    xx = npulse*np.log(vt+0.0000000001)-vt-factor(calc2)
    temp1 = np.exp(xx)
    temp2 = npulse/(vt+0.0000000001)
    diff = .0
    ratio = 1000.0
    if (vt>=npulse):
    	while (ratio>=eps):
    		diff = diff+1.0
    		calc1 = calc1*(calc2-diff)/vt
    		sumnew = sumold+calc1
    		ratio = sumnew/sumold
    		sumold = sumnew
    	value = 1.0-temp1*sumnew*temp2
    	return value
    else:
    	diff = 0.
    	sumold = 1.
    	ratio = 1000.
    	calc1 = 1.
    	while (ratio>=eps):
    		diff = diff+1.0
    		calc1 = calc1*vt/(calc2+diff)
    		sumnew = sumold+calc1
    		ratio = sumnew/sumold
    		sumold = sumnew
    	value = temp1*sumnew
    return value 
    
def factor(n):
	# compute the factorial of n using logarithms to avoid overflow
	n = n+9.0
	n2 = n*n
	temp = (n-1)*np.log(n)-n+np.log(np.sqrt(2.0*np.pi*n))+((1.0-(1.0/30.+(1.0/105)/n2)/n2)/12)/n
	val = temp-np.log((n-1)*(n-2)*(n-3)*(n-4)*(n-5)*(n-6)*(n-7)*(n-8))
	return val