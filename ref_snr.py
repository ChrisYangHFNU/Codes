import math
import radarEq

tau = 0.1e-6
R = 120e3
Sigma = 0.2
snr = radarEq.eq_1p60(tau,R,Sigma)
print('snr is:'snr)