import math
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import radarEq

npulse = np.linspace(1,10000,1000)
snrci = radarEq.pulse_integration(4,94e9,47,20,290,20e6,7,10,5.01e3,npulse,1)
snrnci = radarEq.pulse_integration(4,94e9,47,20,290,20e6,7,10,5.01e3,npulse,2)

p1, = plt.semilogx(npulse,snrci,linestyle='-')
p2, = plt.semilogx(npulse,snrnci,linestyle='-.')
plt.legend([p1,p2],["coherent integration","Non_coherent integration"],loc='lower right')
plt.grid()
plt.xlabel('Number of integrated pulses')
plt.ylabel('SNR-dB')

plt.show()