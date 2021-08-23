import numpy as np
import matplotlib.pyplot as plt

logpfa = np.linspace(.01,250,1000)
var = np.power(10,logpfa/10)
vtnorm = np.sqrt(np.log10(var))
plt.semilogx(logpfa,vtnorm)
plt.grid()
plt.show()