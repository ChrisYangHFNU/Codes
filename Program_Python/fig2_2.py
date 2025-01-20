import numpy as np
from scipy.stats import norm,rayleigh
import matplotlib.pyplot as plt

xg = np.linspace(-6,6,1500)
xr = np.linspace(0,6,1500)
mu = 0
sigma = 1.5
ynorm = norm.pdf(xg,mu,sigma) #Gaussian pdf
yray = rayleigh.pdf(xr,sigma) #rayleigh pdf

plt.figure()
p1, = plt.plot(xg,ynorm,color='k',linestyle='-')
p2, = plt.plot(xr,yray,color='k',linestyle='-.')
plt.legend([p1,p2],["Gaussian pdf","Rayleigh pdf"],loc='upper right')
plt.xlabel('x')
plt.ylabel('Probability density')
plt.show()