import math
import numpy as np
import matplotlib as mpl

zl = 100
z0 = 50
N = 3
A = math.pow(2,-N)*(zl-z0)/(zl+z0)
Gamma_m = 0.32
bw = 2-4/math.pi*math.acos(0.5*math.pow(Gamma_m/np.abs(A),1/N))

print(bw)

Z_next = zeros(1,N)
for ii=1:N;
