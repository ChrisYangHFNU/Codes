import math
import numpy as np

zin = 100 # unit:Ohm
zout = 50
freqLow = 1e9
quarterwl = 0.25*187 # 波导波长 70.7Ohm, microstrip quarter wavelength, epsilon 3.38, frequency,1GHz
a = math.log(zout/zin)/quarterwl
x = np.linspace(0,quarterwl,8)
z = np.zeros([1,8])
exp(a*x)
print(zin*math.exp(a*x))

print(z)