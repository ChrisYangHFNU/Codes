import numpy as np
import matplotlib.pyplot as plt
import radarEq

pfa1 = 1.0e-12
pfa2 = 1.0e-12
pfa3 = 1.0e-12
pfa4 = 1.0e-12
pd1 = .5
pd2 = .8
pd3 = .95
pd4 = .99
Index = 0
L1 = np.zeros([1000,1])
L2 = np.zeros([1000,1])
L3 = np.zeros([1000,1])
L4 = np.zeros([1000,1])
for npulse in (np.arange(1,1001,1)):
	I1 = radarEq.improv_fac(npulse,pfa1,pd1)
	i1 = np.power(10,0.1*I1)
	L1[Index] = -1*10*np.log10(i1/npulse)
	I2 = radarEq.improv_fac(npulse,pfa2,pd2)
	i2 = np.power(10,0.1*I2)
	L2[Index] = -1*10*np.log10(i2/npulse)
	I3 = radarEq.improv_fac(npulse,pfa3,pd3)
	i3 = np.power(10,0.1*I3)
	L3[Index] = -1*10*np.log10(i3/npulse)
	I4 = radarEq.improv_fac(npulse,pfa4,pd4)
	i4 = np.power(10,0.1*I4)
	L4[Index] = -1*10*np.log10(i4/npulse)
	Index = Index+1
npulse = np.arange(1,1001,1)

plt.figure()
p1, = plt.semilogx(npulse,L1,linestyle='-')
p2, = plt.semilogx(npulse,L2,linestyle='--')
p3, = plt.semilogx(npulse,L3,linestyle='-.')
p4, = plt.semilogx(npulse,L4,linestyle=':')
plt.xlabel('Number of pulses')
plt.ylabel('Integration loss - dB')
plt.legend([p1,p2,p3,p4],["pd=.5,nfa=e+12","pd=.8,nfa=e+12","pd=.95,nfa=e+12","pd=.99,nfa=e+12"],loc='lower right')
plt.grid(which='both')
plt.show()
	