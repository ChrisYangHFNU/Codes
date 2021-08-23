import numpy as np
import matplotlib.pyplot as plt
import radarEq

pfa1 = 1.0e-2
pfa2 = 1.0e-6
pfa3 = 1.0e-10
pfa4 = 1.0e-13
pd1 = .5
pd2 = .8
pd3 =.95
pd4 = .999
Index = 0
I1 = np.zeros([1000,1])
I2 = np.zeros([1000,1])
I3 = np.zeros([1000,1])
I4 = np.zeros([1000,1])
for npulse in (np.linspace(1,1000,1000)):
	I1[Index] = radarEq.improv_fac(npulse,pfa1,pd1)
	I2[Index] = radarEq.improv_fac(npulse,pfa2,pd2)
	I3[Index] = radarEq.improv_fac(npulse,pfa3,pd3)
	I4[Index] = radarEq.improv_fac(npulse,pfa4,pd4)
	Index = Index+1
npulse = np.linspace(1,1000,1000)

plt.figure()
p1, = plt.semilogx(npulse,I1,linestyle='-')
p2, = plt.semilogx(npulse,I2,linestyle='--')
p3, = plt.semilogx(npulse,I3,linestyle='-.')
p4, = plt.semilogx(npulse,I4,linestyle=':')
plt.xlabel('Number of Pulses')
plt.ylabel('Improvement factor I-dB')
plt.legend([p1,p2,p3,p4],["pd=.5,nfa=e+2","pd=.8,nfa=e+6","pd=.95,nfa=e+10","pd=.999,nfa=e+13"],loc='lower right')
plt.grid(which='both')
plt.show()