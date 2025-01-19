import numpy as np
import matplotlib.pyplot as plt
import radarEq

ii = 0
val1 = np.zeros(200)
val2 = np.zeros(200)
val3 = np.zeros(200)
val4 = np.zeros(200)
for x in np.arange(0,20,0.1):
	val1[ii] = radarEq.incomplete_gamma(x,1)
	val2[ii] = radarEq.incomplete_gamma(x,3)
	val = radarEq.incomplete_gamma(x,6)
	val3[ii] = val
	val = radarEq.incomplete_gamma(x,10)
	val4[ii] = val
xx = np.arange(0,20,.1)

plt.figure()
p1, = plt.plot(xx,val1,color='k',linestyle='-')
p2, = plt.plot(xx,val2,color='k',linestyle=':')
p3, = plt.plot(xx,val3,color='k',linestyle='--')
p4, = plt.plot(xx,val4,color='k',linestyle='-.')
plt.legend([p1,p2,p3,p4],["N=1","N=3","N=6","N=10"])
plt.xlabel('x')
plt.ylabel('Incomplete Gamma function(x,N)')
plt.grid()
plt.show()