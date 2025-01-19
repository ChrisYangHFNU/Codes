import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax
from matplotlib.pyplot import MultipleLocator
import radarEq
#import ipdb

sequence = np.array([2,4,6,8,10,12])
plt.figure()
for nfa in sequence:
	b = np.sqrt(-2.0*np.log(np.power(float(10),-nfa)))
	Index = 0
	pro = np.zeros([181,1])
	for snr in np.linspace(0,18,181):
		a = np.sqrt(2.0*np.power(10,0.1*snr))
		#pdb.set_trace() #代码调试点
		pro[Index][0] = radarEq.marcumsq(a,b)
		Index = Index+1
	x = np.linspace(0,18,181)
	plt.loglog(x,pro,color='k')
	my_x_ticks = np.arange(1,18,1)
	my_y_ticks = np.array([.1,.2,.3,.4,.5,.6,.7,.8,.85,.9,.95,.9999])
	#ax = plt.gca()
	#ax.xaxis.set_major_locator(x_major_locator) # 把x轴的主刻度设置为1的倍数
	#ax.yaxis.set_major_locator(y_major_locator)
	#ax.Axes.set_xticks(ticks='my_x_ticks')
	plt.xticks(my_x_ticks,rotation='45')
	plt.yticks(my_y_ticks) # 设置y轴主刻度为指定的序列
	plt.xlim(0.9,18)
	plt.ylim(0.1,1)
plt.xlabel('Single pulse SNR -dB')
plt.ylabel('Probability of Detection')
plt.grid(which='both')
plt.show()
	
