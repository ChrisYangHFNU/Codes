from csv import reader
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

'''
filename1 = 'origin.csv'
with open(filename1,'rt') as raw_data:
    readers = reader(raw_data,delimiter=',')
    x = list(readers)
    data1 = np.array(x).astype('float')
    print(data1.shape)

filename2 = 'sincurvefit.csv'
with open(filename2,'rt') as raw_data:
    readers = reader(raw_data,delimiter=',')
    x = list(readers)
    data2 = np.array(x).astype('float')
    print(data2.shape)

filename3 = 'sincurveslot.csv'
with open(filename3,'rt') as raw_data:
    readers = reader(raw_data,delimiter=',')
    x = list(readers)
    data3 = np.array(x).astype('float')
    print(data3.shape)    

# Create plots with pre-defined labels
ax = plt.subplot()
ax.plot(data1[10:50,0],data1[10:50,1],'k--',label='普通型')
ax.plot(data2[10:50,0],data2[10:50,1],'k:',label='正弦曲线')
ax.plot(data3[10:50,0],data3[10:50,1],'k',label='正弦曲线+缺陷地')

legend = ax.legend(loc='upper center',shadow=True,fontsize='x-large')
plt.grid(color='k',linestyle='-',linewidth=1)
plt.title('不同情况下的耦合器方向性')
plt.xlabel('Freq(GHz)')
plt.ylabel('dB(S(3,1))-dB(S(4,1))')
# Put a nicer background color on the legend
legend.get_frame().set_facecolor('C0')

plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.show()
'''

filename = 'couplerSL.csv'
with open(filename,'rt') as raw_data:
    readers = reader(raw_data,delimiter=',')
    x = list(readers)
    data = np.array(x).astype('float')
    print(data.shape)

ax = plt.subplot()
ax.plot(data[:,0],data[:,1],'k--',label='dB(S(1,1))')
ax.plot(data[:,0],data[:,2],'k:',label='dB(S(2,1))')
ax.plot(data[:,0],data[:,3],'k-^',label='dB(S(3,1))')
ax.plot(data[:,0],data[:,4],'k-*',label='dB(S(4,1))')
legend = ax.legend(loc='upper right',shadow=True,fontsize='x-large')
plt.title('原始2阶耦合器的性能')
plt.xlabel('Freq(GHz)')
plt.ylabel('dB')
plt.rcParams['font.sans-serif'] = ['Microsoft YaHei']
plt.show()