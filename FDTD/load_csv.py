import os,sys
os.chdir(sys.path[0])
sys.path.append(os.getcwd())

import numpy as np
from matplotlib import pyplot as plt 

filename = 'SRR_CLS.csv'
with open(filename,'rt') as raw_data:
    data = np.loadtxt(raw_data,delimiter=',')
    print(data.shape)

print('current path:',os.getcwd())
print('sys path:',sys.path[0])

plt.rcParams['font.size'] = 12
plt.figure(figsize=(8, 3.5))
plt.plot(data, color='k', linewidth=1)
plt.show()
