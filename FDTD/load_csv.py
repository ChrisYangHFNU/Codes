import numpy as np
from matplotlib import pyplot as plt 

filename = 'SRR_CLS.csv'
with open(filename,'rt') as raw_data:
    data = np.loadtxt(raw_data,delimiter=',')
    print(data.shape)

plt.rcParams['font.size'] = 12
plt.figure(figsize=(8, 3.5))
plt.plot(data, color='k', linewidth=1)
plt.show()