import numpy as np

# Mode
m = 1
n = 0

# Input parameters
f = 14e9
a = 0.01829
b = 0.00815

lightC = 3e8
lambdaWork = lightC/f
lambdaCutoff = 2/np.sqrt(np.power(m/a,2)+np.power(n/b,2))
print(lambdaCutoff)