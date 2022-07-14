import numpy as np
import math

Zs = 127 # slotline特性阻抗
d = 0.203 # substrate厚度 mm
epsilon_r = 3.55 # substrate relative epsilon
wavelength = 12.5 # 工作波长 mm
wavelength_sl = wavelength*math.sqrt(2/(epsilon_r+1)) # slotline wavelength mm
u = math.sqrt(epsilon_r-pow(wavelength/wavelength_sl,2))
v = math.sqrt(pow(wavelength/wavelength_sl,2)-1)
q = 2*math.pi*u*d/wavelength+math.atan(u/v)
N = math.cos(2*math.pi*u*d/wavelength)-1/math.tan(q)*math.sin(2*math.pi*u*d/wavelength)
Zm = Zs*pow(N,2) # 交叉点microstrip输入阻抗
print(Zm)