# -*- coding: utf-8 -*-
"""
Created on Sun Apr 13 19:36:49 2025

@author: chris
"""

import numpy as np
import matplotlib.pyplot as plt

# Parameters configuration
theta_min = np.deg2rad(5)
theta_max = np.deg2rad(30)
G0 = 1 # 基准增益

# Generate mwshgrid data
distance = np.linspace(1,50,500)
height = np.linspace(0,10,500)
D,H = np.meshgrid(distance,height)

# Calculate elevation angle
theta = np.arctan(H/D)

gain = np.zeros_like(D)

# Using csc second order model
valid_angles = (theta >= theta_min) & (theta <= theta_max)
gain[valid_angles] = G0/(np.sin(theta[valid_angles]) **2)

# generate graphics
plt.figure(figsize=(12,8))
contour = plt.contourf(D, H,gain,cmap='jet',extend='max')

levels = np.linspace(0, 20,50)

# add contour line
CS = plt.contour(D, H,gain,levels=[1,5,10,15,20],colors='black',linewidths=0.5)
plt.clabel(CS,inline=True,fontsize=10)

# label beam boundry
plt.plot(D[0], H[0]*np.tan(theta_min)/H[0],'w--',label=f'{np.rad2deg(theta_min):.0f}°仰角边界')
plt.plot(D[0], H[0]*np.tan(theta_max)/H[0],'w--',label=f'{np.rad2deg(theta_max):.0f}°仰角边界')

# graphic decorate
plt.colorbar(contour,label='Normalized Gain(dB)')
plt.xlabel('Horizontal Distance (km)',fontsize=12)
plt.ylabel('Altitude (km)',fontsize=12)
plt.title('Cosecant-Squared Beam Coverage Pattern (Blake Chart)',fontsize=14)
plt.legend()
plt.grid(True,alpha=0.3)
plt.xlim(0,50)
plt.ylim(0,50)
plt.show()