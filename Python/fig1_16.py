import math
import numpy as np
import matplotlib.pyplot as plt
import radarEq

tsc = 2.5
sigma = 0.1
te = 900
snr = 15
nf = 6
loss = 7 # radar losses in dB
az_angle = 2 # search volume azimuth extent in degrees
el_angle = 2;
distance = np.linspace(20e3,250e3,1000);
pap1 = radarEq.power_aperture(snr,tsc,sigma/10,distance,te,nf,loss,az_angle,el_angle)
pap2 = radarEq.power_aperture(snr,tsc,sigma,distance,te,nf,loss,az_angle,el_angle)
pap3 = radarEq.power_aperture(snr,tsc,sigma*10,distance,te,nf,loss,az_angle,el_angle)

plt.figure()
distancekm = distance/1000;
plt.plot(distancekm,pap1,color='b',marker='+')
p1, = plt.plot(distancekm,pap1,color='b',marker='+')
p2, = plt.plot(distancekm,pap2,color='r',marker='+')
p3, = plt.plot(distancekm,pap3,color='k',marker='+')
plt.legend([p1,p2,p3],["$\sigma=-20dBsm$","$\sigma=-10dBsm$","$\sigma=0dBsm$"],loc='upper right')
plt.grid()
plt.xlabel('Detection range in km')
plt.ylabel('Power aperture product in dB')

wavelength = 0.03
G = 45
ae = np.linspace(1,25,1000) #aperture size 1 to 25 meter squared,1000 points
Ae = 10*np.log10(ae)
distance = 250e3
pap1 = radarEq.power_aperture(snr,tsc,sigma/10,distance,te,nf,loss,az_angle,el_angle)
pap2 = radarEq.power_aperture(snr,tsc,sigma,distance,te,nf,loss,az_angle,el_angle)
pap3 = radarEq.power_aperture(snr,tsc,sigma*10,distance,te,nf,loss,az_angle,el_angle)
Pav1 = pap1-Ae
Pav2 = pap2-Ae
Pav3 = pap3-Ae
plt.figure()
p1, = plt.plot(ae,Pav1,color='k',linestyle='-')
p2, = plt.plot(ae,Pav2,color='k',linestyle='-.')
p3, = plt.plot(ae,Pav3,color='k',linestyle=':')
plt.legend([p1,p2,p3],["$\sigma$=-20dBsm","$\sigma$=-10dBsm","$\sigma$=0dBsm"],loc='upper right')
plt.xlabel('Aperture size in square meters')
plt.ylabel('Pav in dB')
plt.grid()

plt.show()
