import math
import numpy as np
import matplotlib.pyplot as plt
import radarEq

pt = 1.5e6
freq = 5.6e9
g = 45
sigma = 0.1
te = 290  # effective noise temperature in Kelvins
b = 5.0e6
nf = 3.0
loss = 6.0
distance = np.linspace(25e3, 165e3, 1000)  # target range 25-165km,1000 points
snr1 = radarEq.eq_1p56(pt, freq, g, sigma, te, b, nf, loss, distance)
snr2 = radarEq.eq_1p56(pt, freq, g, sigma / 10, te, b, nf, loss, distance)
snr3 = radarEq.eq_1p56(pt, freq, g, sigma * 10, te, b, nf, loss, distance)

plt.figure()
distancekm = distance / 1000
p1, = plt.plot(distancekm, snr3, color='r', marker='+')
p2, = plt.plot(distancekm, snr1, color='k', marker='+')
p3, = plt.plot(distancekm, snr2, color='b', marker='+')
plt.legend([p1, p2, p3], ["\sigma=0dBsm", "\sigma=-10dBsm",
                          "\sigma=-20dBsm"], loc='upper right')
plt.grid()
plt.xlabel('Dectection range -Km')
plt.ylabel('SNR -dB')

plt.figure()
snr1 = radarEq.eq_1p56(pt, freq, g, sigma, te, b, nf, loss, distance)
snr2 = radarEq.eq_1p56(pt * 4, freq, g, sigma / 10, te, b, nf, loss, distance)
snr3 = radarEq.eq_1p56(pt * 1.8, freq, g, sigma *
                       10, te, b, nf, loss, distance)
p1, = plt.plot(distancekm, snr3, color='r', marker='+')
p2, = plt.plot(distancekm, snr1, color='k', marker='+')
p3, = plt.plot(distancekm, snr2, color='b', marker='+')
plt.legend([p1, p2, p3], ["Pt=2.16 MW", "Pt=1.5 MW",
                          "Pt=0.6MW"], loc='upper right')
plt.grid()
plt.xlabel('Dectection range -Km')
plt.ylabel('SNR -dB')

plt.show()
