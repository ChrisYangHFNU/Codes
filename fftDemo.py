import numpy as np
from scipy.fftpack import fft,ifft
import matplotlib.pyplot as plt
from matplotlib.pylab import mpl

mpl.rcParams['font.sans-serif'] = ['SimHei']  # 显示中文
mpl.rcParams['axes.unicode_minus'] = False  # 显示负号

# 采样点选择1400个
x = np.linspace(0,1,1400)
y = 7*np.sin(2*np.pi*200*x)+5*np.sin(2*np.pi*400*x)+3*np.sin(2*np.pi*600*x)

fft_y = fft(y)

N = 1400
x = np.arange(N) #频率个数
half_x = x[range(int(N/2))] #取一半区间
abs_y = np.abs(fft_y) # 双边频谱
angle_y = np.angle(fft_y) #复数角度
normalization_y = abs_y/N #归一化处理
normalization_half_y = normalization_y[range(int(N/2))]

plt.subplot(231)
plt.plot(x,y)
plt.title('原始波形')

plt.subplot(232)
plt.plot(x,fft_y,'black')
plt.title('双边振幅谱（未求振幅绝对值）',fontsize=9,color='black')

plt.subplot(233)
plt.plot(x,abs_y,'r')
plt.title('双边振幅谱（未归一化）',fontsize=9,color='red')

plt.subplot(234)
plt.plot(x,angle_y,'violet')
plt.title('双边相位谱（未归一化）',fontsize=9,color='violet')

plt.subplot(235)
plt.plot(x,normalization_y,'g')
plt.title('双边幅度谱（归一化）',fontsize=9,color='green')

plt.subplot(236)
plt.plot(half_x,normalization_half_y,'blue')
plt.title('单边幅度谱（归一化）',fontsize=9,color='blue')

plt.show()

print(len(fft_y))
print(fft_y[0:5])
