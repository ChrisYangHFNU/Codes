import matplotlib.pyplot as plt
import numpy as np
import csv

#theta = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350, 360])
#r = np.array([7.0844, 6.8243, 5.9143, 4.4119, 2.2474, -0.8388, -4.8587, -6.8045, -4.7941, -4.0153, -6.6739,  -15.3826,   -8.1930,   -4.0148,   -4.2558,   -8.6621,  -14.5864,   -8.2264,   -6.4739,   -9.2014, -15.0949,   -7.8333,   -4.2802,   -4.6870,   -9.3699,  -11.1716,   -5.5722,   -4.0728,   -5.1163,   -6.0394, -3.8577,   -0.6952,    1.9525,    4.0240,    5.6048,    6.6596,    7.0844])
with open('D:/01.yangjing/Codes/GitHub/Microwave_Antenna_RF/12.csv','r') as csvfile:
    reader = csv.reader(csvfile)
    column0 = [row[0] for row in reader] # 读取csv文件里的某一列
    print(column0)
    #for line in reader:
    #    print(line)
with open('D:/01.yangjing/Codes/GitHub/Microwave_Antenna_RF/12.csv','r') as csvfile:
    reader = csv.reader(csvfile)
    column1 = [row[1] for row in reader]
    print(column1)
with open('D:/01.yangjing/Codes/GitHub/Microwave_Antenna_RF/12.csv','r') as csvfile:
    reader = csv.reader(csvfile)
    column2 = [row[2] for row in reader]
    print(column1)

#theta = list(map(int,column0)) 
theta = [int(i) for i in column0] # 将list中的string转换为int
rho0 = [float(i) for i in column1]
rho1 = [float(i) for i in column2]
theta = np.array(theta) # 将list转换为array
theta = theta*np.pi/180
rho0 = np.array(rho0) # convert python list to numpy array class
rho1 = np.array(rho1)

plt.polar(theta,rho0,theta,rho1)
plt.ylim(rho0.min(),rho0.max())
plt.yticks([-15, -10, -5, 0, 5, 10])
plt.show()