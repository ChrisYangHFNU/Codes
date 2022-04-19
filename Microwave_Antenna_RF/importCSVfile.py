import matplotlib.pyplot as plt
import numpy as np
import csv

csvFile = open('D:/01.yangjing/Codes/GitHub/Microwave_Antenna_RF/12.csv','r')
reader = csv.reader(csvFile)

# 建立空字典
result = {}
for item in reader:
    # 忽略第一行
    if reader.line_num == 1:
        continue
    result[item[0]] = item[1]

csvFile.close()
print(result)