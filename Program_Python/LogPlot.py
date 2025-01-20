import matplotlib.pyplot as plt
import numpy as np

# 设置随机数种子以便结果可复现
np.random.seed(0)

# 设置图形尺寸
plt.figure(figsize=(10, 6))

# 设置横坐标和纵坐标的范围
x = np.arange(0.15, 0.35, 0.025)
y = np.logspace(-4, -1, len(x))

# 定义颜色、线型和标记
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
linestyles = ['-', '--', '-.', ':']
markers = ['o', 's', '^', '*', 'D']

# 生成几条随机折线
for i in range(5):
     plt.plot(x, y * np.random.rand(len(x)), label=f'Line {i+1}', color=colors[i % len(colors)], linestyle=linestyles[i % len(linestyles)], marker=markers[i % len(markers)], linewidth=2)
     

# 设置纵坐标为对数坐标
plt.yscale('log')

# 添加网格
plt.grid(True, which="both", ls="-", color='0.7')

# 设置标题和坐标轴标签
plt.title('Random Lines in Logarithmic Scale')
plt.xlabel('X-axis')
plt.ylabel('Y-axis (log scale)')

# 设置图例
plt.legend()

# 显示图形
plt.show()