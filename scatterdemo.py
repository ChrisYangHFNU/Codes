import matplotlib.pyplot as plt

height = [161,170,182,175,173,165]
weight = [50,58,80,70,69,55]

# 画图
plt.scatter(height,weight,alpha=0.7)
# 设置X轴标签
plt.xlabel('height')
# 设置Y轴标签
plt.ylabel('weight')
# 添加图的标题
plt.title('scatter demo')
# 图形显示
plt.show()
# 增加一个tag
