import matplotlib.pyplot as plt

x = [1,2,3,4]
y = [1.2,2.5,4.5,7.3]
x1 = [1,2,3,4,5]
y1 = [2.3,3.4,1.2,6.6,7.0]

# plot函数作图
plt.plot(x,y)
plt.scatter(x1,y1,color='r',marker='+')

# show函数展示该图，如果无该代码，则程序完成绘图但无显示
plt.show()