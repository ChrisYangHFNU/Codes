function DemoPlot()
% 定义一个简单的绘图函数

x = 1:100;
y = log(x);
for i = 1:6
    plot(x, y + 3*i, 'LineWidth', 5);
    hold on
end
hold off
end