export = csvread('*.csv',1,0); % 将HFSS的结果*.csv文件读进MATLAB,要和本脚本文件放在同一文件夹下
% ====== 读入的数据分别放到不同的变量中 ======
freq = num(:,1);
S11re = num(:,3);
S11im = num(:,2);

% ======== 绘图 =========
plot(freq,S11re,'k-',freq,S11im,'k--','Linewidth',2); % 一幅图里绘制两条曲线
title('S11');
xlabel('freq/GHz');
ylabel('dB');
legend('Re','Im','Location','NorthEast');
set(gca,'FontSize',15)
axis tight
box off