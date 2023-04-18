clear all
close all

% 读取csv文件数据
num = csvread('SRR_CLS2.csv',1,0); % SRR_CLS2.csv和本脚本放在同一个文件夹内

freq = num(:,1);
S11re = num(:,3);
S11im = num(:,2);
S21re = num(:,5);
S21im = num(:,4);

% 组合S参数
r = S11re+1j*S11im; % 反射系数S11
t1 = S21re+1j*S21im; % 传输系数S21

% 常量
d = 1.6e-3;% 立方体边长，单位m
m=0;
c=3e8;

% 利用S11和S21求解 折射率、介电常数和磁导率
k0=2*pi*freq*1e9/c; % 真空中波数
k0d=(k0*d);
r=conj(r);t1=conj(t1); % 两种傅氏变换之间转换 exp(j*omega*t) 和 exp(-j*omega*t)
z=sqrt(((1+r).^2-t1.^2)./((1-r).^2-t1.^2)); % 计算相对波阻抗
check_sign=-1*(real(z)<0)+(real(z)>=0);
z = check_sign.*z; % 确保z的实部为正
temp = acos((1-(r.^2-t1.^2))/2./t1)./(k0d);
n = temp+2*pi*m./(k0d); % 折射率n
x_axis = freq *1e9*d/c; % 频率归一化

% 绘制图形
subplot(2,2,1);
% plot(x_axis,real(z),'-b','linewidth',2);
plot(freq,real(z),'b-','linewidth',2);
hold on;
% plot(x_axis,imag(z),'-.r','linewidth',2);
plot(freq,imag(z),'r-.','linewidth',2);
title('Impedance z');
legend('Re','Im','Location','NorthEast');
legend('boxoff','fontsize',10);
%axis([0 13  -10 10]);

subplot(2,2,2);
% plot(x_axis,real(n),'-b','linewidth',2);
plot(freq,real(n),'b-','linewidth',2);
hold on;
% plot(x_axis,imag(n),'-.r','linewidth',2);
plot(freq,imag(n),'r-.','linewidth',2);
title('Refraction Index n');
legend('Re','Im','Location','NorthEast');
legend('boxoff','fontsize',10);
%axis([0 13 -10 10]);
epsilon=n./z;
mu=n.*z;

% =======================>以该段绘图代码为准<=========================
subplot(2,2,3);
plot(freq,real(epsilon),'k-',freq,imag(epsilon),'k--','Linewidth',2);
hold on;
plot(freq,zeros(length(freq),1),'b-','linewidth',2);
title('permittivity \epsilon');
xlabel('freq/GHz');
ylabel('');
legend('Re','Im','Location','NorthEast');
set(gca,'FontSize',15)
axis tight
box off
%axis([0 13 -10 10]);
% ==================================================================

subplot(2,2,4);
plot(freq,real(mu),'-b','linewidth',2);
hold on;
plot(freq,imag(mu),'-.r','linewidth',2);
plot(freq,zeros(length(freq),1),'k-','linewidth',1);
title('permeability \mu');
legend('Re','Im','Location','NorthEast');
legend('boxoff','fontsize',10);
%axis([0 13 -10 10]);
