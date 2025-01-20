% 采用 Woodward-Lawson 抽样法综合余割平方波束方向图
% author:yangjing 2024/10/30
clear 

% constant parameters
eleNum = 32;
M = eleNum/2; % 取样点相关数
n = -M:1:M; % 采取奇数点抽样
xn = 2*n/eleNum;
theta_n = rad2deg(acos(xn)); % 抽样点角度
sampleDegree = -theta_n+90; % 将抽样角度范围转换为-90 ~ +90
an = [0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.8 1 0.94 0.86 0.78 0.69 0.59 0.01 0.01]; % 从目标函数中抽样得到相应的激励系数,即在对应角度的方向函数值 16个采样角度
an2 = [0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.8 1 1 0.97 0.94 0.9 0.86 0.82 0.78 0.74 0.69 0.65 0.59 0.54 0.01 0.01 0.01]; % 33 个角度抽样点

m = 1:eleNum;
zm = (m-(eleNum+1)/2)*0.5; % 确定阵列里的单元位置,单元间距二分之一波长
% 计算阵列各阵元的激励电流
Im = zeros(1,eleNum);
% for ii=1:length(m)
%     for jj=1:length(n)
%         tmp = an(jj)*exp(-1j*2*pi*zm(ii)*cosd(theta_n(jj)));
%         Im(ii) = Im(ii)+tmp;
%     end
% end
% --------- 利用矩阵运算实现上述二重循环 ----------
Im_2 = an2*exp(-1j*2*pi*zm.'*cosd(theta_n)).';
% -------------------------------------------

Im = Im/eleNum;
Im_2 = Im_2/eleNum;
amp = abs(Im_2);
pha = angle(Im_2);
PhaAmp = [pha amp];
thetaRange = -90:90;
fm = linearArrayfactor(eleNum,thetaRange,0.5,PhaAmp);