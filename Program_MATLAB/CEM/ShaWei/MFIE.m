%  MFIE.m
%  本程序用于磁场积分方程求解圆柱散射
%  编程人    沙威(Wei Sha) 安徽大学(Anhui University) ws108@ahu.edu.cn

clear;clc;
%%%% 1.常数

f=300*10^8;                 %  TM 波频率
v=3*10^8;                   %  波速
anglemin=pi/128;            %  最小离散角度步长
amplitude=1;                %  电场振幅

wavelength=v/f;             %  真空波长
k0=2*pi/wavelength;         %  真空波数
a=4*wavelength;             %  圆柱半径
count=2*pi/anglemin;        %  选配点数
circlelength=a*anglemin;    %  弧长

%%%%  2.坐标定位

%  圆柱中心建立坐标系
%  圆弧坐标

anglevalue=-anglemin/2;  %  角度初值

for m=1:count;
    anglevalue=anglevalue+anglemin;  %  角度递增
    point(m)=anglevalue;             %  角度赋值
end;



%%%%  3.矩阵赋值

for m=1:count;                     %  场循环
    for n=1:count;                 %  源循环
        
        field=point(m);            %  场点位置
        source=point(n);           %  源点位置
        
        if (m==n)
            L(m,n)=-1/2;           %  奇异点处理
        else
            quad_value=2*k0*a*abs( sin( (field-source)/2 ) );
            L(m,n)=(j*k0/4)*circlelength*(besselj(1,quad_value)-j*bessely(1,quad_value))*...
                (abs( sin( (field-source)/2 ) ));
        end;
    end;
end;

for t=1:count;  %  电压矩阵
    y(t)=-amplitude*exp(-j*a*cos(point(t))*k0);
end;



%%%%  4.求逆

result=L\y.';
figure(1);
plot(abs(result));
title('矩量法电流');



%%%%  5.rcs计算

t=0;
for s=0:pi/180:pi
    t=t+1;
    sum=0;
    for r=1:count;
        nr=cos(point(r)-s);  %  N.R,即远场与源点方向
        sum=sum+nr*result(r,1)*circlelength*...
            exp(j*k0*a*cos(point(r))*cos(s)+j*k0*a*sin(point(r))*sin(s));
    end;
    rcs1(t)=(k0/4)*abs(sum).^2;
end;

rcs=10.*log10(rcs1/wavelength);  %  波长归一化

figure(2);
plot(rcs);
title('双站雷达散射截面');



% %%%%  6.解析解电流
% 
% for m=1:count;
%     sum=0;
%     for n=-30:30;
%         t=(1/2)*(besselh(n+1,2,k0*a)-besselh(n-1,2,k0*a));  %  bessel函数求导
%         sum=sum+((-j*2*amplitude)/(pi*k0*a))*j^(-n)*exp(j*n*(point(m)))/t;
%     end;
%     current(m)=sum;
% end;
% 
% figure(3);
% title('电流对比');
% plot(abs(current),'r');
% title('电流分布对比');
% hold on;
% plot(abs(result),'bx');
% legend('解析','矩量')
% 
% 
% 
% %%%%  7.解析解RCS
% 
% RR=10.^6;  %  足够大的R
% t=0;
% for angle=0:pi/60:pi;
%     sum=0;
%     t=t+1;
%     for n=-30:30
%         cons2=(1/2)*(besselh(n+1,2,k0*a)-besselh(n-1,2,k0*a));
%         cons1=(1/2)*(besselj(n+1,k0*a)-besselj(n-1,k0*a));
%         A=cons1/cons2;
%         sum=sum+amplitude*j.^(-n)*A*besselh(n,2,k0*RR)*exp(j*n*angle);
%     end;
%     es(t)=sum;
% end;
% 
% RCS_ANA=2*pi*RR*abs(es).^2;
% RCS=10*log10(RCS_ANA/wavelength);
 
% figure(4);
% title('雷达散射截面对比');
% plot(RCS,'r');
% hold on;
% plot(rcs,'bx');
% legend('解析解','矩量解');
