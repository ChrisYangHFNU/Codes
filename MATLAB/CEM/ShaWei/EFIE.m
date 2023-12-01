%  EFIE.m
%  本程序用于电场积分方程求解圆柱散射
%  编程人    沙威(Wei Sha) 安徽大学(Anhui University) ws108@ahu.edu.cn


tic
clc;clear;

%  常数定义
f=500*10^6;  %  TM 波频率
v=3*10^8;    %  波速
S0=1/(4*pi*9*10^9); %  真空介电常数
U0=4*pi*10^(-7);    %  磁导率
Z=sqrt(U0/S0);      %  真空中波阻抗
wavelength=v/f;     %  波长
k=2*pi/wavelength;  %  波数
EULER=1.78107;      %  欧拉常数
EX=exp(1);          %  指数常数


amplitude=1;         %  电场振幅
a=wavelength;        %  圆柱半径，为简了简便计算取 k*a=1
anglemin=pi/80;      %  最小离散角度步长
count=2*pi/anglemin; %  选配点数

circlelength=a*anglemin;  %  弧长

%  坐标定位(角度)
anglevalue=-anglemin/2;  %  角度初值
for i=1:count
    anglevalue=anglevalue+anglemin;  %  角度递增
    point(i)=anglevalue;  %  角度赋值
end



%  阻抗矩阵l_mn确定
for m=1:count  %  场点循环
    for n=1:count  %  源点循环  
        if (m==n)
            l_mn(m,n)=(1-(2*1i/pi)*log(EULER*k*circlelength/(4*EX)))*circlelength*(k*Z/4);
            %  奇异点处理
        else
            fieldpoint=point(m);   %  场点
            currentpoint=point(n); %  源点
   
            % %  被积函数
            integerfunction1=2*k*a*abs(sin((currentpoint-fieldpoint)/2));
            integerfunction2=besselj(0,integerfunction1)-1i*bessely(0,integerfunction1);
            % % 元素赋值
            l_mn(m,n)=(circlelength*k*Z/4)*integerfunction2; %  调用积分函数求解阻抗矩阵l_mn
        end  
    end  
end 
     
  

%  电压矩阵g_m确定
for ii=1:count
    g_m(ii)=amplitude*exp(-1i*k*a*cos(point(ii)) );  %  往X方向传播的TM波
end



%  求解电流分布
currentresult=l_mn\g_m.';  %  矩阵除法



%  解析解求解
initialcurrent=0;  %  求和初始化
for m=1:count
    for n=-20:20
        initialcurrent=initialcurrent+(2*amplitude/(2*pi*f*U0*pi*a))*...
                       1i^(-n)*exp(1i*n*point(m))/(besselj(n,k*a)-1i*bessely(n,k*a));  %  解析解电流
    end
        analysis(m)=initialcurrent;  %  赋值
        initialcurrent=0;  %  置零
end
    


%   RCS求解
ii=0;  %  矩阵下标
for si=0:pi/180:2*pi  %  散射角
    
    ii=ii+1;
    sum=0; %  求和初值
    
    for seg=1:count
       sum=sum+currentresult(seg,1)*circlelength*exp(1i*k*( a*cos(point(seg))*cos(si)+a*sin(point(seg))*sin(si) ) );
    end
    
    rcs(ii)=abs(sum).^2*k*Z^2/4;  %  RCS计算
end

rcs=10*log10(rcs*k/(2*pi));  %  归一化RCS



%  图形显示
figure(1)
plot(point*180/pi,abs(currentresult.'),'rx');  %  矩量法计算出的电流
hold on;
plot(point*180/pi,abs(analysis),'b');  %  解析法计算出的电流

% axis([0,360,min(abs(analysis))-10,max(abs(analysis))+10]); %  坐标轴调整
title('MOM APPLIED IN CYLINDER SCATTERING PROBLEM');  %  题目
xlabel('angle');  %  x坐标轴
ylabel('current');  %  y坐标轴
legend('mom result','analysis result');  %  标识
hold off
        
figure(2)
plot([0:1:180],rcs(1:181));
xlabel('scattered angle（/^{。}）');  %  y坐标轴
ylabel('Radar Cross Section /dB');  %  x坐标轴
title('RCS OF CYLINDER WITH SIDE-DIMENSION OF DOUBLE WAVELENGTHS');  %  题目


toc;

