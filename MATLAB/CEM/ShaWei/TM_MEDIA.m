%  TM_MEDIA.m

%  本程序用于等效原理和矩量法求解介质圆柱电流和磁流分布，及双站雷达散射截面。
%  本程序入射波为TM波,Z方向极化，X方向传播；介质柱轴向沿Z轴。
%  对于等效原理,采用内外等效相互结合和电场积分(EFIE)方程求解。
%  此程序没有考虑内谐振问题，若想消除内谐振，采用组合场积分(CFIE)方程或PMCHWT方程处理。
%  此程序采用“归一化”电流处理，避免了截断误差。但对内部问题必须注意。J外=J*Z0;J内=J外*Z1/Z0。
%  对于最小离散精度，每波长取10个点，但注意内外波长不同，要以最小波长为标准。
%  此程序当相对介电常数（epsilon）趋近无穷时，介质壁转化为金属壁。（可以验证）。
%  此程序当相对磁导率（u_r）趋近无穷时,介质壁转化为磁壁。若看成金属壁，相当于TE波的情况。(可以验证）。
%  此程序相对介电常数（epsilon）和相对磁导率（u_r）可以为复数，即此程序可以计算有耗介质圆柱的散射。
%  编程人沙威，安徽大学电磁场与微波专业
%  Email: dr.weisha@gmail.com
%  编程时间 2004-5-18

tic;
clc;clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.常数定义

f=300*10^6;                     %  TM 波频率
v=3*10^8;                       %  波速
S0=1/(4*pi*9*10^9);             %  真空介电常数
U0=4*pi*10^(-7);                %  真空磁导率
EULER=1.78107;                  %  欧拉常数
EX=exp(1);                      %  指数常数
anglemin=pi/80;                 %  最小离散角度步长
amplitude=1;                    %  电场振幅

wavelength=v/f;                 %  真空波长
k0=2*pi/wavelength;             %  真空波数
a=wavelength;                   %  介质柱半径
count=2*pi/anglemin;            %  选配点数
circlelength=a*anglemin;        %  弧长
Z0=sqrt(U0/S0);                 %  真空中波阻抗

epsilon=4;                      %  介质柱相对介电常数
u_r=1;                          %  介质柱相对磁导率
Z1=sqrt(u_r)*Z0/sqrt(epsilon);  %  介质柱波阻抗
k1=k0*sqrt(epsilon*u_r);        %  介质柱波数


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.坐标定位

%  圆柱中心建立坐标系


%  圆弧坐标
anglevalue=-anglemin/2;              %  角度初值

for m=1:count;
    anglevalue=anglevalue+anglemin;  %  角度递增
    point(m)=anglevalue;             %  角度赋值
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  3.阻抗矩阵赋值


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  A.场点在外圆弧

%(1)L0赋值,0表示0区域，即外等效。以下同样。

for m=1:count;                     %  场循环
    for n=1:count;                 %  源循环
        
        field0=point(m);           %  场点位置
        source0=point(n);          %  源点位置
        
        if (m==n)
            L0(m,n)=(k0/4)*circlelength*( 1-(j*2/pi)*log(EULER*k0*circlelength/(4*EX)) );
        else
            quad_value=2*k0*a*abs( sin( (field0-source0)/2 ) );
            L0(m,n)=(k0/4)*circlelength*(besselj(0,quad_value)-j*bessely(0,quad_value));
        end;
        
    end;
end;


%(2)M0赋值

for m=1:count;                     %  场循环
    for n=1:count;                 %  源循环

        field0=point(m);           %  场点位置
        source0=point(n);          %  源点位置
        
        if (m==n)                  
            M0(m,n)=1/2;           %  奇点处理
        else
            quad_value=2*k0*a*abs( sin( (field0-source0)/2 ) );
            M0(m,n)=-(j*k0/4)*circlelength*(besselj(1,quad_value)-j*bessely(1,quad_value))*...
                      abs( sin( (field0-source0)/2 ) );  %  角度N.R=|sin((si-si')/2)|
        end;
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  B.场点在内圆弧

%(1)L1赋值,1表示1区域，为内等效。以下同样。

for m=1:count;                     %  场循环
    for n=1:count;                 %  源循环
        
        field1=point(m);           %  场点位置
        source1=point(n);          %  源点位置
         
        if (m==n)
            L1(m,n)=(Z1/Z0)*(k1/4)*circlelength*( 1-(j*2/pi)*log(EULER*k1*circlelength/(4*EX)) );
        else
            quad_value=2*k1*a*abs( sin( (field1-source1)/2 ) );
            L1(m,n)=(Z1/Z0)*(k1/4)*circlelength*(besselj(0,quad_value)-j*bessely(0,quad_value));
        end;
        
    end;
end;


%(2)M1赋值

for m=1:count;                     %  场循环
    for n=1:count;                 %  源循环

        field1=point(m);           %  场点位置
        source1=point(n);          %  源点位置
        
        if (m==n)                 
            M1(m,n)=-1/2;          %  奇点处理
        else
            quad_value=2*k1*a*abs( sin( (field1-source1)/2 ) );
            M1(m,n)=-(j*k1/4)*circlelength*(besselj(1,quad_value)-j*bessely(1,quad_value))*...
                      abs( sin( (field1-source1)/2 ) );    %  角度N.R=|sin((si-si')/2)|
        end;
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  C.阻抗矩阵合并

IMPEDANCE_MATRIX=[ L0 ,M0 ;
                   L1 ,M1 ];                           
            
                            
                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  4.激励矩阵赋值

for t=1:count
    E(t)=amplitude*exp(-j*k0*a*cos(point(t)));  %  电场
end;


%%%%%矩阵合并
ZERO_MATRIX=zeros(1,count);
EXCITATION_MATRIX=[E,ZERO_MATRIX].';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  5.求解电流和磁流

RESULT=IMPEDANCE_MATRIX\EXCITATION_MATRIX;

figure(1);
plot(abs(RESULT(1:count)),'r');            %  电流
hold on;
plot(abs(RESULT(count+1:2*count)),'g');    %  磁流
legend('电流','磁流');
title('电磁流分布');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  6. RCS计算
 
t=0;  %  电流和磁流产生总场
for s=0:pi/60:pi
    
    t=t+1;
    sum=0;
    
    for r=1:count;
        coeff=cos(s-point(r));  %  磁流差乘时产生的角度
        sum=sum+(RESULT(r)-RESULT(r+count)*coeff)*circlelength*...
                exp(j*k0*a*cos(point(r))*cos(s)+j*k0*a*sin(point(r))*sin(s));  
                %  电磁流结合在一起算出来的双站雷达散射截面
    end;
    
    rcs(t)=((k0/4)*(abs(sum).^2))/wavelength; %  归一化 
end;

rcs=10.*log10(rcs);

figure(2);
plot(0:pi/60:pi,rcs,'b');
title('电磁流结合算出来的双站雷达散射截面');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  7.解析解  (哈林登正弦电磁场283页)

t=0;
for angle=0:pi/60:pi;   %  散射角
    
    sum=0;
    t=t+1;
    
    for n=-30:30;
        
        cons1=(1/2)*(besselj(n-1,k1*a)-besselj(n+1,k1*a));      %  贝塞尔函数求导
        cons2=(1/2)*(besselh(n-1,2,k0*a)-besselh(n+1,2,k0*a));  %  汉开尔函数求导 
        cons3=(1/2)*(besselj(n-1,k0*a)-besselj(n+1,k0*a));      %  贝塞尔函数求导
        
        d1=Z0*besselj(n,k0*a)*cons1-Z1*besselj(n,k1*a)*cons3;
        d2=Z1*besselj(n,k1*a)*cons2-Z0*besselh(n,2,k0*a)*cons1;
        A=d1/d2;  %  系数An

        sum=sum+A*exp(j*n*angle);  %  大自变量展开后的散射场
    end;
    
    es(t)=abs(sum).^2*4/k0;
end;

RCS=10*log10(es/wavelength);  %  归一化

figure(3);
plot(0:pi/60:pi,RCS,'R');
hold on;
plot(0:pi/60:pi,rcs,'bx');
legend('解析解','矩量解');
title('双站雷达散射截面对比');


