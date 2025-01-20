%   BEM_SIMPLE.m
%   本程序用边界元方法求解正方形柱体内电位分布
%   编程人    沙威(Wei Sha) 安徽大学(Anhui University) ws108@ahu.edu.cn

clear;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  1.常数定义

a=6;          %  正方形长
N=3;          %  每边点数
minstep=a/N;  %  最小离散步长
TOTAL=N*4;    %  所有点数
C=1/2;        %  常数定义
NN=100;       %  积分离散精度 
V_L=300;      %  已知电压矩阵
xx=a/2;       %  方形内部任意一点X坐标
yy=a/2;       %  方形内部任意一点Y坐标

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  2.坐标定位
%  以方柱左下角为坐标原点建立坐标系
%  匹配点采用逆时针方向从(0,a/6)坐标点依次编号
%  方柱左右两边X为常数,方柱上下两边Y为常数

value_b=-minstep/2;   %  下侧初值
value_r=-minstep/2;   %  右侧初值
value_t=a+minstep/2;  %  上测初值
value_l=a+minstep/2;  %  左侧初值

for i=1:TOTAL;
    
    if (i>0 & i<N+1)          %  下侧
        value_b=value_b+minstep;
        point(1,i)=value_b;
        point(2,i)=0;
    
    elseif (i>N & i<2*N+1)    %  右侧
        value_r=value_r+minstep;
        point(1,i)=a;
        point(2,i)=value_r;
    
    elseif (i>2*N & i<3*N+1)  %  上侧
        value_t=value_t-minstep;
        point(1,i)=value_t;
        point(2,i)=a;

    else                      %  左侧
        value_l=value_l-minstep;
        point(1,i)=0;
        point(2,i)=value_l;
    end;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  3.H矩阵h_st确定

for s=1:TOTAL      %  场点循环
    for t=1:TOTAL  %  源点循环  
        
        if (s==t)  %  奇异点处理
            h_st(s,t)=C;
        else
            fieldpoint_x=point(1,s);   %  场点X坐标
            currentpoint_x=point(1,t); %  源点X坐标
            fieldpoint_y=point(2,s);   %  场点Y坐标
            currentpoint_y=point(2,t); %  源点Y坐标 
            
            current_x=linspace(currentpoint_x-minstep/2,currentpoint_x+minstep/2,NN); %  X积分变量离散
            current_y=linspace(currentpoint_y-minstep/2,currentpoint_y+minstep/2,NN); %  Y积分变量离散
            
            
            if (t>0 & t<N+1)|(t>2*N & t<3*N+1)  %  上下侧
                
                quad=abs(fieldpoint_y-currentpoint_y)./...
                    ((fieldpoint_x-current_x).^2+(currentpoint_y-fieldpoint_y).^2);
                h_st(s,t)=-(1/(2*pi))*trapz(current_x,quad);
                
                
            else    %  左右侧
                
                quad=abs(fieldpoint_x-currentpoint_x)./...
                    ((fieldpoint_x-currentpoint_x).^2+(current_y-fieldpoint_y).^2);
                h_st(s,t)=-(1/(2*pi))*trapz(current_y,quad);
                
            end;   
        end;
    end;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  4.K矩阵k_st确定

for s=1:TOTAL       %  场点循环
    for t=1:TOTAL   %  源点循环  
        
         if (s==t)  %  奇异点处理
             k_st(s,t)=-(log(minstep/2)-1)*minstep/(2*pi);
             
         else
            fieldpoint_x=point(1,s);   %  场点X坐标
            currentpoint_x=point(1,t); %  源点X坐标
            fieldpoint_y=point(2,s);   %  场点Y坐标
            currentpoint_y=point(2,t); %  源点Y坐标 
            
            current_x=linspace(currentpoint_x-minstep/2,currentpoint_x+minstep/2,NN); %  X积分变量离散
            current_y=linspace(currentpoint_y-minstep/2,currentpoint_y+minstep/2,NN); %  Y积分变量离散
            
            
            if ((t>0 & t<N+1)|(t>2*N & t<3*N+1))   %  上下侧
                 
                quad=log( ( (fieldpoint_x-current_x).^2 + ...
                            (currentpoint_y-fieldpoint_y).^2 ).^(1/2) );
                k_st(s,t)=-(1/(2*pi))*trapz(current_x,quad);
                
            else  %  左右侧
                
                quad=log( ( (fieldpoint_x-currentpoint_x).^2 + ...
                            (current_y-fieldpoint_y).^2 ).^(1/2) );
                k_st(s,t)=-(1/(2*pi))*trapz(current_y,quad);
                
            end;
        end;
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  5.矩阵整序
%   上下侧电荷分布已知
%   左右侧电压分布已知

H_K1=[h_st(:,[1:N]),-k_st(:,[N+1:2*N]),h_st(:,[2*N+1:3*N]),-k_st(:,[3*N+1:4*N])];
H_K2=[k_st(:,[1:N]),-h_st(:,[N+1:2*N]),k_st(:,[2*N+1:3*N]),-h_st(:,[3*N+1:4*N])];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  6.已知部分电压和电荷矩阵g_u确定

for u=1:TOTAL;
    
    if ( (u>3*N) & (u<4*N+1) )   %  上侧
        g_u(u)=V_L;    
    else
        g_u(u)=0;
    end;
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  7.求解剩下电荷电压分布并显示
charge_voltage=(H_K1)^(-1)*H_K2*g_u.';

disp('下侧电位从左到右:');
disp(charge_voltage(1:N));
disp('上侧电位从左到右:')
disp(charge_voltage(3*N:-1:2*N+1));
disp('右侧电荷从上到下:');
disp(charge_voltage(2*N:-1:N+1));
disp('左侧电荷从上到下:');
disp(charge_voltage(3*N+1:4*N));

voltage=[charge_voltage(1:N);g_u(N+1:2*N)';charge_voltage(2*N+1:3*N);g_u(3*N+1:4*N)'];
charge= [g_u(1:N)';charge_voltage(N+1:2*N);g_u(:,[2*N+1:3*N])';charge_voltage(3*N+1:4*N)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  8.方形内部测试点H1矩阵h_st1确定


for t=1:TOTAL  %  源点循环  
    
    fieldpoint_x=xx;           %  场点X坐标
    currentpoint_x=point(1,t); %  源点X坐标
    fieldpoint_y=yy;           %  场点Y坐标
    currentpoint_y=point(2,t); %  源点Y坐标 
    
    current_x=linspace(currentpoint_x-minstep/2,currentpoint_x+minstep/2,NN); %  X积分变量离散
    current_y=linspace(currentpoint_y-minstep/2,currentpoint_y+minstep/2,NN); %  Y积分变量离散
    
    
    if (t>0 & t<N+1)|(t>2*N & t<3*N+1) %  上下侧
        
        quad=abs(fieldpoint_y-currentpoint_y)./...
            ((fieldpoint_x-current_x).^2+(currentpoint_y-fieldpoint_y).^2);
        h_st1(t)=-(1/(2*pi))*trapz(current_x,quad);
        
    else %  左右侧
        
        quad=abs(fieldpoint_x-currentpoint_x)./...
            ((fieldpoint_x-currentpoint_x).^2+(current_y-fieldpoint_y).^2);
        h_st1(t)=-(1/(2*pi))*trapz(current_y,quad);
        
    end;   
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  9.方形内部测试点K1矩阵k_st1确定

for t=1:TOTAL  %  源点循环  
    
    fieldpoint_x=xx;           %  场点X坐标
    currentpoint_x=point(1,t); %  源点X坐标
    fieldpoint_y=yy;           %  场点Y坐标
    currentpoint_y=point(2,t); %  源点Y坐标 
    
    current_x=linspace(currentpoint_x-minstep/2,currentpoint_x+minstep/2,NN); %  X积分变量离散
    current_y=linspace(currentpoint_y-minstep/2,currentpoint_y+minstep/2,NN); %  Y积分变量离散
    
    
    if ((t>0 & t<N+1)|(t>2*N & t<3*N+1))   %  上下侧
        
        quad=log( ( (fieldpoint_x-current_x).^2 + ...
            (currentpoint_y-fieldpoint_y).^2 ).^(1/2) );
        k_st1(t)=-(1/(2*pi))*trapz(current_x,quad);
        
    else  %  左右侧
        
        quad=log( ( (fieldpoint_x-currentpoint_x).^2 + ...
            (current_y-fieldpoint_y).^2 ).^(1/2) );
        k_st1(t)=-(1/(2*pi))*trapz(current_y,quad);
        
    end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  10.求解内部测试点电位与解析解
resolve=k_st1*charge-h_st1*voltage;  %  代入离散化泊松公式
show=[xx,yy];

disp('在方柱内部电位值');
disp('    x=    y=');
disp(show);
disp('BEM方法为:');
disp(resolve);

analysis=V_L*(a-xx)/a;
disp('解析解为:')
disp(analysis)




