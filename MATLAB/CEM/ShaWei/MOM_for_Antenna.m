%  MOM_for_Antenna.m
%  本函数用于实现用矩量法和海伦积分方程球解对称天线电流分布
%  编程人    沙威(Wei Sha) 安徽大学(Anhui University) ws108@ahu.edu.cn


%  wave_length_value 波长
%  measurement 物体电尺寸
%  v0 电压常数
%  count_point 待求点数
%  range 求解的天线归一化长度
%  a 导体半径
%  k 波数
%  l 天线长度
%  step 点匹配间距

clear;clc;

range=1/4;
point=3;
wave_length_value=1;
measurement=7.022*(10)^(-3);
v0=1;

k=(2*pi)/wave_length_value;
a=measurement*wave_length_value;
l=range*wave_length_value;
step=l/(point-1);

i_point=1:point;  %  测试点赋值
matrix_wavelength(i_point)=step*(i_point-1);  

z=linspace(-l,l,100);  %  积分离散

%  求解A向量
for i_point=1:point;
  r=((matrix_wavelength(i_point)-z).^2 + a^2).^(1/2);  %  场源距离离散
  g=exp(-j*k*r)./r;  %  格林函数离散
  A(i_point)=trapz(z,cos(k*z).*g);  %  A元素确定
end

%  求解B向量
for i_point=1:point;
  r=((matrix_wavelength(i_point)-z).^2 + a^2).^(1/2);  %  场源距离离散
  g=exp(-j*k*r)./r;  %  格林函数离散
  B(i_point)=trapz(z,sin(2*k*abs(z)).*g);  %  B元素确定
end

%  求解C向量
for i_point=1:point;
   C(i_point)=cos(k*matrix_wavelength(i_point)); %  B元素确定
end

%  阻抗矩阵确定
impedance_matrix=[A.',B.',C.'];

%  电压矩阵确定
for i_point=1:point;
   voltage_matrix(i_point)=(-j*v0/60)*sin(k*abs(matrix_wavelength(i_point))); %  B元素确定
end

%  求解a1,a2,C_contant;
current=impedance_matrix\voltage_matrix';

%  图形表示
z_distribute=linspace(0,l,100);
current_function=current(1,1)*sin(k*(l-abs(z_distribute))) +...
                 current(2,1)*sin(2*k*(l-abs(z_distribute)));  %  离散化电流分布

%  电流实虚部
current_re=real(current_function);
current_im=imag(current_function);


%  绘图
plot(current_re,z_distribute,'r');
hold on;
plot(current_im,z_distribute,'g');

xlabel('current distribution');
ylabel('unitary distance')
title('antenna current distribution plot');
legend('real current','imag current',2);







      





