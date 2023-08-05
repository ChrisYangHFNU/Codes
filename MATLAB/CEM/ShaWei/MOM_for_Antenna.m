%  MOM_for_Antenna.m
%  ����������ʵ���þ������ͺ��׻��ַ������Գ����ߵ����ֲ�
%  �����    ɳ��(Wei Sha) ���մ�ѧ(Anhui University) ws108@ahu.edu.cn


%  wave_length_value ����
%  measurement �����ߴ�
%  v0 ��ѹ����
%  count_point �������
%  range �������߹�һ������
%  a ����뾶
%  k ����
%  l ���߳���
%  step ��ƥ����

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

i_point=1:point;  %  ���Ե㸳ֵ
matrix_wavelength(i_point)=step*(i_point-1);  

z=linspace(-l,l,100);  %  ������ɢ

%  ���A����
for i_point=1:point;
  r=((matrix_wavelength(i_point)-z).^2 + a^2).^(1/2);  %  ��Դ������ɢ
  g=exp(-j*k*r)./r;  %  ���ֺ�����ɢ
  A(i_point)=trapz(z,cos(k*z).*g);  %  AԪ��ȷ��
end

%  ���B����
for i_point=1:point;
  r=((matrix_wavelength(i_point)-z).^2 + a^2).^(1/2);  %  ��Դ������ɢ
  g=exp(-j*k*r)./r;  %  ���ֺ�����ɢ
  B(i_point)=trapz(z,sin(2*k*abs(z)).*g);  %  BԪ��ȷ��
end

%  ���C����
for i_point=1:point;
   C(i_point)=cos(k*matrix_wavelength(i_point)); %  BԪ��ȷ��
end

%  �迹����ȷ��
impedance_matrix=[A.',B.',C.'];

%  ��ѹ����ȷ��
for i_point=1:point;
   voltage_matrix(i_point)=(-j*v0/60)*sin(k*abs(matrix_wavelength(i_point))); %  BԪ��ȷ��
end

%  ���a1,a2,C_contant;
current=impedance_matrix\voltage_matrix';

%  ͼ�α�ʾ
z_distribute=linspace(0,l,100);
current_function=current(1,1)*sin(k*(l-abs(z_distribute))) +...
                 current(2,1)*sin(2*k*(l-abs(z_distribute)));  %  ��ɢ�������ֲ�

%  ����ʵ�鲿
current_re=real(current_function);
current_im=imag(current_function);


%  ��ͼ
plot(current_re,z_distribute,'r');
hold on;
plot(current_im,z_distribute,'g');

xlabel('current distribution');
ylabel('unitary distance')
title('antenna current distribution plot');
legend('real current','imag current',2);







      





