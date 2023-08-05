%  TM_MEDIA.m

%  ���������ڵ�Чԭ��;�����������Բ�������ʹ����ֲ�����˫վ�״�ɢ����档
%  ���������䲨ΪTM��,Z���򼫻���X���򴫲���������������Z�ᡣ
%  ���ڵ�Чԭ��,���������Ч�໥��Ϻ͵糡����(EFIE)������⡣
%  �˳���û�п�����г�����⣬����������г�񣬲�����ϳ�����(CFIE)���̻�PMCHWT���̴���
%  �˳�����á���һ�����������������˽ض��������ڲ��������ע�⡣J��=J*Z0;J��=J��*Z1/Z0��
%  ������С��ɢ���ȣ�ÿ����ȡ10���㣬��ע�����Ⲩ����ͬ��Ҫ����С����Ϊ��׼��
%  �˳�����Խ�糣����epsilon����������ʱ�����ʱ�ת��Ϊ�����ڡ���������֤����
%  �˳�����Դŵ��ʣ�u_r����������ʱ,���ʱ�ת��Ϊ�űڡ������ɽ����ڣ��൱��TE���������(������֤����
%  �˳�����Խ�糣����epsilon������Դŵ��ʣ�u_r������Ϊ���������˳�����Լ����кĽ���Բ����ɢ�䡣
%  �����ɳ�������մ�ѧ��ų���΢��רҵ
%  Email: dr.weisha@gmail.com
%  ���ʱ�� 2004-5-18

tic;
clc;clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1.��������

f=300*10^6;                     %  TM ��Ƶ��
v=3*10^8;                       %  ����
S0=1/(4*pi*9*10^9);             %  ��ս�糣��
U0=4*pi*10^(-7);                %  ��մŵ���
EULER=1.78107;                  %  ŷ������
EX=exp(1);                      %  ָ������
anglemin=pi/80;                 %  ��С��ɢ�ǶȲ���
amplitude=1;                    %  �糡���

wavelength=v/f;                 %  ��ղ���
k0=2*pi/wavelength;             %  ��ղ���
a=wavelength;                   %  �������뾶
count=2*pi/anglemin;            %  ѡ�����
circlelength=a*anglemin;        %  ����
Z0=sqrt(U0/S0);                 %  ����в��迹

epsilon=4;                      %  ��������Խ�糣��
u_r=1;                          %  ��������Դŵ���
Z1=sqrt(u_r)*Z0/sqrt(epsilon);  %  ���������迹
k1=k0*sqrt(epsilon*u_r);        %  ����������


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2.���궨λ

%  Բ�����Ľ�������ϵ


%  Բ������
anglevalue=-anglemin/2;              %  �Ƕȳ�ֵ

for m=1:count;
    anglevalue=anglevalue+anglemin;  %  �Ƕȵ���
    point(m)=anglevalue;             %  �Ƕȸ�ֵ
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  3.�迹����ֵ


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  A.��������Բ��

%(1)L0��ֵ,0��ʾ0���򣬼����Ч������ͬ����

for m=1:count;                     %  ��ѭ��
    for n=1:count;                 %  Դѭ��
        
        field0=point(m);           %  ����λ��
        source0=point(n);          %  Դ��λ��
        
        if (m==n)
            L0(m,n)=(k0/4)*circlelength*( 1-(j*2/pi)*log(EULER*k0*circlelength/(4*EX)) );
        else
            quad_value=2*k0*a*abs( sin( (field0-source0)/2 ) );
            L0(m,n)=(k0/4)*circlelength*(besselj(0,quad_value)-j*bessely(0,quad_value));
        end;
        
    end;
end;


%(2)M0��ֵ

for m=1:count;                     %  ��ѭ��
    for n=1:count;                 %  Դѭ��

        field0=point(m);           %  ����λ��
        source0=point(n);          %  Դ��λ��
        
        if (m==n)                  
            M0(m,n)=1/2;           %  ��㴦��
        else
            quad_value=2*k0*a*abs( sin( (field0-source0)/2 ) );
            M0(m,n)=-(j*k0/4)*circlelength*(besselj(1,quad_value)-j*bessely(1,quad_value))*...
                      abs( sin( (field0-source0)/2 ) );  %  �Ƕ�N.R=|sin((si-si')/2)|
        end;
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  B.��������Բ��

%(1)L1��ֵ,1��ʾ1����Ϊ�ڵ�Ч������ͬ����

for m=1:count;                     %  ��ѭ��
    for n=1:count;                 %  Դѭ��
        
        field1=point(m);           %  ����λ��
        source1=point(n);          %  Դ��λ��
         
        if (m==n)
            L1(m,n)=(Z1/Z0)*(k1/4)*circlelength*( 1-(j*2/pi)*log(EULER*k1*circlelength/(4*EX)) );
        else
            quad_value=2*k1*a*abs( sin( (field1-source1)/2 ) );
            L1(m,n)=(Z1/Z0)*(k1/4)*circlelength*(besselj(0,quad_value)-j*bessely(0,quad_value));
        end;
        
    end;
end;


%(2)M1��ֵ

for m=1:count;                     %  ��ѭ��
    for n=1:count;                 %  Դѭ��

        field1=point(m);           %  ����λ��
        source1=point(n);          %  Դ��λ��
        
        if (m==n)                 
            M1(m,n)=-1/2;          %  ��㴦��
        else
            quad_value=2*k1*a*abs( sin( (field1-source1)/2 ) );
            M1(m,n)=-(j*k1/4)*circlelength*(besselj(1,quad_value)-j*bessely(1,quad_value))*...
                      abs( sin( (field1-source1)/2 ) );    %  �Ƕ�N.R=|sin((si-si')/2)|
        end;
    end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  C.�迹����ϲ�

IMPEDANCE_MATRIX=[ L0 ,M0 ;
                   L1 ,M1 ];                           
            
                            
                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  4.��������ֵ

for t=1:count
    E(t)=amplitude*exp(-j*k0*a*cos(point(t)));  %  �糡
end;


%%%%%����ϲ�
ZERO_MATRIX=zeros(1,count);
EXCITATION_MATRIX=[E,ZERO_MATRIX].';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  5.�������ʹ���

RESULT=IMPEDANCE_MATRIX\EXCITATION_MATRIX;

figure(1);
plot(abs(RESULT(1:count)),'r');            %  ����
hold on;
plot(abs(RESULT(count+1:2*count)),'g');    %  ����
legend('����','����');
title('������ֲ�');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  6. RCS����
 
t=0;  %  �����ʹ��������ܳ�
for s=0:pi/60:pi
    
    t=t+1;
    sum=0;
    
    for r=1:count;
        coeff=cos(s-point(r));  %  �������ʱ�����ĽǶ�
        sum=sum+(RESULT(r)-RESULT(r+count)*coeff)*circlelength*...
                exp(j*k0*a*cos(point(r))*cos(s)+j*k0*a*sin(point(r))*sin(s));  
                %  ����������һ���������˫վ�״�ɢ�����
    end;
    
    rcs(t)=((k0/4)*(abs(sum).^2))/wavelength; %  ��һ�� 
end;

rcs=10.*log10(rcs);

figure(2);
plot(0:pi/60:pi,rcs,'b');
title('���������������˫վ�״�ɢ�����');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  7.������  (���ֵ����ҵ�ų�283ҳ)

t=0;
for angle=0:pi/60:pi;   %  ɢ���
    
    sum=0;
    t=t+1;
    
    for n=-30:30;
        
        cons1=(1/2)*(besselj(n-1,k1*a)-besselj(n+1,k1*a));      %  ������������
        cons2=(1/2)*(besselh(n-1,2,k0*a)-besselh(n+1,2,k0*a));  %  ������������ 
        cons3=(1/2)*(besselj(n-1,k0*a)-besselj(n+1,k0*a));      %  ������������
        
        d1=Z0*besselj(n,k0*a)*cons1-Z1*besselj(n,k1*a)*cons3;
        d2=Z1*besselj(n,k1*a)*cons2-Z0*besselh(n,2,k0*a)*cons1;
        A=d1/d2;  %  ϵ��An

        sum=sum+A*exp(j*n*angle);  %  ���Ա���չ�����ɢ�䳡
    end;
    
    es(t)=abs(sum).^2*4/k0;
end;

RCS=10*log10(es/wavelength);  %  ��һ��

figure(3);
plot(0:pi/60:pi,RCS,'R');
hold on;
plot(0:pi/60:pi,rcs,'bx');
legend('������','������');
title('˫վ�״�ɢ�����Ա�');


