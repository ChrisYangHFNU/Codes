%  MFIE.m
%  ���������ڴų����ַ������Բ��ɢ��
%  �����    ɳ��(Wei Sha) ���մ�ѧ(Anhui University) ws108@ahu.edu.cn

clear;clc;
%%%% 1.����

f=300*10^8;                 %  TM ��Ƶ��
v=3*10^8;                   %  ����
anglemin=pi/128;            %  ��С��ɢ�ǶȲ���
amplitude=1;                %  �糡���

wavelength=v/f;             %  ��ղ���
k0=2*pi/wavelength;         %  ��ղ���
a=4*wavelength;             %  Բ���뾶
count=2*pi/anglemin;        %  ѡ�����
circlelength=a*anglemin;    %  ����

%%%%  2.���궨λ

%  Բ�����Ľ�������ϵ
%  Բ������

anglevalue=-anglemin/2;  %  �Ƕȳ�ֵ

for m=1:count;
    anglevalue=anglevalue+anglemin;  %  �Ƕȵ���
    point(m)=anglevalue;             %  �Ƕȸ�ֵ
end;



%%%%  3.����ֵ

for m=1:count;                     %  ��ѭ��
    for n=1:count;                 %  Դѭ��
        
        field=point(m);            %  ����λ��
        source=point(n);           %  Դ��λ��
        
        if (m==n)
            L(m,n)=-1/2;           %  ����㴦��
        else
            quad_value=2*k0*a*abs( sin( (field-source)/2 ) );
            L(m,n)=(j*k0/4)*circlelength*(besselj(1,quad_value)-j*bessely(1,quad_value))*...
                (abs( sin( (field-source)/2 ) ));
        end;
    end;
end;

for t=1:count;  %  ��ѹ����
    y(t)=-amplitude*exp(-j*a*cos(point(t))*k0);
end;



%%%%  4.����

result=L\y.';
figure(1);
plot(abs(result));
title('����������');



%%%%  5.rcs����

t=0;
for s=0:pi/180:pi
    t=t+1;
    sum=0;
    for r=1:count;
        nr=cos(point(r)-s);  %  N.R,��Զ����Դ�㷽��
        sum=sum+nr*result(r,1)*circlelength*...
            exp(j*k0*a*cos(point(r))*cos(s)+j*k0*a*sin(point(r))*sin(s));
    end;
    rcs1(t)=(k0/4)*abs(sum).^2;
end;

rcs=10.*log10(rcs1/wavelength);  %  ������һ��

figure(2);
plot(rcs);
title('˫վ�״�ɢ�����');



% %%%%  6.���������
% 
% for m=1:count;
%     sum=0;
%     for n=-30:30;
%         t=(1/2)*(besselh(n+1,2,k0*a)-besselh(n-1,2,k0*a));  %  bessel������
%         sum=sum+((-j*2*amplitude)/(pi*k0*a))*j^(-n)*exp(j*n*(point(m)))/t;
%     end;
%     current(m)=sum;
% end;
% 
% figure(3);
% title('�����Ա�');
% plot(abs(current),'r');
% title('�����ֲ��Ա�');
% hold on;
% plot(abs(result),'bx');
% legend('����','����')
% 
% 
% 
% %%%%  7.������RCS
% 
% RR=10.^6;  %  �㹻���R
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
% title('�״�ɢ�����Ա�');
% plot(RCS,'r');
% hold on;
% plot(rcs,'bx');
% legend('������','������');
