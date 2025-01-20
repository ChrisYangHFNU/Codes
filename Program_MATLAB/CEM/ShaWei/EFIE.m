%  EFIE.m
%  ���������ڵ糡���ַ������Բ��ɢ��
%  �����    ɳ��(Wei Sha) ���մ�ѧ(Anhui University) ws108@ahu.edu.cn


tic
clc;clear;

%  ��������
f=500*10^6;  %  TM ��Ƶ��
v=3*10^8;    %  ����
S0=1/(4*pi*9*10^9); %  ��ս�糣��
U0=4*pi*10^(-7);    %  �ŵ���
Z=sqrt(U0/S0);      %  ����в��迹
wavelength=v/f;     %  ����
k=2*pi/wavelength;  %  ����
EULER=1.78107;      %  ŷ������
EX=exp(1);          %  ָ������


amplitude=1;         %  �糡���
a=wavelength;        %  Բ���뾶��Ϊ���˼�����ȡ k*a=1
anglemin=pi/80;      %  ��С��ɢ�ǶȲ���
count=2*pi/anglemin; %  ѡ�����

circlelength=a*anglemin;  %  ����

%  ���궨λ(�Ƕ�)
anglevalue=-anglemin/2;  %  �Ƕȳ�ֵ
for i=1:count
    anglevalue=anglevalue+anglemin;  %  �Ƕȵ���
    point(i)=anglevalue;  %  �Ƕȸ�ֵ
end



%  �迹����l_mnȷ��
for m=1:count  %  ����ѭ��
    for n=1:count  %  Դ��ѭ��  
        if (m==n)
            l_mn(m,n)=(1-(2*1i/pi)*log(EULER*k*circlelength/(4*EX)))*circlelength*(k*Z/4);
            %  ����㴦��
        else
            fieldpoint=point(m);   %  ����
            currentpoint=point(n); %  Դ��
   
            % %  ��������
            integerfunction1=2*k*a*abs(sin((currentpoint-fieldpoint)/2));
            integerfunction2=besselj(0,integerfunction1)-1i*bessely(0,integerfunction1);
            % % Ԫ�ظ�ֵ
            l_mn(m,n)=(circlelength*k*Z/4)*integerfunction2; %  ���û��ֺ�������迹����l_mn
        end  
    end  
end 
     
  

%  ��ѹ����g_mȷ��
for ii=1:count
    g_m(ii)=amplitude*exp(-1i*k*a*cos(point(ii)) );  %  ��X���򴫲���TM��
end



%  �������ֲ�
currentresult=l_mn\g_m.';  %  �������



%  ���������
initialcurrent=0;  %  ��ͳ�ʼ��
for m=1:count
    for n=-20:20
        initialcurrent=initialcurrent+(2*amplitude/(2*pi*f*U0*pi*a))*...
                       1i^(-n)*exp(1i*n*point(m))/(besselj(n,k*a)-1i*bessely(n,k*a));  %  ���������
    end
        analysis(m)=initialcurrent;  %  ��ֵ
        initialcurrent=0;  %  ����
end
    


%   RCS���
ii=0;  %  �����±�
for si=0:pi/180:2*pi  %  ɢ���
    
    ii=ii+1;
    sum=0; %  ��ͳ�ֵ
    
    for seg=1:count
       sum=sum+currentresult(seg,1)*circlelength*exp(1i*k*( a*cos(point(seg))*cos(si)+a*sin(point(seg))*sin(si) ) );
    end
    
    rcs(ii)=abs(sum).^2*k*Z^2/4;  %  RCS����
end

rcs=10*log10(rcs*k/(2*pi));  %  ��һ��RCS



%  ͼ����ʾ
figure(1)
plot(point*180/pi,abs(currentresult.'),'rx');  %  ������������ĵ���
hold on;
plot(point*180/pi,abs(analysis),'b');  %  ������������ĵ���

% axis([0,360,min(abs(analysis))-10,max(abs(analysis))+10]); %  ���������
title('MOM APPLIED IN CYLINDER SCATTERING PROBLEM');  %  ��Ŀ
xlabel('angle');  %  x������
ylabel('current');  %  y������
legend('mom result','analysis result');  %  ��ʶ
hold off
        
figure(2)
plot([0:1:180],rcs(1:181));
xlabel('scattered angle��/^{��}��');  %  y������
ylabel('Radar Cross Section /dB');  %  x������
title('RCS OF CYLINDER WITH SIDE-DIMENSION OF DOUBLE WAVELENGTHS');  %  ��Ŀ


toc;

