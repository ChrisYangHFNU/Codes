%  �����ݶȿ��ٸ���Ҷ�任���������ַ���
%  ����� ɳ�� ��۴�ѧ������ӹ���ϵ
%  Email: wsha@eee.hku.hk

function fft_2d_VIE;
clc;clear

%  ����
wavelength=1;  %--wavelength
k0=2*pi/wavelength; %--wave number in free space
r=wavelength;  %--radius ****
z=120*pi;      %--wave impedance in free space
epr=8;         % ---permittivity  ****
deltax=wavelength/20/sqrt(epr);  %--grid size  ****�Ƽ�20����1�����ʲ���
a=sqrt(deltax*deltax/pi);        %--effective area
const=-j*z*epr/(k0*(epr-1));

%  ��ģ\�����\���䳡
M=round(2*r/(deltax))+10;  %--��Ҫ����ż������������10��
N=M;
Ein=zeros(M,N);  %--���䳡����
model=zeros(M,N);  %--�ж��ǲ��ǽ��ʻ��ǿ���
Z=zeros(M,N);  %--�迹����
xc=-M/2*deltax+deltax/2;
yc=N/2*deltax-deltax/2;
epr_arr=const*ones(M,N);  %--��糣������

for f_n=1:N;  %  ��
    for f_m=1:M;  %  ��
        x=xc+(f_n-1)*deltax;
        y=yc-(f_m-1)*deltax;
        flag=cylinder(x,y,r);  %  �ǲ���Բ���ڲ�����ģ��Ҫ�޸�cylinder������

        if (flag)
            model(f_m,f_n)=1;           %  ��������1������0
            Ein(f_m,f_n)=exp(-j*k0*x);  %  ����糡
        end
        
        R=sqrt((x-xc)^2+(y-yc)^2);
        Z(f_m,f_n)=z*pi*a/2*besselj(1,k0*a)*besselh(0,2,k0*R); % impedance matrix
    end
end
Z(1,1)=z*pi*a/2*besselh(1,2,k0*a);

%  ��������
Zf1=Z(:,end:-1:2);
Zf2=Z(end:-1:2,:);
Zf3=Zf2(:,end:-1:2);
Z=[Z,Zf1;Zf2,Zf3];

%  ������
J=zeros(M,N);  % current
r0=-Ein;
p1=-fft_mv_trans(Z,r0,M,N,model,const,epr_arr);
q1=-p1;
change3=q1.*conj(q1);
con3=sum(sum(change3));
changex=Ein.*conj(Ein);
conx=sum(sum(changex));

%  �����ݶ�--���ٸ���Ҷ�㷨
for n=1:10^8
    %  A^{*}r_{n-1}
    change1=change3;
    con1=con3;
    %  Ap_{n}
    Ap=fft_mv(Z,p1,M,N,model,const,epr_arr);
    change2=Ap.*conj(Ap);
    con2=sum(sum(change2));
    %  alpha
    alpha=con1/con2;
    %  Update J
    J=J+alpha*p1;
    %  Update r
    r0=r0+alpha*Ap;
    change4=r0.*conj(r0);
    
    %  Beta
    q1=fft_mv_trans(Z,r0,M,N,model,const,epr_arr);
    change3=q1.*conj(q1);
    con3=sum(sum(change3));
    beta=con3/con1;
    %  p1
    p1=-q1+beta*p1;
    %  �ж�
    con4=sum(sum(change4))/conx;
    
    %  ÿ��ʮ����ʾ���
    if (mod(n,10)==0)
        con4
    end
    
    %  error truncation
    if con4<1*10^(-6)
        n  %һ�����ٲ�����
        break;
    end
end

figure(1)
pcolor(abs(J));
shading interp
axis equal
colorbar

%  RCS����
t=0;  %  �����ʹ��������ܳ�
for s=0:pi/180:pi
    t=t+1;
    sum1=0;
    for f_n=1:N;  %  ��
        for f_m=1:M;  %  ��
            x=xc+(f_n-1)*deltax;
            y=yc-(f_m-1)*deltax;
            sum1=sum1+J(f_m,f_n)*exp(j*k0*x*cos(s)+j*k0*y*sin(s));  %  ��������
        end;
    end
    rcs(t)=((k0/4)*(abs(sum1).^2))*(a*wavelength*z*besselj(1,k0*a))^2; %  ��һ�� 
end;

rcs=10.*log10(rcs/wavelength);

%  ������
t=0;
a=r;              % radius of the cylinder
Z0=z;             % wave impedance in free space
Z1=Z0/sqrt(epr);  % wave impedance in dielectric
k1=k0*sqrt(epr);  % wave number in dielectric
NN=round(2*(k1*a)); % �ض�����

for angle=0:pi/180:pi;   %  ɢ���
    
    sumx=0;
    t=t+1;
  
    for n=-NN:NN;
        
        cons1=(1/2)*(besselj(n-1,k1*a)-besselj(n+1,k1*a));      %  ������������
        cons2=(1/2)*(besselh(n-1,2,k0*a)-besselh(n+1,2,k0*a));  %  ������������ 
        cons3=(1/2)*(besselj(n-1,k0*a)-besselj(n+1,k0*a));      %  ������������
        
        d1=Z0*besselj(n,k0*a)*cons1-Z1*besselj(n,k1*a)*cons3;
        d2=Z1*besselj(n,k1*a)*cons2-Z0*besselh(n,2,k0*a)*cons1;
        A=d1/d2;  %  ϵ��An

        sumx=sumx+A*exp(j*n*angle);  %  ���Ա���չ�����ɢ�䳡
    end
    
    es(t)=abs(sumx).^2*4/k0;
end;

RCS=10*log10(es/wavelength);  %  ��һ��

figure(2);
hold on
plot(0:1:180,rcs,'b');
plot(0:1:180,RCS,'rx-');
xlabel('Angle (degree)')
ylabel('RCS (dB)')
legend('VIE-FFT','Analytical')
title('˫վ�״�ɢ�����');

%  Z*j (FFT ʵ�־���������)
function current=fft_mv(Z,J,M,N,model,const,epr_arr)
current=ifft2(fft2(Z).*fft2(J,2*M-1,2*N-1));    %  ���ټ���
current=current([1:M],[1:N]);                   %  FFT������Ч
current=epr_arr.*J+current;                     %  D*J+K**J
current=current.*model;                         %  ģ�ͣ����㣩

%  Z'*j  (FFT ʵ�־�����ת��������)
function current=fft_mv_trans(Z,J,M,N,model,const,epr_arr)
current=ifft2(fft2(conj(Z)).*fft2(J,2*M-1,2*N-1));  %  ����ת��;
current=current([1:M],[1:N]);                       %  FFT������Ч
% current=current.';                                %  ת��
current=current+conj(epr_arr).*J;                   %  D*J+K**J
current=current.*model;                             %  ģ�� �����㣩  

% cylinder ���ж��ǲ���Բ���ڲ���
function flag=cylinder(x,y,r)
flag=0;
if (x^2+y^2<=r^2)
    flag=1;
end

 
