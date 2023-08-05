%  ������ʵ���ϸ���ϲ��㷨����TM�����䵽���ڹ�դ
%  �ṹ������
%  �ο����ף�M.G.Moharam, Formulation for stable and efficient implementation of
%  the rigorous coupled-wave analysis of binary gratings, J. Opt. Soc. AM. A, 
%  Vol. 12, No. 5, May 1995
%  ����ˣ�ɳ�� (Wei E.I. Sha)����۴�ѧ������ӹ���ѧϵ��2009��11��26��
%  Email: wsha@eee.hku.hk

function RCW_TM_Multilayer_Adv;

clc;clear
theta=0:1:89;  %  ����

for index=1:length(theta)
    disp(length(theta)-index)
    [ret_t,ret_r,ret_a]=Compute(theta(index));  %  ����0�׷��䴫��ϵ��
    T(index)=(ret_t);
    R(index)=(ret_r);
    A(index)=(ret_a);
end

plot(R)

function [ret_t,ret_r,ret_a]=Compute(theta);
%  ��������
v0=3*10^8;          %  ���ɿռ䲨��
lambda0=2;          %  ���ɿռ䲨��
f=v0/lambda0;       %  Ƶ��
k0=2*pi/lambda0;    %  ���ɿռ䲨��
n1=1;               %  ��һ����ʳ��������ȣ�
n2=1;               %  ���һ����ʳ��������ȣ�
theta=theta*pi/180;  %  ����Ƕ�
Per=0.6;             %  ���ڴ�С
N=21;                %  ��ϲ���Ŀ
 
%  ��������ֵ����������

%  ��ϲ�����
cw(1)=0;
cw(2:2:N-1)=1:1:(N-1)/2;
cw(3:2:N)=-1:-1:-(N-1)/2;

kx_i=k0*(n1*sin(theta)-cw*(lambda0/Per));  %  �������ڲ���
Kx=diag(kx_i/k0);  %  �����ԽǾ���

%  ������
kz1=zeros(1,N);
kz2=zeros(1,N);
for index=1:N

    %  ������
    kz1(index)=k0*sqrt(n1^2-(kx_i(index)/k0)^2);
    kz2(index)=k0*sqrt(n2^2-(kx_i(index)/k0)^2);
    
    %  �����˲���improper wave��
    if (imag(kz1(index))>0)
        kz1(index)=-kz1(index);
    end
    %  �����˲���improper wave��
    if (imag(kz2(index))>0)
        kz2(index)=-kz2(index);
    end    
end
Z1=diag(kz1/(k0*n1^2));
Z2=diag(kz2/(k0*n2^2));

%  �����������
Casade=[eye(N,N);j*Z2];
%  �������
RR=[eye(N,N);-j*Z1];

%  ����м��(�Ӻ���ǰ˳�����һ���Ӧ����ϵ��)
nrd=[(4.0-0.1*j)^(1/2) (3.0-0.2*j)^(1/2)];   %  ÿһ��ridge���ʳ���
ngr=[1 1];             %  ÿһ��groove���ʳ���
fra=[0.5 1];           %  ÿһ��ridge��ռ����
d=[0.5 0.5];           %  ÿһ����

T_Norm=eye(N,N); %  ��������һ��

for index=1:length(d)
    [Casade,ret]=periodic(k0,Kx,cw,nrd(index),ngr(index),fra(index),d(index),Casade);  %  ���ڼ����������
    T_Norm=T_Norm*ret;
end

%  ���վ�������
ZZ=[Casade,-RR];
VV=[1;zeros((N-1),1);j*cos(theta)/n1;zeros((N-1),1)];
Result=ZZ\VV;            %  �������棨����ϵ��R������ϵ��T��
Tra=Result([1:N],1);     %  ��ȡ����ϵ��
Tra=T_Norm*Tra;          %  ����ϵ����һ��
Ref=Result([N+1:2*N],1); %  ��ȡ����ϵ��

DE_T=Tra.*conj(Tra).*real((n1*kz2/(k0*n2^2*cos(theta)))).';  %  ��һ������ϵ��
DE_R=Ref.*conj(Ref).*real(kz1/(k0*n1*cos(theta))).';  %  ��һ������ϵ��

%  ������֤(��������Ľ��ʣ������غ㣬����1)
sum(DE_T+DE_R);
ret_a=1-sum(DE_T+DE_R);

%  ����ϵ��
ret_t=DE_T(1,1);
ret_r=DE_R(1,1);

%  cw��ϲ�����
%  nrd/ngr--deilectric constant
%  fra--ratio of nrd to Periodic
function [Casade,Ret]=periodic(k0,Kx,cw,nrd,ngr,fra,d,T)

N=length(cw);        %  ��ϲ���Ŀ
four=(1-N):1:(N-1);  %  ����Ҷ����
epr_h=(nrd^2-ngr^2)*sin(pi*four*fra)./(pi*four);  %  ���ڽ��ʳ����ĸ���Ҷ����
epr_h(N)=nrd^2*fra+(1-fra)*ngr^2;                 %  ���ڽ��ʳ�����ƽ��ֵ��ֱ��������

E=zeros(N,N);

%  ���ڽ��ʳ�������
for ii=1:N
    for p=1:N
        E(ii,p)=epr_h(cw(ii)-cw(p)-four(1)+1);
    end
end

% INV_E=inv(E);
xx = E\Kx;
B=Kx*xx-eye(N,N);
% B=Kx*INV_E*Kx-eye(N,N);
A=E*B;            %  ��ϲ����������������ֵ������������

[W,Qq] = eig(A);  %  ����ֵ������
Q=sqrt(Qq);       %  ���������ľ�����

%  ��������
X=diag(exp(-k0*d*diag(Q)));
xx = E\W;
% V=INV_E*W*Q;
V=xx*Q;

S=[W,W;V,-V];
% U=inv(S)*T;
U=S\T;

a=zeros(N,N);
b=zeros(N,N);
a=U([1:N],[1:N]);
b=U([N+1:2*N],[1:N]);

Ret=a\X;

A=X*b*Ret;
Casade=[W*(eye(N,N)+A);V*(eye(N,N)-A)];






