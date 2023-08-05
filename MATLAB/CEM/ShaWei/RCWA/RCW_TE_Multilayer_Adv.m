%  ������ʵ���ϸ���ϲ��㷨����TE�����䵽���ڹ�դ
%  �ṹ������
%  �ο����ף�M.G.Moharam, Formulation for stable and efficient implementation of
%  the rigorous coupled-wave analysis of binary gratings, J. Opt. Soc. AM. A, 
%  Vol. 12, No. 5, May 1995
%  ����ˣ�ɳ�� (Wei E.I. Sha)����۴�ѧ��2009��11��26��
%  Email: wsha@eee.hku.hk

function RCW_TE_Multilayer_Adv;

clc;clear
theta=0*pi/180:pi/180:70*pi/180;  % ����Ƕ�

for index=1:length(theta)
    disp(length(theta)-index)
    [ret_t,ret_r]=Compute(theta(index));  %  ����0�׷��䴫��ϵ��
    T(index)=(ret_t);
    R(index)=(ret_r);
end
plot(0:1:70,T,'bO-')
hold on;
plot(0:1:70,R,'rO-')


function [ret_t,ret_r]=Compute(theta);

%  ��������
v0=3*10^8;          %  ���ɿռ䲨��
lambda0=1;          %  ���ɿռ䲨��
f=v0/lambda0;       %  Ƶ��
k0=2*pi/lambda0;    %  ���ɿռ䲨��
n1=1;               %  ��һ����ʳ��������ȣ�
n2=1;               %  ���һ����ʳ��������ȣ�
% theta=0*pi/180;     %  ����Ƕ�
Per=0.6;            %  ���ڴ�С
N=21;               %  ��ϲ���Ŀ

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
Z1=diag(kz1/(k0));
Z2=diag(kz2/(k0));

%  �����������
Casade=[eye(N,N);j*Z2];
%  �������
RR=[eye(N,N);-j*Z1];

%  ����м��(�Ӻ���ǰ˳�����һ���Ӧ����ϵ��)
nrd=[3,4];     %  ÿһ��ridge���ʳ���
ngr=[1 1];     %  ÿһ��groove���ʳ���
fra=[0.5 1];   %  ÿһ��ridge��ռ����
d=[0.1 0.2];   %  ÿһ����

T_Norm=eye(N,N);       %  ��������һ��

for index=1:length(d)
    [Casade,ret]=periodic(k0,Kx,cw,nrd(index),ngr(index),fra(index),d(index),Casade);  %  ���ڼ����������
    T_Norm=T_Norm*ret;
end

%  ���վ�������
ZZ=[Casade,-RR];
VV=[1;zeros((N-1),1);j*cos(theta)*n1;zeros((N-1),1)];
Result=ZZ\VV;            %  �������棨����ϵ��R������ϵ��T��
Tra=Result([1:N],1);     %  ��ȡ����ϵ��
Tra=T_Norm*Tra;          %  ����ϵ����һ��
Ref=Result([N+1:2*N],1); %  ��ȡ����ϵ��

DE_T=Tra.*conj(Tra).*real(kz2/(k0*n1*cos(theta))).';  %  ��һ������ϵ��
DE_R=Ref.*conj(Ref).*real(kz1/(k0*n1*cos(theta))).';  %  ��һ������ϵ��

%  ������֤(��������Ľ��ʣ������غ㣬����1)
sum(DE_T+DE_R);

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

A=Kx^2-E;         %  ��ϲ����������������ֵ������������
[W,Qq] = eig(A);  %  ����ֵ������
Q=sqrt(Qq);       %  ���������ľ�����


%  ��������
X=diag(exp(-k0*d*diag(Q)));
V=W*Q;

S=[W,W;V,-V];
U=inv(S)*T;

a=zeros(N,N);
b=zeros(N,N);
a=U([1:N],[1:N]);
b=U([N+1:2*N],[1:N]);

% [Uu,Ss,Vv]=svd(a);
% V_s=diag(Ss);
% V_s=diag(1./V_s);
% a=Vv*V_s*Uu';

a=inv(a);

Ret=a*X;
A=X*b*Ret;
Casade=[W*(eye(N,N)+A);V*(eye(N,N)-A)];






