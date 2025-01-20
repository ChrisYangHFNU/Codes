%  本程序实现严格耦合波算法分析TM波入射到周期光栅
%  结构任意多层
%  参考文献：M.G.Moharam, Formulation for stable and efficient implementation of
%  the rigorous coupled-wave analysis of binary gratings, J. Opt. Soc. AM. A, 
%  Vol. 12, No. 5, May 1995
%  编程人：沙威 (Wei E.I. Sha)，香港大学电机电子工程学系，2009年11月26日
%  Email: wsha@eee.hku.hk

function RCW_TM_Multilayer_Adv;

clc;clear
theta=0:1:89;  %  波长

for index=1:length(theta)
    disp(length(theta)-index)
    [ret_t,ret_r,ret_a]=Compute(theta(index));  %  计算0阶反射传输系数
    T(index)=(ret_t);
    R(index)=(ret_r);
    A(index)=(ret_a);
end

plot(R)

function [ret_t,ret_r,ret_a]=Compute(theta);
%  参数定义
v0=3*10^8;          %  自由空间波速
lambda0=2;          %  自由空间波长
f=v0/lambda0;       %  频率
k0=2*pi/lambda0;    %  自由空间波数
n1=1;               %  第一层介质常数（均匀）
n2=1;               %  最后一层介质常数（均匀）
theta=theta*pi/180;  %  入射角度
Per=0.6;             %  周期大小
N=21;                %  耦合波数目
 
%  周期特征值和特征向量

%  耦合波向量
cw(1)=0;
cw(2:2:N-1)=1:1:(N-1)/2;
cw(3:2:N)=-1:-1:-(N-1)/2;

kx_i=k0*(n1*sin(theta)-cw*(lambda0/Per));  %  切向周期波数
Kx=diag(kx_i/k0);  %  参数对角矩阵

%  纵向波数
kz1=zeros(1,N);
kz2=zeros(1,N);
for index=1:N

    %  纵向波数
    kz1(index)=k0*sqrt(n1^2-(kx_i(index)/k0)^2);
    kz2(index)=k0*sqrt(n2^2-(kx_i(index)/k0)^2);
    
    %  不适宜波（improper wave）
    if (imag(kz1(index))>0)
        kz1(index)=-kz1(index);
    end
    %  不适宜波（improper wave）
    if (imag(kz2(index))>0)
        kz2(index)=-kz2(index);
    end    
end
Z1=diag(kz1/(k0*n1^2));
Z2=diag(kz2/(k0*n2^2));

%  级联传输矩阵
Casade=[eye(N,N);j*Z2];
%  反射矩阵
RR=[eye(N,N);-j*Z1];

%  求解中间层(从后向前顺序，最后一层对应传输系数)
nrd=[(4.0-0.1*j)^(1/2) (3.0-0.2*j)^(1/2)];   %  每一层ridge介质常数
ngr=[1 1];             %  每一层groove介质常数
fra=[0.5 1];           %  每一层ridge所占比例
d=[0.5 0.5];           %  每一层厚度

T_Norm=eye(N,N); %  传输矩阵归一化

for index=1:length(d)
    [Casade,ret]=periodic(k0,Kx,cw,nrd(index),ngr(index),fra(index),d(index),Casade);  %  周期级联传输矩阵
    T_Norm=T_Norm*ret;
end

%  最终矩阵生成
ZZ=[Casade,-RR];
VV=[1;zeros((N-1),1);j*cos(theta)/n1;zeros((N-1),1)];
Result=ZZ\VV;            %  矩阵求逆（传输系数R，反射系数T）
Tra=Result([1:N],1);     %  抽取传输系数
Tra=T_Norm*Tra;          %  传输系数归一化
Ref=Result([N+1:2*N],1); %  抽取反射系数

DE_T=Tra.*conj(Tra).*real((n1*kz2/(k0*n2^2*cos(theta)))).';  %  归一化传输系数
DE_R=Ref.*conj(Ref).*real(kz1/(k0*n1*cos(theta))).';  %  归一化反射系数

%  收敛验证(对于无损耗介质，能量守恒，等于1)
sum(DE_T+DE_R);
ret_a=1-sum(DE_T+DE_R);

%  返回系数
ret_t=DE_T(1,1);
ret_r=DE_R(1,1);

%  cw耦合波向量
%  nrd/ngr--deilectric constant
%  fra--ratio of nrd to Periodic
function [Casade,Ret]=periodic(k0,Kx,cw,nrd,ngr,fra,d,T)

N=length(cw);        %  耦合波数目
four=(1-N):1:(N-1);  %  傅里叶级数
epr_h=(nrd^2-ngr^2)*sin(pi*four*fra)./(pi*four);  %  周期介质常数的傅里叶级数
epr_h(N)=nrd^2*fra+(1-fra)*ngr^2;                 %  周期介质常数的平均值（直流分量）

E=zeros(N,N);

%  周期介质常数矩阵
for ii=1:N
    for p=1:N
        E(ii,p)=epr_h(cw(ii)-cw(p)-four(1)+1);
    end
end

% INV_E=inv(E);
xx = E\Kx;
B=Kx*xx-eye(N,N);
% B=Kx*INV_E*Kx-eye(N,N);
A=E*B;            %  耦合波矩阵（用于求解特征值和特征向量）

[W,Qq] = eig(A);  %  特征值和向量
Q=sqrt(Qq);       %  特征向量的均方根

%  参数矩阵
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






