% 应用反应积分方程(1)，采用分段正弦基函数伽辽金法获取对称振子的输入阻抗，沿线电流分布与方向图
clear;clc;
global beta a zzm zzn Deltz;
C = 3.0e8;
f = 3.0e8;
beta = 2*pi*f/C; % 相位常数
a = 0.0005; % 导线半径
L = [0.1:.01:2.0]; % 对称振子电长度 L=2l ==>0.1*lambda~2*lambda，间隔0.01*lambda
len = length(L);
ZA = zeros(1,len); % 输入阻抗向量初始值置零

for num=1:len
    if L(1,num)<=1.0  % 如果L<=1.0*lambda
       N = 11; % 振子均匀分割为12个子段，11个分段正弦基函数跨越这些子段
    else
       N = 23;  % 振子均匀分割为24个子段，23个分段正弦基函数跨越这些子段
    end

    Deltz = L(1,num)/(N+1); %子段长度
    zm = [(N-1)/2:-1:-(N-1)/2]*Deltz; % 振子中心位于坐标原点，沿z轴放置，分段节点坐标（不含端点）
    % =======================================================
    % 获取广义互阻抗矩阵[Zmn]
    Zmn = zeros(N,N); % [Zmn]初始化
    zzm = zm(1,1); % 不含端点，最上端分段节点坐标

    for n=1:N
        zzn = zm(1,n); % 不含端点，从上到下，各分段节点坐标
        Zmn(1,n) = -1j*30/(sin(beta*Deltz)^2)*(integral(@Fun_Zmn01,zzn,zzn+Deltz)+integral(@Fun_Zmn02,zzn-Deltz,zzn));
    end

    Zmn = toeplitz(Zmn(1,:),Zmn(1,:)); % [Zmn]为Toeplitz矩阵
    % ============================
    % 获取广义电压向量[Vm]
    VA = 1.0; % 外加激励电压
    Vm = zeros(N,1); % [Vm]初始化
    % 缝隙δ激励模型(slice generator)
    Vm((N+1)/2,1) = -1/sin(beta*Deltz)*(VA/Deltz)*(integral(@Fun_Vm01,0,Deltz)+integral(@Fun_Vm02,-Deltz,0));
    % ===============================
    % 获取广义电流向量 [In]
    In = zeros(N,1); % [In]初始化
    In = Zmn\Vm;
    % ===============================
    ZA(1,num) = VA/In((N+1)/2,1);
end

figure(1)
plot(L,real(ZA),'r.-'); % 绘制输入电阻随对称振子电长度L的变化曲线
axis([0 2.0 0 2500]);
hold on;

figure(2)
plot(L,imag(ZA),'b.-'); % 绘制输入电抗随对称振子电长度L的变化曲线
axis([0 2.0 -1500 1200]);
hold on;
% =================================
% 电长度L=2.0时，绘制对称振子沿线电流分布图
% 绘制构成振子沿线电流分布的分段正弦基函数图
pw = 50; % 振子均匀分割为N+1子段，各子段均匀细分为50分段
Sn = zeros((N+1)*pw,N); % 为便于获得振子沿线电流分布，特构造分段正弦电流基函数矩阵，初值置零
for n=1:N
    cn = zm(1,n); % 第n节点坐标Zn
    zp = linspace(cn-Deltz,cn+Deltz,2*pw)'; % 含端点，从下到上，各子段细分节点坐标向量
    % 填充电流基函数矩阵第n列
    Sn((n-1)*pw+1:(n+1)*pw,n) = In(n,1)*sin(beta*(Deltz-abs(zp-cn)))/sin(beta*Deltz);
end
% ====================================

figure(3)  % 绘制分段正弦电流基函数图
for n=1:N
    plot([1:(N+1)*pw]*(Deltz/pw),abs(Sn(:,n)),'k-'); % 绘制第n分段正弦电流基函数图
    hold on;
end
I = zeros((N+1)*pw,1); % 对称振子沿线电流分布向量，初值置零
I = sum(Sn,2); % 由分段正弦电流基函数矩阵，获取对称振子沿线电流分布向量
plot([pw:N*pw]*(Deltz/pw),abs(I(pw:N*pw)),'r.-'); % 绘制对称振子沿线电流分布图
set(gca,'XDir','reverse'); % 横坐标反转
% ====================================

% 电长度L=2.0时，绘制对称振子方向图，通过比较，说明传输线模型（正弦驻波电流分布）是对称振子远场特性分析（方向图）的良好近似
Theta = 0:pi/180:2*pi; % 为绘制E面方向图，取子午角 0<theta<2*pi
Len = length(Theta);
F = zeros(1,Len); % 矩量法方向函数向量，初值置零
Fp = zeros(1,Len); % 传输线模型方向图函数向量，初值置零
% 获取矩量法方向函数向量，为此，将原对称振子视为由分段正弦电流基函数短对称振子构成的不均匀直线阵列，其方向函数的求取参见教材
for n=1:N % 计算第n分段正弦电流基函数短对称振子的方向函数值
    cn = zm(1,n); % 第n分段正弦电流基函数短对称振子的中点（视为馈电点）坐标
    for m=1:Len % 计算第n分段正弦电流基函数短对称振子在Theta处的方向函数值
        theta = (m-1)*pi/180;
        if sin(theta)~=0 % 此式隐含： 若sin(theta)=0,则F(theta)=0
        F(1,m) = F(1,m)+1j*(60/sin(beta*Deltz))*((cos(beta*Deltz*cos(theta))-cos(beta*Deltz))/sin(theta))*In(n,1)*exp(1j*beta*cn*cos(theta));
        end
    end
end
% ===== 获取传输线模型方向函数向量，参见教材(l=L/2=1.0lambda) ======
for m=1:Len
    theta = (m-1)*pi/180;
    if sin(theta)~=0
       Fp(1,m) = Fp(1,m)+(cos(beta*cos(theta))-cos(beta))/sin(theta);
    end
end
figure(4)
polar(Theta,abs(Fp)/max(abs(Fp)),'ko'); % 绘制传输线模型归一化方向图
hold on;
polar(Theta(1,1:3:Len),abs(F(1,1:3:Len))/max(abs(F)),'r-'); % 绘制矩量法方向图
hold on;

function y = Fun_Zmn01(z)
global beta a zzm zzn Deltz;
j = sqrt(-1);
R0 = sqrt((z-zzm).^2+a^2);
R1 = sqrt((z-(zzm-Deltz)).^2+a^2);
R2 = sqrt((z-(zzm+Deltz)).^2+a^2);
y = sin(beta*(zzn+Deltz-z)).*(exp(-j*beta*R1)./R1-2*cos(beta*Deltz)*exp(-j*beta*R0)./R0+exp(-j*beta*R2)./R2);
end

function y = Fun_Zmn02(z)
global beta a zzm zzn Deltz;
j = sqrt(-1);
R0 = sqrt((z-zzm).^2+a^2);
R1 = sqrt((z-(zzm-Deltz)).^2+a^2);
R2 = sqrt((z-(zzm+Deltz)).^2+a^2);
y = sin(beta*(z-(zzn-Deltz))).*(exp(-j*beta*R1)./R1-2*cos(beta*Deltz)*exp(-j*beta*R0)./R0+exp(-j*beta*R2)./R2);
end

function y = Fun_Vm01(z)
global beta a Deltz;
j = sqrt(-1);
y = sin(beta*(Deltz-z));
end

function y = Fun_Vm02(z)
global beta a Deltz;
j = sqrt(-1);
y = sin(beta*(z-(-Deltz)));
end