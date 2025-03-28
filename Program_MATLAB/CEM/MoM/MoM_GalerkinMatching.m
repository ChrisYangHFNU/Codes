clear;clc;
n = input('Input the N \n');
%Galerkin Method
syms x
g = 1 + 4*x^2; % 已知源
fa = 5/6*x -1/2*x^2 -1/3*x^4;    %实际解
% 分配内存
% fn = zeros(1,n);
% wn = zeros(1,n);
% lfn = zeros(1,n);
% lmn = zeros(n,n);
% gm = zeros(1,n);
for ii = 1:n
    fn(ii) = x - x^(ii+1); % 基函数构造，利用边界条件
    wn(ii) = fn(ii);              %伽略金法检验函数等于基函数
    lfn(ii) = -diff(diff(fn(ii))); % 将算子作用到由基函数线性组合而成的原函数近似值上
end
for ii = 1:n
    for jj = 1:n
        lmn(ii,jj) = int(wn(ii)*lfn(jj),x,0,1); % 求内积
    end
    gm(ii) = int(wn(ii)*g,x,0,1); % 求内积
end
a = lmn\gm'; % 矩阵左除
f = fn*a;   % 因为矩阵形状所以相应矩阵乘法调换顺序
val = 0:0.01:1;
plot(val,subs(f,val),'p') % 矩量法解的图形
hold on
plot(val,subs(fa,val),'r') % 解析解图形
str=['N = ',num2str(n)];
title(str);
legend('矩量法解','精确值');