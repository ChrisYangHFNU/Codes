clear;clc
n = input('Input the N...\n');
%Point matching
syms x;
g = 1 + 4*x^2;
fa = 5/6*x -1/2*x^2 -1/3*x^4;   %实际解
for ii = 1:n
    fn(ii) = x - x^(ii+1);
end
for ii = 1:n
    for jj = 1:n
        lmn(ii,jj) = -diff(diff(fn(jj))); % 算子作用到 基函数的线性组合上
        lmn(ii,jj) = subs(lmn(ii,jj),x,ii*(1/(n+1))); % 点选配方法:将原 lmn 表达式中连续的 x 用等间隔的离散点 m/(N+1)来表示
    end
    gm(ii) = subs(g,x,ii*(1/(n+1))); % 在源函数上离散化 点选配
end
a = lmn\gm';
f = fn*a;
val = 0:0.01:1; % 函数的定义域
plot(val,subs(f,val),'p')
hold on
plot(val,subs(fa,val),'r')
str=['N = ',num2str(n)];
title(str);
legend('矩量法解','精确值');