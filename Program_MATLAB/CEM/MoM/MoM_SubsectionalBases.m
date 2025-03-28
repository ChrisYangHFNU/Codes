clear;clc;
n = input('Input the N...\n');
%subsectional bases
syms x;
g = 1 + 4*x^2;
fa = 5/6*x -1/2*x^2 -1/3*x^4;   %实际解
for ii = 1:n
    %基函数为三角形函数，检验函数为脉冲函数
    xm = ii * 1/(n+1);
    %两种定义基函数的方法同样有用，在这里选择内置piecewise函数
    %fn(i) = (1-abs(x-xm)*(n+1)) .* ( x - xm < (1/(n+1)) &  x - xm > -(1/(n+1)));
    fn(ii) = piecewise(abs(x - xm) < 1/(n+1) ,1-abs(x-xm)*(n+1),0);
end
lmn = ones(n,n);
gm = ones(n,1);
for ii = 1:n
    for jj = 1:n
        if ii == jj
            lmn(ii,jj) = 2*(n+1);
        elseif ii - jj == 1 | jj - ii == 1
            lmn(ii,jj) = - (n+1);
        else
            lmn(ii,jj) = 0;
        end
    end
    gm(ii) = 1/(n+1) * (1+(4*ii^2+1/3)/(n+1)^2);
end
a = lmn\gm;
val = 0:0.01:1;
f = fn * a;
plot(val,eval(subs(f,x,val)));
hold
plot(val,subs(fa,x,val));
str=['N = ',num2str(n)];
title(str);
legend('矩量法解','精确值');