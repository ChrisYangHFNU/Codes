clear;
clc;
n = input('Input the N \n');
% Galerkin Method
syms x 
g = 1+4*x^2;
fa = 5/6*x-1/2*x^2-1/3*x^4; %ʵ�ʽ�

for ii=1:n
    fn(ii) = x-x^(ii+1);
    wn(ii) = fn(ii);
    lfn(ii) = -diff(diff(fn(ii)));
end
for ii=1:n
    for jj=1:n
        lmn(ii,jj) = int(wn(ii)*lfn(jj),x,0,1); %���ڻ�
    end
    gm(ii) = int(wn(ii)*g,x,0,1);
end

a = lmn\gm';
f = fn*a; % ��Ϊ������״������Ӧ����˷�����˳��

val = 0:.01:1;
plot(val,subs(f,val),'p')
hold on
plot(val,subs(fa,val),'r')
str = ['N = ',num2str(n)];
title(str);
legend('��������','��ȷֵ')