function expfitdemo
%   Copyright 2012 Cleve Moler and The MathWorks, Inc.
%   Revised by 张志涌（Zhiyong Zhang）,2022

t = (0:.1:2)';
y = [5.8955 3.5639 2.5173 1.9790 1.8990 1.3938 1.1359 ...
1.0096 1.0343 0.8435 0.6856 0.6100 0.5392 0.3946 ...
0.3903 0.5474 0.3459 0.1370 0.2211 0.1704 0.2636]';
clf
shg
set(gcf,'doublebuffer','on')
h = plot(t,y,'o',t,0*t,'-');
ht=text(1.4,5,' ','HorizontalAlignment','center');
hr=text(1.4,4.5,' ','HorizontalAlignment','center');
axis([0 2 0 6.5])
lambda0 = [3 6]';
Lam10=num2str(lambda0(1),'%.4f');
Lam20=num2str(lambda0(2),'%.4f');
text(1.4,5.5,['\lambda_{10}= ',Lam10,...
    blanks(5),'\lambda_{20}= ',Lam20],...
    'HorizontalAlignment','center');
lambda = fminsearch(@expfitfun,lambda0);
set(h(2),'linewidth',2)
Lam1=num2str(-lambda(1),'%.4f');
Lam2=num2str(-lambda(2),'%.4f');
title(['y \approx',' e^{',Lam1,'} + ','e^{',Lam2,'}'])

function res = expfitfun(lambda)
m = length(t);
n = length(lambda);
X = zeros(m,n);
for j = 1:n
   X(:,j) = exp(-lambda(j)*t);
end
beta = X\y;
z = X*beta;
res = norm(z-y);
set(h(2),'ydata',z);
lam1=num2str(lambda(1),'%.4f');
lam2=num2str(lambda(2),'%.4f');
set(ht,'string',[' \lambda_{1} = ',lam1,...
    blanks(6),'\lambda_{2} = ',lam2])
set(hr,'string',['residuals = ',num2str(res,'%.4f')])
pause(.3)
end
end
