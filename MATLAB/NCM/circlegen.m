function circlegen(h)
%CIRCLEGEN  Generate approximate circles.
%   CIRCLEGEN(h) uses step size h.
%   CIRCLE with no arguments uses h = 0.20906.

%   Copyright 2013 Cleve Moler and The MathWorks, Inc.
%   Revised by 张志涌（Zhiyong Zhang）,2022

if nargin < 1
   h = sqrt(2*(1-cos(2*pi/30)));
end

shg
clf
set(gcf,'menu','none','numbertitle','off',...
    'name','Circlegen','Color','w')
z = exp(2*pi*1i*(0:256)/256);
grey = [.8 .8 .8];
line(real(z),imag(z),'Color',grey)
line(0,0,'Marker','o','Color',grey)
p=animatedline(1,0,'Marker','.','linestyle','none','Color','k');
axis([-2 2 -2 2])
axis square,box on
title(['h = ' num2str(h)])
stop = uicontrol('style','toggle','string','stop');

x = 1;
y = 0;
while ~get(stop,'value')
   x = x + h*y;
   y = y - h*x;
   addpoints(p,x,y)
   drawnow
end

set(stop,'string','close','value',0,'callback','close')
