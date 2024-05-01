function [tfinal,yfinal]=orbit(reltol)
%orbit  the two-body problem
%   Revised by 张志涌（Zhiyong Zhang）,2022

y0=[1;0;0;0.3];
opts=odeset('events',@gstop,'reltol',reltol);
[t,y,te,ye]=ode45(@twobody,[0 2*pi],y0,opts,y0);
format long
tfinal=te(end);
yfinal=ye(end,1:2);
plot(y(:,1),y(:,2),'-',0,0,'ro')
axis([-.1 1.05 -.35 .35])
title(['相对容差为 ',num2str(reltol,'%.1e'),' 时的计算轨道'])

function ydot = twobody(t,y,y0)
r = sqrt(y(1)^2 + y(2)^2);
ydot = [y(3); y(4); -y(1)/r^3; -y(2)/r^3];
end
function [val,isterm,dir] = gstop(t,y,y0)
d = y(1:2)-y0(1:2);
v = y(3:4);
val = d'*v;
isterm = 1;
dir = 1;
end
end
