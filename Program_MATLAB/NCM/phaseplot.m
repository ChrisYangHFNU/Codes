function flag = phaseplot(t,y,job)
%   Copyright 2012 Cleve Moler and The MathWorks, Inc.
%   Revised by 张志涌（Zhiyong Zhang）,2022

persistent p
if isequal(job,'init')
p=animatedline(y(1),y(2),'Marker','o');
axis([-1.2 1.2 -1.2 1.2])
axis square
flag = 0;
elseif isequal(job,'')
addpoints(p,y(1),y(2))
pause(0.2)
flag = 0;
end
