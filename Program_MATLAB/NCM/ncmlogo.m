function ncmlogo
% NCMLOGO
% L-shaped membrane on the cover of the print
% version of Numerical Computing with MATLAB.
% Activate the cameratoolbar and move the camera.

%   Copyright 2013 Cleve Moler and The MathWorks, Inc.
%   Revised by 张志涌（Zhiyong Zhang）,2022


set(gcf,'numbertitle','off','name','ncmlogo', ...
  'colormap',jet(8),'color',[0 0 1/4])

axes('pos',[0.1,0.15,0.8,0.8])%[0,0,1,1]
axis off
daspect([1 1 1])
cc=jet(8);

% Compute MathWorks logo
L = rot90(membranetx(1,32,10,10),2);
% Filled contour plot with transparent lifted patches
b = (1/16:1/8:15/16)';
hold on
for k = 1:8
   c = contourf(L,[b(k) b(k)],'fill','off','LineStyle','none');
   m=length(c)-1;
   x=c(1,2:end);
   y=c(2,2:end);
   z=4*k*ones(size(x));
   hf(k)=fill3(x,y,z,cc(k,:),'EdgeColor','w',...
       'LineWidth',2,'FaceAlpha',.5);
end
hold off
view(12,30)
pause(5)
axis([0 75 0 75 0 40])
NN=12;kk=0;an=30;
while kk<NN
kk=kk+1;    
rotate(hf(1),[0,0,1],an)
rotate(hf(2),[0,0,1],an)
rotate(hf(3),[0,0,1],an)
rotate(hf(4),[0,0,1],an)
rotate(hf(5),[0,0,1],an)
rotate(hf(6),[0,0,1],an)
rotate(hf(7),[0,0,1],an)
rotate(hf(8),[0,0,1],an)
axis image
pause(0.15)
end
