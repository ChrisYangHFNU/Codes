function [x_data,y_data]=Smith_plotGammaCircle(ax,Z,Z0)
% plot SWR circle, lossless TL

gamma=z2gamma(Z,Z0);  
r=abs(gamma); 
alpha=0:2*pi/100:2*pi;                                       
hold all;
sub_hp2=plot(ax,r*cos(alpha),r*sin(alpha),'-','LineWidth',.5,'Color',[1 .2 0])   
x_data=sub_hp2.XData
y_data=sub_hp2.YData

end