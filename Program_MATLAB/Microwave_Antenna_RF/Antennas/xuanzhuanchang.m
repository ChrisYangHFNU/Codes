%% 计算两个电基本振子构成的旋转场天线的方向图，并动画演示其工作过程
clear all;clc;
% 计算旋转场天线的立体方向图
theta = meshgrid(eps:pi/180:pi);
phi = meshgrid(eps:2*pi/180:2*pi)';
f = sqrt(1+cos(theta).^2);
fmax = max(max(f));
[x,y,z] = sph2cart(phi,pi/2-theta,f/fmax);
figure(1);
mesh(x,y,z);title('旋转场天线的立体方向图');
axis([-1 1 -1 1 -1 1])

% 计算旋转场天线的E面和H面方向图
theta = pi/2;
phi = linspace(eps,2*pi,100);
fE = sqrt(1+cos(theta)^2)*ones(1,100);
fEmax = max(max(fE));
figure(2)
subplot(1,2,1)
polar(phi,fE/fEmax);title('E面方向图')
theta = linspace(0,2*pi,100);
fH = sqrt(1+cos(theta).^2);
fHmax = max(fH);
subplot(1,2,2);
polar(theta-pi/2,fH/fHmax);title('H面方向图')

% 旋转场天线的动画演示
omega = pi; % 变化周期
%m = moviein(100); % moviein is no longer needed as of MATLAB Release 11 (5.3)
for ii=1:100
    t = 0.1*ii; 
    F_t = abs(sin(omega*t+phi));
    figure(3);
    polar(phi,F_t,'r');title('旋转场动画演示');
    m(ii) = getframe(gcf);
end
%movie(m,1)

% 将动画写入.avi文件
avi1 = VideoWriter('xuanzhuanchang.avi');
avi1.FrameRate = 5;
open(avi1);
writeVideo(avi1,m);
close(avi1)