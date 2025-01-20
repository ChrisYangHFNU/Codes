close all;clear;clc;

N = 10;
dr = pi/2;

Ntheta = 108;
Nphi = 72;
dtheta = pi/Ntheta;
dphi = 2*pi/Nphi;
theta = linspace(0,pi,Ntheta);
phi = 0:dphi:2*pi-dphi;
theta([1 end]) = [eps pi-eps];
[THETA,PHI] = meshgrid(theta,phi);

Theta = dr*(cos(theta)-1)-pi/N;
E = sin(pi/2/N)*sin(N*Theta/2)./sin(Theta/2);
En = E/max(abs(E));
D = 1/sum(En.^2.*sin(theta)*dtheta/2);
En = repmat(En,Nphi,1);

figure(1),hold on
set(gcf,'render','painter','color','w')
view(0,8),axis image,colormap jet

[x,y,z] = sph2cart([PHI;PHI(1,:)],pi/2-[THETA;THETA(1,:)],abs([En;En(1,:)]));
surf(x,y,z,abs(En),'facea',.8,'edgea',.5);