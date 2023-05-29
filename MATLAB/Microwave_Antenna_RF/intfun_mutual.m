% 边靠边平行对称振子的互阻抗，利用互阻抗的原始积分公式计算 

clear;clc;
global k l dp;
j = sqrt(-1);
C = 0.5772; % Euler constant
c = 3.0e8;
f = 3.0e8;
k = 2*pi*f/c;
l = 0.5;
d = [0:.01:3.0];
len = length(d);
Zr = zeros(1,len);
for num=1:len
    dp = d(1,num);
    Zr(1,num) = j*30*(integral(@myfun01,-l/2,0)+integral(@myfun02,0,l/2));
end
figure(1)
plot(d,real(Zr),'r.-');hold on;
plot(d,imag(Zr),'b.-');hold on;
grid on;
axis([0 3.0 -40 80]);

function y = myfun01(z)
    global k l dp;
    j = sqrt(-1);
    R1 = sqrt(dp^2+(z-l/2).^2);
    r = sqrt(dp^2+z.^2);
    R2 = sqrt(dp^2+(z+l/2).^2);
    y = sin(k*(l/2+z)).*(exp(-j*k*R1)./R1+exp(-j*k*R2)./R2-2*cos(k*l/2)*exp(-j*k*r)./r);
end

function y = myfun02(z)
    global k l dp;
    j = sqrt(-1);
    R1 = sqrt(dp^2+(z-l/2).^2);
    r = sqrt(dp^2+z.^2);
    R2 = sqrt(dp^2+(z+l/2).^2);
    y = sin(k*(l/2-z)).*(exp(-j*k*R1)./R1+exp(-j*k*R2)./R2-2*cos(k*l/2)*exp(-j*k*r)./r);
end