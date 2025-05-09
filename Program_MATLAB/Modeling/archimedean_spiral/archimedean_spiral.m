% Archimedean Spiral
clc;clear;close all;format shortg;

turns = 2.65; % no. of turns of the spiral
w = 0.3; % width of the spiral arms
b1 = 0.3; % initial distance of the spiral from the center
a1 = (2*w)/pi; % spiral factor, determines how fast the spiral grows

ts1 = 0;
te1 = 2*pi*turns;
np = 500;
t1 = linspace(ts1,te1,np);
r1 = b1+w+a1*te1;
te2 = (r1-b1)/a1;
t2 = linspace(ts1,te2,np);

x1 = (b1+a1*t2).*cos(t2);
y1 = (b1+a1*t2).*sin(t2);

x2 = -(b1+a1*t2).*cos(t2);
y2 = -(b1+a1*t2).*sin(t2);

x3 = (b1+w+a1*t1).*cos(t1);
y3 = (b1+w+a1*t1).*sin(t1);

x4 = -(b1+w+a1*t1).*cos(t1);
y4 = -(b1+w+a1*t1).*sin(t1);

t3 = linspace(te1,te2,np);
x5 = r1*cos(t3);
y5 = r1*sin(t3);

x6 = -x5;
y6 = -y5;

figure;plot(x1,y1,'-b');
hold on;plot(x2,y2,'-r')
hold on;plot(x3,y3,'-b');
hold on;plot(x4,y4,'-r');
hold on;plot(x5,y5,'-b');
hold on;plot(x6,y6,'-r');
hold on;plot([b1 b1+w],[0 0],'-b');
hold on;plot(-[b1 b1+w],[0 0],'-r');

n = 4;
xlim([-n n]); ylim([-n,n]);