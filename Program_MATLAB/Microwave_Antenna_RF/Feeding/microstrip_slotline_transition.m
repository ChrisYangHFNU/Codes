clear;

Zs = 90; % slotline特性阻抗
d = 1.6; % 介质板厚度，unit：mm
epsilon_r = 4.4; % 介质板相对介电常数
lambda = 30; % 工作波长，unit:mm
lambda_s = lambda*sqrt(2/(epsilon_r+1)); % 槽线工作波长，unit:mm
u = sqrt(epsilon_r-(lambda/lambda_s)^2);
v = sqrt((lambda/lambda_s)^2-1);
q = 2*pi*u*d/lambda+atan(u/v);
N = cos(2*pi*u*d/lambda)-cot(q)*sin(2*pi*u*d/lambda);
Zm = Zs*N^2; % microstrip和slotline交叉点位置输入阻抗
Zm