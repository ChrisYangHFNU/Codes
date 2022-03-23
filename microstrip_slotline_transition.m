clear;

Zs = 90; % slotline�����迹
d = 1.6; % ���ʰ��ȣ�unit��mm
epsilon_r = 4.4; % ���ʰ���Խ�糣��
lambda = 30; % ����������unit:mm
lambda_s = lambda*sqrt(2/(epsilon_r+1)); % ���߹���������unit:mm
u = sqrt(epsilon_r-(lambda/lambda_s)^2);
v = sqrt((lambda/lambda_s)^2-1);
q = 2*pi*u*d/lambda+atan(u/v);
N = cos(2*pi*u*d/lambda)-cot(q)*sin(2*pi*u*d/lambda);
Zm = Zs*N^2; % microstrip��slotline�����λ�������迹
Zm