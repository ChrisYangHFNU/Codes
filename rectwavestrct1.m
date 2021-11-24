function rectwavestrct1(a,b,d,H0,f,t)
% ���ƾ��β������ṹ�����м��㵥λΪm������Ϊmm
% f ����Ƶ��
% d ��������
% H0 �������
% t tʱ�̵ĳ��ṹͼ

lightC = 3e8;
a = a/1000;
b = b/1000;
lambdaC = 2*a; % TE10��ֹ����
lambda = lightC/f;
mu = 4*pi*1e-7;
if lambda>lambdaC
    return;
else
    clf;
    lambdaG = lambda/sqrt(1-(lambda/lambdaC)^2);
    c = lambdaG;
    beta = 2*pi/lambdaG;
    omega = 2*pi*f;
    x = 0:a/d:a;
    y = 0:b/d:b;
    z = 0:c/d:c;
    [x1,y1,z1] = meshgrid(x,y,z);
    hx = -beta.*a.*H0.*sin(pi./a.*x1).*sin(omega*t-beta.*z1)./pi;
    hz = H0.*cos(pi./a.*x1).*cos(omega*t-beta.*z1);
    hy = zeros(size(y1));
    quiver3(z1,x1,y1,hz,hx,hy,'b');
    hold on
    x2 = x1-0.001;
    y2 = y1-0.001;
    z2 = z1-0.001;
    ex = zeros(size(x2));
    ey = omega.*mu.*H0.*sin(pi./a.*x2).*sin(omega*t-beta.*z2)./pi;
    ez = zeros(size(z2));
    quiver3(z2,x2,y2,ez,ex,ey,'r');
    xlabel('���䷽��')
    ylabel('�������')
    zlabel('����խ��')
    hold off
end
    