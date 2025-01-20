%% 【figure15】：z为复数时三类Bessel函数的区别

% Hankel函数
    [X,Y] = meshgrid(-4:0.002:2,-1.5:0.002:1.5);
    
    H = besselh(0,X+1i*Y);
    
    figure(15)
    subplot(3,2,1)
    contour(X,Y,abs(H),0:0.2:3.2)
    legend('H_0','Location','Best')
    title('Bessel Functions of the Third Kind of H_amplitude for $\nu = 0$','interpreter','latex')
    xlabel('X','interpreter','latex')
    ylabel('Y','interpreter','latex')   
    
    subplot(3,2,2)
    contour(X,Y,rad2deg(angle(H)),-180:10:180)
    legend('H_0','Location','Best')
    title('Bessel Functions of the Third Kind of H_phase for $\nu = 0$','interpreter','latex')
    xlabel('X','interpreter','latex')
    ylabel('Y','interpreter','latex')
    
% 第一类Bessel函数    
    J = besselj(0,X+1i*Y);
    
    figure(15)
    subplot(3,2,3)
    contour(X,Y,abs(J),0:0.2:3.2)
    legend('J_0','Location','Best')
    title('Bessel Functions of the First Kind of J_amplitude for $\nu = 0$','interpreter','latex')
    xlabel('X','interpreter','latex')
    ylabel('Y','interpreter','latex')   
    
    subplot(3,2,4)
    contour(X,Y,rad2deg(angle(J)),-180:10:180)
    legend('J_0','Location','Best')
    title('Bessel Functions of the First Kind of J_phase for $\nu = 0$','interpreter','latex')
    xlabel('X','interpreter','latex')
    ylabel('Y','interpreter','latex')
    
% 第二类Bessel函数    
    YB= bessely(0,X+1i*Y);
    
    figure(15)
    subplot(3,2,5)
    contour(X,Y,abs(YB),0:0.2:3.2)
    legend('YB_0','Location','Best')
    title('Bessel Functions of the Second Kind of YB_amplitude for $\nu = 0$','interpreter','latex')
    xlabel('X','interpreter','latex')
    ylabel('Y','interpreter','latex')   
    
    subplot(3,2,6)
    contour(X,Y,rad2deg(angle(YB)),-180:10:180)
    legend('YB_0','Location','Best')
    title('Bessel Functions of the Second Kind of YB_phase for $\nu = 0$','interpreter','latex')
    xlabel('X','interpreter','latex')
    ylabel('Y','interpreter','latex')  