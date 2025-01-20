%% 【figure18】【figure19】：z为复数，缩放对修正Bessel函数的的影响

% 定义域和scale系数
    x = -10:0.3:10;
    y = x';
    z = x + 1i*y;
    scale = 1;
    nu = 2

%第一类
    I = besseli(nu,z);
    Is = besseli(nu,z,scale);

%第二类
    K = besselk(nu,z);
    Ks = besselk(nu,z,scale);

% 虚部对比
    figure(18)
    subplot(2,2,1)
    surf(x,y,imag(I))
    title('Bessel Function of the First Kind, imag_I','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')

    subplot(2,2,2)
    surf(x,y,imag(Is))
    title('Scaled Bessel Function of the First Kind, imag_Is','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')

    subplot(2,2,3)
    surf(x,y,imag(K))
    title('Bessel Function of the Second Kind, imag_K','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')

    subplot(2,2,4)
    surf(x,y,imag(Ks))
    title('Scaled Bessel Function of the Second Kind, imag_Ks','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')


% 实部对比
    figure(19)
    subplot(2,2,1)
    surf(x,y,real(I))
    title('Bessel Function of the First Kind, real_I','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')

    subplot(2,2,2)
    surf(x,y,real(Is))
    title('Scaled Bessel Function of the First Kind, real_Is','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')

    subplot(2,2,3)
    surf(x,y,real(K))
    title('Bessel Function of the Second Kind, real_K','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')

    subplot(2,2,4)
    surf(x,y,real(Ks))
    title('Scaled Bessel Function of the Second Kind, real_Ks','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')
