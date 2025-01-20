%% 【figure16】：z为复数时，缩放对三类Bessel函数输出的影响

% 定义域和scale系数
    x = -10:0.3:10;
    y = x';
    z = x + 1i*y;
    scale = 1;
    nu = 2;

%第一类
    J = besselj(nu,z);
    Js = besselj(nu,z,scale);

%第二类
    Y = bessely(nu,z);
    Ys = bessely(nu,z,scale);

%第三类k=1
    H1 = besselh(nu,1,z);
    H1s = besselh(nu,1,z,scale);

%第三类k=2
    H2 = besselh(nu,2,z);
    H2s = besselh(nu,2,z,scale);

% 虚部对比
    figure(16)
    subplot(4,2,1)
    surf(x,y,imag(J))
    title('Bessel Function of the First Kind, imag_J','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')

    subplot(4,2,2)
    surf(x,y,imag(Js))
    title('Scaled Bessel Function of the First Kind, imag_Js','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')

    subplot(4,2,3)
    surf(x,y,imag(Y))
    title('Bessel Function of the Second Kind, imag_Y','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')

    subplot(4,2,4)
    surf(x,y,imag(Ys))
    title('Scaled Bessel Function of the Second Kind, imag_Ys','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')

    subplot(4,2,5)
    surf(x,y,imag(H1))
    title('Bessel Function of the Third Kind, imag_H1','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')

    subplot(4,2,6)
    surf(x,y,imag(H1s))
    title('Scaled Bessel Function of the Third Kind, imag_H1s','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')

    subplot(4,2,7)
    surf(x,y,imag(H2))
    title('Bessel Function of the Third Kind, imag_H2','interpreter','latex')
    xlabel('imag(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')

    subplot(4,2,8)
    surf(x,y,imag(H2s))
    title('Scaled Bessel Function of the Third Kind, imag_H2s','interpreter','latex')
    xlabel('imag(z)','interpreter','latex')
    ylabel('imag(z)','interpreter','latex')

% 实部对比
    figure(17)
    subplot(4,2,1)
    surf(x,y,real(J))
    title('Bessel Function of the First Kind, real_J','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')

    subplot(4,2,2)
    surf(x,y,real(Js))
    title('Scaled Bessel Function of the First Kind, real_Js','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')

    subplot(4,2,3)
    surf(x,y,real(Y))
    title('Bessel Function of the Second Kind, real_Y','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')

    subplot(4,2,4)
    surf(x,y,real(Ys))
    title('Scaled Bessel Function of the Second Kind, real_Ys','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')

    subplot(4,2,5)
    surf(x,y,real(H1))
    title('Bessel Function of the Third Kind, real_H1','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')

    subplot(4,2,6)
    surf(x,y,real(H1s))
    title('Scaled Bessel Function of the Third Kind, real_H1s','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')

    subplot(4,2,7)
    surf(x,y,real(H2))
    title('Bessel Function of the Third Kind, real_H2','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')

    subplot(4,2,8)
    surf(x,y,real(H2s))
    title('Scaled Bessel Function of the Third Kind, real_H2s','interpreter','latex')
    xlabel('real(z)','interpreter','latex')
    ylabel('real(z)','interpreter','latex')