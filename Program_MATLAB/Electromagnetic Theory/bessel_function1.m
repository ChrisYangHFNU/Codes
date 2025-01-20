%% 【figure11】：z为实数，三类Bessel函数的输出对比，第三类输出是复数
z = 0:0.1:20;

% 观察\niu的影响
    m =1                            % niu的上界，观察niu的影响用，注意niu的定义下界为0
    J = zeros(m+1,201);          % 创建一个二维输出矩阵
    Y = zeros(m+1,201); 
    H1 = zeros(m+1,201);
    H2 = zeros(m+1,201);
    scale = 1
    for i = 0:m
        J(i+1,:) = besselj(i,z);
        Y(i+1,:) = bessely(i,z);
        H1(i+1,:) = besselh(i,1,z);
        H2(i+1,:) = besselh(i,2,z);
    end

% besselj返回的数值与z同维
    figure(11)
    subplot(1,4,1)
    plot(z,J)
    grid on
    legend('J_0','J_1','J_2','J_3','J_4','Location','Best')
    title('Bessel Functions of the First Kind for $\nu \in [0, 1]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$J_\nu(z)$','interpreter','latex')
    
    subplot(1,4,2)
    plot(z,Y)
    grid on
    legend('Y_0','Y_1','Y_2','Y_3','J_4','Location','Best')
    title('Bessel Functions of the Second Kind for $\nu \in [0, 1]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$Y_\nu(z)$','interpreter','latex')
    
    subplot(1,4,3)
    plot(z,sqrt(real(H1).^2 + imag(H1).^2))
    grid on
    legend('H1_0','H1_1','H1_2','H1_3','H1_4','Location','Best')
    title('Bessel Functions of the Third Kind for $\nu \in [0, 1]$ and $k=1$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$H1_\nu(z)$','interpreter','latex')
    
    subplot(1,4,4)
    plot(z,sqrt(real(H2).^2 + imag(H2).^2))
    grid on
    legend('H2_0','H2_1','H2_2','H2_3','H2_4','Location','Best')
    title('Bessel Functions of the Third Kind for $\nu \in [0, 1]$ and $k=1$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$H2_\nu(z)$','interpreter','latex')