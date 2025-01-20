%% 【figure12】：z为实数，两类Bessel修正函数的输出对比
z = 0:0.1:20;

% 观察\niu的影响
    m =1                            % niu的上界，观察niu的影响用，注意niu的定义下界为0
    I= zeros(m+1,201);          % 创建一个二维输出矩阵
    K = zeros(m+1,201); 
    scale = 1
    for i = 0:m
        I(i+1,:) = besseli(i,z);
        K(i+1,:) = besselk(i,z);
    end

% besselj返回的数值与z同维
    figure(12)
    subplot(3,1,1)
    plot(z,I)
    grid on
    legend('I_0','I_1','I_2','I_3','I_4','Location','Best')
    title('Bessel Functions of the First Kind for $\nu \in [0, 1]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$I_\nu(z)$','interpreter','latex')
    
    figure(12)
    subplot(2,1,1)
    plot(z,I)
    grid on
    legend('Y_0','Y_1','Y_2','Y_3','J_4','Location','Best')
    title('Bessel Functions of the Second Kind for $\nu \in [0, 1]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$Y_\nu(z)$','interpreter','latex')
    
    figure(12)
    subplot(2,1,2)
    plot(z,K)
    grid on
    legend('Y_0','Y_1','Y_2','Y_3','J_4','Location','Best')
    title('Bessel Functions of the Second Kind for $\nu \in [0, 1]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$Y_\nu(z)$','interpreter','latex')