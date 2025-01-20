%% 【figure14】：z为实数，scale对第一、二类Bessel修正函数输出影响的对比
    z = 0:0.1:20;

% 观察\niu的影响
    m = 2                               % niu的上界，观察niu的影响用，注意niu的定义下界为0
    Is = zeros(m+1,201);          % 创建一个二维输出矩阵
    Ks = zeros(m+1,201); 
    scale = 1
    for i = 0:m
        Is(i+1,:) = besseli(i,z,scale);
        Ks(i+1,:) = besselk(i,z,scale);
    end

    I = zeros(m+1,201);          % 创建一个二维输出矩阵
    K = zeros(m+1,201); 
    scale = 1
    for i = 0:m
        I(i+1,:) = besseli(i,z);
        K(i+1,:) = besselk(i,z);
    end

% besselj返回的数值与z同维
    figure(14)
    subplot(2,2,1)
    plot(z,Is)
    grid on
    legend('Is_0','Is_1','Is_2','Is_3','Is_4','Location','Best')
    title('Scaled Modified Bessel Functions of the First Kind for $\nu \in [0, 2]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$Is_\nu(z)$','interpreter','latex')
    
    figure(14)
    subplot(2,2,2)
    plot(z,Ks)
    grid on
    legend('Ks_0','Ks_1','Ks_2','Ks_3','Ks_4','Location','Best')
    title('Scaled Modified Bessel Functions of the Second Kind for $\nu \in [0, 2]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$Ks_\nu(z)$','interpreter','latex')
    
    figure(14)
    subplot(2,2,3)
    plot(z,I)    
    legend('I_0','I_1','I_2','I_3','I_4','Location','Best')
    title('Modified Bessel Functions of the First Kind for $\nu \in [0, 2]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$I_\nu(z)$','interpreter','latex')
    
    figure(14)
    subplot(2,2,4)
    plot(z,K)
    grid on
    legend('K_0','K_1','K_2','K_3','K_4','Location','Best')
    title('Modified Bessel Functions of the Second Kind for $\nu \in [0, 2]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$K_\nu(z)$','interpreter','latex')