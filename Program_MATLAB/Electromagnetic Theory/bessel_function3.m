%% 【figure13】：z为实数，scale对第一、二类Bessel函数输出影响的对比
    z = 0:0.1:20;

% 观察\niu的影响
    m = 2                               % niu的上界，观察niu的影响用，注意niu的定义下界为0
    Js = zeros(m+1,201);          % 创建一个二维输出矩阵
    Ys = zeros(m+1,201); 
    scale = 1
    for i = 0:m
        Js(i+1,:) = besselj(i,z,scale);
        Ys(i+1,:) = bessely(i,z,scale);
    end

    J = zeros(m+1,201);          % 创建一个二维输出矩阵
    Y = zeros(m+1,201); 
    scale = 1
    for i = 0:m
        J(i+1,:) = besselj(i,z);
        Y(i+1,:) = bessely(i,z);
    end

% besselj返回的数值与z同维
    figure(13)
    subplot(2,2,1)
    plot(z,Js)
    grid on
    legend('Js_0','Js_1','Js_2','Js_3','Js_4','Location','Best')
    title('Scaled Bessel Functions of the First Kind for $\nu \in [0, 2]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$Js_\nu(z)$','interpreter','latex')
    
    figure(13)
    subplot(2,2,2)
    plot(z,Ys)
    grid on
    legend('Ys_0','Ys_1','Ys_2','Ys_3','Ys_4','Location','Best')
    title('Scaled Bessel Functions of the Second Kind for $\nu \in [0, 2]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$Ys_\nu(z)$','interpreter','latex')
    
    figure(13)
    subplot(2,2,3)
    plot(z,J)    
    legend('J_0','J_1','J_2','J_3','J_4','Location','Best')
    title('Bessel Functions of the First Kind for $\nu \in [0, 2]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$J_\nu(z)$','interpreter','latex')
    
    figure(13)
    subplot(2,2,4)
    plot(z,Y)
    grid on
    legend('Y_0','Y_1','Y_2','Y_3','J_4','Location','Best')
    title('Bessel Functions of the Second Kind for $\nu \in [0, 2]$','interpreter','latex')
    xlabel('z','interpreter','latex')
    ylabel('$Y_\nu(z)$','interpreter','latex')