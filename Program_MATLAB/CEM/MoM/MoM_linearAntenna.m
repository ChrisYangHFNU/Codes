clear;clc
%n = input('Input the number of delta C\n');
n = 64
mu = 4*pi*10^-7;      %真空中的磁导率
epsilon = 8.854187817e-12;   %真空中的介电常数
lambda = 1;           %波长
c = 3*10^8;           %真空中光速
w = 2*pi*c/lambda;    %角速度
k = 2*pi / lambda;    %波数
for index_y = 1:250
    l = index_y*lambda*0.01;           %天线总长度

    a = l/148.4;            %线的半径 l/2a = 74.2
    delta_l = l/(n+1);      %因为是线天线，所以分段均匀，并且两端看做带有零电流的一小段

    zmn = zeros(n);     %阻抗矩阵
    v = zeros(n,1);     %电压矩阵
    v(50) = 1;          %在中点加激励电压
    psi = zeros(1:5);
    l_points(:,1) = 0 : l/(2*n) : l;
    
    for index_i = 1:n
        for index_j = 1:n
            psi = cal_psi8(index_i,index_j,a,k,delta_l,l_points);
            zmn(index_i,index_j) = j*w*mu*delta_l*delta_l*psi(1)+1/(j*w*epsilon)*(psi(2) - psi(3) - psi(4) + psi(5));
        end
    end
    ymn = inv(zmn);         
    alpha = ymn*v;
    zz(index_y) = ymn(32,32);
end
xrange = 0.01:0.01:2.5;
plot(xrange,real(zz));
title('Real\_Y')
xlabel('l/lambda')
figure();
plot(xrange,imag(zz));
title('Imag\_Y')
xlabel('l/lambda')

function psi = cal_psi8(m,n,a,k,delta_l,l_points)

        if n == m
                R_mn = sqrt((l_points(2*m-1) - l_points(2*n+1)) * (l_points(2*m-1) - l_points(2*n+1)) + a * a);
                
                psi(1) = 1 / (2 * pi * delta_l) * log(delta_l / a) - 1j * k / (4 * pi); 
                psi(2) = 1 / (2 * pi * delta_l) * log(delta_l / a) - 1j * k / (4 * pi); 
                psi(3) = exp(-1j * k * abs(R_mn)) / (4 * pi * abs(R_mn));
                psi(4) = exp(-1j * k * abs(R_mn)) / (4 * pi * abs(R_mn));
                psi(5) = 1 / (2 * pi * delta_l) * log(delta_l / a) - 1j * k / (4 * pi); 
                
        elseif m - n == 1   
                R_mn(:,1) = sqrt((l_points(2*m) - l_points(2*n)) * (l_points(2*m) - l_points(2*n)) + a * a);
                R_mn(:,2) = sqrt((l_points(2*m+1) - l_points(2*n-1)) * (l_points(2*m+1) - l_points(2*n-1)) + a * a);
                
                psi(1) = exp(-1j * k * abs(R_mn(1))) / (4 * pi * abs(R_mn(1))); 
                psi(2) = exp(-1j * k * abs(R_mn(1))) / (4 * pi * abs(R_mn(1)));
                psi(3) = exp(-1j * k * abs(R_mn(2))) / (4 * pi * abs(R_mn(2)));
                psi(4) = 1 / (2 * pi * delta_l) * log(delta_l / a) - 1j * k / (4 * pi);
                psi(5) = exp(-1j * k * abs(R_mn(1))) / (4 * pi * abs(R_mn(1)));
        
        elseif m - n == -1
                R_mn(:,1) = sqrt((l_points(2*m) - l_points(2*n)) * (l_points(2*m) - l_points(2*n)) + a * a);
                R_mn(:,2) = sqrt((l_points(2*m-1) - l_points(2*n+1)) * (l_points(2*m-1) - l_points(2*n+1)) + a * a);
                
                psi(1) = exp(-1j * k * abs(R_mn(1))) / (4 * pi * abs(R_mn(1)));
                psi(2) = exp(-1j * k * abs(R_mn(1))) / (4 * pi * abs(R_mn(1)));
                psi(3) = 1 / (2 * pi * delta_l) * log(delta_l / a) - 1j * k / (4 * pi);
                psi(4) = exp(-1j * k * abs(R_mn(2))) / (4 * pi * abs(R_mn(2)));
                psi(5) = exp(-1j * k * abs(R_mn(1))) / (4 * pi * abs(R_mn(1)));
        else
                R_mn(:,1) = sqrt((l_points(2*m) - l_points(2*n)) * (l_points(2*m) - l_points(2*n)) + a * a);
                R_mn(:,2) = sqrt((l_points(2*m+1) - l_points(2*n-1)) * (l_points(2*m+1) - l_points(2*n-1)) + a * a);
                R_mn(:,3) = sqrt((l_points(2*m-1) - l_points(2*n+1)) * (l_points(2*m-1) - l_points(2*n+1)) + a * a);
                
                psi(1) = exp(-1j * k * abs(R_mn(1))) / (4 * pi * abs(R_mn(1))); 
                psi(2) = exp(-1j * k * abs(R_mn(1))) / (4 * pi * abs(R_mn(1))); 
                psi(3) = exp(-1j * k * abs(R_mn(2))) / (4 * pi * abs(R_mn(2)));
                psi(4) = exp(-1j * k * abs(R_mn(3))) / (4 * pi * abs(R_mn(3)));
                psi(5) = exp(-1j * k * abs(R_mn(1))) / (4 * pi * abs(R_mn(1)));
                
        end
end