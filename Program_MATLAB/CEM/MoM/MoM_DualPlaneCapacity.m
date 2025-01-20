clear;clc
n = input('Input the N of 2N\n')
syms x y
%计算平行金属板的电容的l矩阵分为ltt、ltb、lbb、lbt，其中t表示top，b表示bottom，为一个2N阶方阵
%ltt = lbb = lmn，lmn与计算单个导体板的lmn相同
%根据对称原理：ltb = lbt
%平板为单位面积平板
num = n^0.5;
v = 1;
a = 1 / 2;
b = a / num;
epsilon = 8.854187817e-12;   %真空中的介电常数

syms x y;         
%lmn是delta_Sn上单位振幅的均匀电荷密度在delta_Sm的中心处产生的电位
for index_d = 1:100
d = index_d;    
for i = 1:n
    for j = 1:n
        if i == j   % m=n的情况
            lmn(i,j) = 2*b/(pi*epsilon) *0.8814;
            ltb(i,j) = 0.282/epsilon*2*b*((1+pi/4*(d/b)^2)^0.5 - pi^0.5*d/(2*b));
        else        % m~=n的情况
            if mod(i,num) == 0
                xm = fix(i/num);
                ym = num;
            else
                xm = fix(i/num) +1;
                ym = mod(i,num);
            end
            if mod(j,num) == 0
                xn = j/num;
                yn = num;
            else
                xn = fix(j/num) + 1;
                yn = mod(j,num);
            end
            deltax = 2*b*(xm-xn);   %Xm与Xn之间的距离
            deltay = 2*b*(ym-yn);   %Ym与Yn之间的距离
            %disp(['i=',num2str(i),'j=',num2str(j),' ',num2str([xm,ym,xn,yn])])
            %验证编号问题，成功
            lmn(i,j) = b^2 / (pi*epsilon* (deltax^2 + deltay^2)^0.5 );
            ltb(i,j) = b^2 /( pi*epsilon* (deltax^2 + deltay^2 +d^2)^0.5 );
        end
    end
end
gm = ones(n,1);
alpha = (lmn-ltb)\gm;
c(index_d) = 4*b*b*sum(sum(inv(lmn-ltb)));
end
plot(c);
str=['N = ',num2str(n)];
title(str);
xlabel('d');
ylabel('C');