% 采用切比雪夫综合法 给出低副瓣阵列天线的幅相值
% author: yangjing 2024/11/09
% x = cos(u)
% T_m(x) = cos(mu) = cos(m*arccosx)
% T_0(x) = cos(0); T_1(x) = cos(u) = x
% T_m+1(x) = 2x*T_m(x)-T_m-1(x)
% 

clear 
close all

% constant parameters
eleNum = 16;
eleSpa = 0.6;
sll_dB = 30; % 注意符号
sll = 10^(sll_dB/20);

% 调用 chebyshevpolynomial.m 计算切比雪夫多项式系数
% S_even(u) =
% I1*cos(u)+I2*cos(3u)+I3*cos(5u)+I4*cos(7u)+I5*cos(9u)+I6*cos(11u)+I7*cos(13u)+I8*cos(15u)
% T_1(x) = cos(u) = x
% T_3(x) = cos(3u) = 4x^3-3x
% T_5(x) = cos(5u) = 16x^5-20x^3+5x
% T_7(x) = cos(7u) = 64x^7-112x^5+56x^3-7x
% T_9(x) = cos(9u) = 256x^9-576x^7+432x^5-120x^3+9x
% T_11(x) = cos(11u) = 1024x^11-2816x^9+2816x^7-1232x^5+220x^3-11x
% T_13(x) = cos(13u) =
% 4096x^13-13312x^11+16640x^9-9984x^7+2912x^5-364x^3+13x
% T_15(x) = cos(15u) =
% 16384x^15-61140x^13+92160x^11-70400x^9+28800x^7-6048x^5+560x^3-15x

x0 = 0.5*((sll+sqrt(sll^2-1))^(1/(eleNum-1))+(sll-sqrt(sll^2-1))^(1/(eleNum-1)));

I1_coef = [1 0];
I2_coef = [4 0 -3 0];
I3_coef = [16 0 -20 0 5 0];
I4_coef = [64 0 -112 0 56 0 -7 0];
I5_coef = [256 0 -576 0 432 0 -120 0 9 0];
I6_coef = [1024 0 -2816 0 2816 0 -1232 0 220 0 -11 0];
I7_coef = [4096 0 -13312 0 16640 0 -9984 0 2912 0 -364 0 13 0];
I8_coef = [16384 0 -61140 0 92160 0 -70400 0 28800 0 -6048 0 560 0 -15 0];

I_coef = I5_coef;

% I8 = I_coef(1)*(x0^(eleNum-1))/I8_coef(1);
% I7 = (I_coef(3)*(x0^(eleNum-3))-I8_coef(3)*I8)/I7_coef(1);
% I6 = (I_coef(5)*(x0^(eleNum-5))-I8_coef(5)*I8-I7_coef(3)*I7)/I6_coef(1);
% I5 = (I_coef(7)*(x0^(eleNum-7))-I8_coef(7)*I8-I7_coef(5)*I7-I6_coef(3)*I6)/I5_coef(1);
% I4 = (I_coef(9)*(x0^(eleNum-9))-I8_coef(9)*I8-I7_coef(7)*I7-I6_coef(5)*I6-I5_coef(3)*I5)/I4_coef(1);
% I3 = (I_coef(11)*(x0^(eleNum-11))-I8_coef(11)*I8-I7_coef(9)*I7-I6_coef(7)*I6-I5_coef(5)*I5-I4_coef(3)*I4)/I3_coef(1);
% I2 = (I_coef(13)*(x0^(eleNum-13))-I8_coef(13)*I8-I7_coef(11)*I7-I6_coef(9)*I6-I5_coef(7)*I5-I4_coef(5)*I4-I3_coef(3)*I3)/I2_coef(1);
% I1 = (I_coef(15)*(x0^(eleNum-15))-I8_coef(15)*I8-I7_coef(13)*I7-I6_coef(11)*I6-I5_coef(9)*I5-I4_coef(7)*I4-I3_coef(5)*I3-I2_coef(3)*I2)/I1_coef(1);
% I = [I8 I7 I6 I5 I4 I3 I2 I1 I1 I2 I3 I4 I5 I6 I7 I8];

% I5 = I_coef(1)*x0^(eleNum-1)/I5_coef(1);
% I4 = (I_coef(3)*x0^(eleNum-3)-I5_coef(3)*I5)/I4_coef(1);
% I3 = (I_coef(5)*x0^(eleNum-5)-I5_coef(5)*I5-I4_coef(3)*I4)/I3_coef(1);
% I2 = (I_coef(7)*x0^(eleNum-7)-I5_coef(7)*I5-I4_coef(5)*I4-I3_coef(3)*I3)/I2_coef(1);
% I1 = (I_coef(9)*x0^(eleNum-9)-I5_coef(9)*I5-I4_coef(7)*I4-I3_coef(5)*I3-I2_coef(3)*I2)/I1_coef(1);
% I = [I5 I4 I3 I2 I1 I1 I2 I3 I4 I5];

% =========== Barbiere Method ===========
I_right = zeros(1,eleNum/2); % 偶数阵列
for ii=1:eleNum/2
    for jj=ii:eleNum/2
        I_right(ii) = I_right(ii)+(-1)^(eleNum/2-jj)*x0^(2*jj-1)*factorial(jj+eleNum/2-2)*(2*eleNum/2-1)...
            /(factorial(jj-ii)*factorial(jj+ii-1)*factorial(eleNum/2-jj));
    end
end
% =======================================

I_left = fliplr(I_right);
I = [I_left I_right];

I = I/max(I);
x = [angle(I) abs(I)];

fm = linearArrayfactor(eleNum,-90:0.5:90,eleSpa,x);