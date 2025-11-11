function fm = linearArrayfactor(eleNum,thetaRange,eleSpace,x)
%
% fm = linearArrayfactor(eleNum,thetaRange,eleSpace,x)
% eleNum:阵元数
% thetaRange:方向图观察角度范围
% eleSpace:阵元间距，归一化波长
% x:阵元上激励电流的幅度和相位，相位在前rad，幅度在后线性
%{
4.8GHz 16阵元线阵
pha1 = [8.23 -88.56 172.1 83.03 -3.11 -93.2 178.3 81.17 -10.37 -98.64 -178.19 93.72 -10.89 -86.26 -171.88 80.49];
mag1 = [0.202 0.197 0.175 0.184 0.179 0.214 0.246 0.25 0.262 0.247 0.228 0.266 0.214 0.188 0.178 0.214];
x1 = [deg2rad(pha1) mag1];
4.32GHz 16阵元线阵
pha2 = [14.93 -59.47 -152.49 122.76 51.08 -58.74 -118.02 161.44 69.66 -10.99 -84.2 -166.68 105.09 28.6 -51.03 -145.89];
mag2 = [0.142 0.227 0.205 0.276 0.199 0.175 0.221 0.255 0.263 0.213 0.165 0.233 0.214 0.176 0.189 0.157];
x2 = [deg2rad(pha2) mag2];
%}


fm = zeros(1,length(thetaRange));

for ii=1:eleNum
    fm = fm+x(ii+eleNum)*exp(1j*(ii-1)*2*pi*eleSpace*sind(thetaRange)+1j*x(ii));
end
 fm = abs(fm)./max(abs(fm)); %归一化
 %figure
 plot(thetaRange,20*log10(fm),'k-.','LineWidth',2) % decibal value
 %plot(thetaRange,fm,'b-')
hold on

%{
D_theta = zeros(1,length(thetaRange));
for ii=1:length(thetaRange)
    if thetaRange(ii)<=-2
        D_theta(ii) = 0.001;
    elseif thetaRange(ii)<0
        D_theta(ii) = 0.4*thetaRange(ii)+0.8;
    elseif thetaRange(ii)<4
        D_theta(ii) = 0.05*thetaRange(ii)+0.8;
    elseif thetaRange(ii)<8
        D_theta(ii) = 1;
    elseif thetaRange(ii)<58
        D_theta(ii) = -0.01*thetaRange(ii)+1.08;
    else
        D_theta(ii) = 0.001;
    end
end
 plot(thetaRange,20*log10(D_theta),'r-')
 ylim([-50 0])
 plot(thetaRange,D_theta,'r-')

end
%}

% optimiztionX = [16.2333 172.3717 -117.6549 -162.0345 -24.2158 -37.9609 -102.4292 -99.0187...
%     -6.3424 -88.8941 68.3311 70.0410 131.9329 74.5378 26.5732 76.2868...
%     0.3696 0.7798 0.8206 0.6955 1.0000 1.0000 0.5534 0.4261...
%     0.6564 0.3807 0.0455 0.1865 0.3125 0.1814 0.1559 0.1640];