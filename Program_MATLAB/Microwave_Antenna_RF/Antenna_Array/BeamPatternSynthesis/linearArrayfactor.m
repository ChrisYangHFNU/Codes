function fm = linearArrayfactor(eleNum,thetaRange,eleSpace,x)
%
% fm = linearArrayfactor(eleNum,thetaRange,eleSpace,x)
% eleNum:阵元数
% thetaRange:方向图观察角度范围
% eleSpace:阵元间距
% x:阵元上激励电流的幅度和相位，相位在前，幅度在后

fm = zeros(1,length(thetaRange));

for ii=1:eleNum
    fm = fm+x(ii+eleNum)*exp(1j*(ii-1)*2*pi*eleSpace*sind(thetaRange)+1j*x(ii));
end
 fm = abs(fm)./max(abs(fm)); %归一化
 figure
 plot(thetaRange,20*log10(fm),'k','LineWidth',2) % decibal value
 plot(thetaRange,fm,'b--')
hold on

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

% optimiztionX = [16.2333 172.3717 -117.6549 -162.0345 -24.2158 -37.9609 -102.4292 -99.0187...
%     -6.3424 -88.8941 68.3311 70.0410 131.9329 74.5378 26.5732 76.2868...
%     0.3696 0.7798 0.8206 0.6955 1.0000 1.0000 0.5534 0.4261...
%     0.6564 0.3807 0.0455 0.1865 0.3125 0.1814 0.1559 0.1640];