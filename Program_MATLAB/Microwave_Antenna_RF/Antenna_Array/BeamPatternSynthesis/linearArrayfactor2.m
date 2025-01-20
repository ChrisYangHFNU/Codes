function fm = linearArrayfactor2(eleNum,thetaRange,eleSpace,x)
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
plot(thetaRange,20*log10(fm),'b--') % decibal value
ylim([-40 0])
end