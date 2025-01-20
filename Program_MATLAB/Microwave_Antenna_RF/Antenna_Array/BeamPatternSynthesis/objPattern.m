function D_theta = objPattern(thetaRange)

D_theta = zeros(1,length(thetaRange));
for ii=1:length(thetaRange)
    if thetaRange(ii)<=-2
        D_theta(ii) = 0.01;
    elseif thetaRange(ii)<0
        D_theta(ii) = 0.4*thetaRange(ii)+0.8;
    elseif thetaRange(ii)<4
        D_theta(ii) = 0.05*thetaRange(ii)+0.8;
    elseif thetaRange(ii)<8
        D_theta(ii) = 1;
    elseif thetaRange(ii)<58
        D_theta(ii) = -0.01*thetaRange(ii)+1.08;
    else
        D_theta(ii) = 0.01;
    end
end
plot(thetaRange,D_theta)

end