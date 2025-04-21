% 波束端口4 入射 天线端口的幅度相位分布,10GHz

tmp11 = csvread('phaB1.csv',1,0); % unit:rad
pha1 = tmp11(11,2:end);
tmp12 = csvread('magB1.csv',1,0);% unit:dB
mag1_dB = tmp2(11,2:end);
mag1 = 10.^(mag1_dB/10);
x1 = [pha1 mag1];

tmp21 = csvread('phaB2.csv',1,0);
pha2 = tmp21(11,2:end);
tmp22 = csvread('magB2.csv',1,0);
mag2_dB = tmp22(11,2:end);
mag2 = 10.^(mag2_dB/10);
x2 = [pha2 mag2];

tmp31 = csvread('phaB3.csv',1,0);
pha3 = tmp31(11,2:end);
tmp32 = csvread('magB3.csv',1,0);
mag3_dB = tmp32(11,2:end);
mag3 = 10.^(mag3_dB/10);
x3 = [pha3 mag3];

tmp41 = csvread('phaB4.csv',1,0);
pha4 = tmp41(11,2:end);
tmp42 = csvread('magB4.csv',1,0);
mag4_dB = tmp42(11,2:end);
mag4 = 10.^(mag4_dB/10);
x4 = [pha4 mag4];

eleNum = 8;
elespace = 0.5;
thetaRange = -90:90;
fm1 = linearArrayfactor(eleNum,thetaRange,elespace,x1);
fm2 = linearArrayfactor(eleNum,thetaRange,elespace,x2);
fm3 = linearArrayfactor(eleNum,thetaRange,elespace,x3);
fm4 = linearArrayfactor(eleNum,thetaRange,elespace,x4);
plot(thetaRange,20*log10(fm1),'k-.',thetaRange,20*log10(fm2),'k--',thetaRange,20*log10(fm3),'k-',thetaRange,20*log10(fm4),'k-*','LineWidth',2)
xlim([-90 90]);
ylim([-40 0]);