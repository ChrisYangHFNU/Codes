function y = Fun_Zmn02(z)
global beta a zzm zzn Deltz;
j = sqrt(-1);
R0 = sqrt((z-zzm).^2+a^2);
R1 = sqrt((z-(zzm-Deltz)).^2+a^2);
R2 = sqrt((z-(zzm+Deltz)).^2+a^2);
y = sin(beta*(z-(zzn-Deltz))).*(exp(-j*beta*R1)./R1-2*cos(beta*Deltz)*exp(-j*beta*R0)./R0+exp(-j*beta*R2)./R2);

