close all;
clear;

import1 = csvread('S11.csv',1,0);
import1_1 = csvread('S11An2.csv',1,0);
% import2 = csvread('Eplane.csv',1,0);
% import2_1 = csvread('EplaneAn2.csv',1,0);
% import3 = csvread('Hplane.csv',1,0);
% import3_1 = csvread('HplaneAn2.csv',1,0);
import4 = csvread('gainAn1.csv',1,0);
import5 = csvread('gainAn2.csv',1,0);
import6 = csvread('inputZ1.csv',1,0);
import7 = csvread('inputZ2.csv',1,0);

freq = import1(:,1);
S11 = import1(:,2);
S11An2 = import1_1(:,2);

%{ 
phi_deg = import2(:,1);phi_deg2 = import2_1(:,1);
phi_rad = deg2rad(phi_deg);phi_rad2 = deg2rad(phi_deg2);
gain3g_db = import2(:,2);gain3g_db2 = import2_1(:,2);
gain6g_db = import2(:,3);gain6g_db2 = import2_1(:,3);
gain9g_db = import2(:,4);gain9g_db2 = import2_1(:,4);

theta_deg = import3(:,1);theta_deg2 = import3_1(:,1);
theta_rad = deg2rad(theta_deg);theta_rad2 = deg2rad(theta_deg2);
gain3g_dbH = import3(:,2);gain3g_db2H = import3_1(:,2);
gain6g_dbH = import3(:,3);gain6g_db2H = import3_1(:,3);
gain9g_dbH = import3(:,4);gain9g_db2H = import3_1(:,4);
%}

f1 = import4(:,3);
realisedGainAn1 = import4(:,4);
f2 = import5(:,3);
realisedGainAn2 = import5(:,4);

fre1 = import6(:,1);fre2 = import7(:,1);
R1 = import6(:,3);X1 = import6(:,2);
R2 = import7(:,3);X2 = import7(:,2);

h1 = figure(1);
plot(freq,S11,'-',freq,S11An2,'--','Linewidth',2,'Color',[64/255,105/255,224/255]);hold on;
plot(freq,ones(1,length(freq))*(-10),'k','LineWidth',2);hold on;
% plot(2.4*ones(1,length(freq(5:end))),S11(5:end),'k-.','LineWidth',1.5);

set(gca,'FontSize',10);

xlabel('频率/GHz','interpreter','latex');
ylabel('反射系数/dB','Interpreter','latex');

legend('天线1','天线2','Orientation','horizontal');
set(legend,'Location','Northeast');

xlim([2,10]);
% ylim([-30,-5]);
xticks(2:1:10);

%text(2.4,-10,'2.4GHz','FontSize',15);

%{
h2 = figure(2);
polarplot(phi_rad,gain3g_db,'b-',phi_rad,gain6g_db,'b--',phi_rad,gain9g_db,'b-.','LineWidth',2);
rlim([-20,10])

set(gca,'FontSize',15);
title('天线1的E面场强方向图');

legend('3GHz','6GHz','9GHz','Orientation','vertical');
set(legend,'Location','NortheastOutside')

h2_1 = figure(3);
polarplot(phi_rad2,gain3g_db2,'b-',phi_rad2,gain6g_db2,'b--',phi_rad2,gain9g_db2,'b-.','LineWidth',2);
rlim([-30,15]);

set(gca,'FontSize',15);
title('天线2的E面场强方向图');

legend('3GHz','6GHz','9GHz','Orientation','vertical');
set(legend,'Location','NortheastOutside')

h3 = figure(4);
polarplot(theta_rad,gain3g_dbH,'b-',theta_rad,gain6g_dbH,'b--',theta_rad,gain9g_dbH,'b-.','LineWidth',2);
rlim([-20,10]);

set(gca,'FontSize',15);
title('天线1的H面场强方向图');
legend('3GHz','6GHz','9GHz','Orientation','vertical');
set(legend,'Location','NortheastOutside')

h3_1 = figure(5);
polarplot(theta_rad2,gain3g_db2H,'b-',theta_rad2,gain6g_db2H,'b--',theta_rad2,gain9g_db2H,'b-.','LineWidth',2);
rlim([-20,10]);

set(gca,'FontSize',15);
title('天线2的H面场强方向图');
legend('3GHz','6GHz','9GHz','Orientation','vertical');
set(legend,'Location','NortheastOutside')
%}

h4 = figure(6);
plot(f1,realisedGainAn1,'-',f2,realisedGainAn2,'-.','Linewidth',2,'Color',[64/255,105/255,224/255]);

set(gca,'FontSize',10);
xlabel('频率/GHz','Interpreter','latex');
ylabel('增益/dBi','Interpreter','latex');
legend('天线1','天线2','Orientation','horizontal');
set(legend,'Location','Northeast')
xticks(2:1:10);

h5 = figure(7);
plot(fre1,R1,'-',fre2,R2,'-.','Linewidth',2,'Color',[64/255,105/255,224/255]);
hold on
plot(fre1,X1,'-',fre2,X2,'-.','Linewidth',2,'Color','k');hold on;
plot(freq,ones(1,length(freq))*(50),'r',freq,ones(1,length(freq))*(0),'r','LineWidth',2);

set(gca,'FontSize',10);
xlabel('频率/GHz','Interpreter','latex');
ylabel('输入阻抗/Ω','Interpreter','latex');
legend('天线1阻抗虚部','天线2阻抗虚部','天线1阻抗实部','天线2阻抗实部','Orientation','vertical');
set(legend,'Location','Northeast')
xticks(2:1:10);

% 设置输出图大小
% set(gcf,'PaperUnits','centimeters')
% set(gcf,'PaperSize',[28,11.4])
% set(gcf,'PaperPositionMode','manual')
% set(gcf,'PaperPosition',[0,0,28,11.4])
% set(gcf,'Renderer','painters');

% print fig_s11.eps -depsc2 -r600
print(h1,'fig_s11.jpg','-djpeg','-r600');
% 
% print(h2,'fig_Eplane.jpg','-djpeg','-r600')
% print(h3,'fig_Hplane.jpg','-djpeg','-r600')
% print(h2_1,'fig_Eplane2.jpg','-djpeg','-r600')
% print(h3_1,'fig_Hplane2.jpg','-djpeg','-r600')
print(h4,'fig_Gain.jpg','-djpeg','-r600')
print(h5,'fig_InputZ.jpg','-djpeg','-r600')