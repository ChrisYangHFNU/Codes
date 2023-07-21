clear all;
close all;

import1 = csvread('S11.csv',1,0);
import2 = csvread('Eplane.csv',1,0);
import3 = csvread('Hplane.csv',1,0);

freq = import1(:,1);
S11 = import1(:,2);

phi_deg = import2(:,1);
phi_rad = deg2rad(phi_deg);
gain3g_db = import2(:,2);
gain6g_db = import2(:,3);
gain9g_db = import2(:,4);
gain3g = 10.^(gain3g_db/20);
gain6g = 10.^(gain6g_db/20);
gain9g = 10.^(gain9g_db/20);

theta_deg = import3(:,1);
theta_rad = deg2rad(theta_deg);
gain3g_dbH = import3(:,2);
gain6g_dbH = import3(:,3);
gain9g_dbH = import3(:,4);
gain3gH = 10.^(gain3g_dbH/20);
gain6gH = 10.^(gain6g_dbH/20);
gain9gH = 10.^(gain9g_dbH/20);

h1 = figure(1);
plot(freq,S11,'Linewidth',2,'Color',[64/255,105/255,224/255]);hold on;
plot(freq,ones(1,length(freq))*(-10),'k','LineWidth',2);hold on;
plot(2.4*ones(1,length(freq(5:end))),S11(5:end),'k-.','LineWidth',1.5);

set(gca,'FontSize',15);

xlabel('$Frequency$/GHz','interpreter','latex');
ylabel('$|S_{11}|$/dB','Interpreter','latex');

legend('Antenna1','Orientation','horizontal');
set(legend,'Location','Northeast');

xlim([2,10]);
ylim([-30,-5]);
xticks([2,2.4,3:10]);

%text(2.4,-10,'2.4GHz','FontSize',15);

h2 = figure(2);
polarplot(phi_rad,gain3g,'b-',phi_rad,gain6g,'b--',phi_rad,gain9g,'b-.','LineWidth',2);

set(gca,'FontSize',15);
title('E面场强方向图');
xlabel('$\varphi$','Interpreter','latex');

legend('3GHz','6GHz','9GHz','Orientation','vertical');
set(legend,'Location','NortheastOutside')

h3 = figure(3);
polarplot(theta_rad,gain3gH,'b-',theta_rad,gain6gH,'b--',theta_rad,gain9gH,'b-.','LineWidth',2);

set(gca,'FontSize',15);
title('H面场强方向图');
legend('3GHz','6GHz','9GHz','Orientation','vertical');
set(legend,'Location','NortheastOutside')

% 设置输出图大小
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperSize',[28,11.4])
set(gcf,'PaperPositionMode','manual')
set(gcf,'PaperPosition',[0,0,28,11.4])
set(gcf,'Renderer','painters');

% print fig_s11.eps -depsc2 -r600
% print(h1,'fig_s11.jpg','-djpeg','-r600');

% set(gcf,'PaperUnits','centimeters')
% set(gcf,'PaperSize',[28,11.4])
% set(gcf,'PaperPositionMode','manual')
% set(gcf,'PaperPosition',[0,0,28,11.4])
% set(gcf,'Renderer','painters');

% print(h2,'fig_Eplane.jpg','-djpeg','-r600')
% print(h3,'fig_Hplane.jpg','-djpeg','-r600')