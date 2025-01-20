clc;
close all;
clear ;

f = -10:.01:10;
w = 2*pi*f;
fl = 3;
fk = 4;
wl = fl*2*pi;
wk = fk*2*pi;
xil = 0.8;
xik = 0.9;
Hk = 1i*w.^0./(w.^2-wk^2+2*1i*xik*w);
Hl = 1i*w.^0./(w.^2-wl^2+2*1i*xil*w);
fun1 = Hk.*conj(Hl);
wu = 4.5;
n2 = 1000-wu/0.01;
figure
plot(f,abs(fun1),'LineWidth',2,'Color','k')
hold on;fill([f(n2:end-n2),f(end-n2),f(n2)],[abs(fun1(n2:end-n2)),0,0],[0.8,0.8,0.8])
set(gca,'FontSize',15); % 横纵轴标注字体大小
xlabel('Frequency','FontSize',15,'FontWeight','bold');
ylabel('H(\omega)','FontSize',15,'FontWeight','bold');
set(gca,'XTick',[-4.5,0,4.5]); % X轴的记号位置
set(gca,'XTickLabel',{'-\omega_u','0','\omega_u'}); % X轴的记号
set(gca,'YTick',[0]); % Y轴的记号位置
set(gca,'YTicklabel',{'0'}); % Y轴的记号
set(gca,'YLim',[0 1.1*max(abs(fun1))]); % Y轴的数据显示范围
text(6,0.3*max(abs(fun1)),'I(\omega_u)','FontSize',15);
hold on;plot([4,6],[0.2*max(abs(fun1)),0.25*max(abs(fun1))],'k');