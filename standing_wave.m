clear all
close all

C = 3e8;
f = 10e9;
lambda = C/f;
l = lambda;
omega = 2*pi*f;
beta = 2*pi/lambda;
d = linspace(0,l,100);
t = linspace(0,2*pi/omega,100);
LenT = length(t);
vMatrix_sw = [];
iMatrix_sw = [];
vMatrix_tw = [];
iMatrix_tw = [];
for k=1:LenT
    v_sw = -2*sin(beta*d)*sin(omega*t(k));
    i_sw = 2*cos(beta*d)*cos(omega*t(k));
    v_tw = cos(omega*t(k)-beta*d);
    i_tw = cos(omega*t(k)-beta*d);
    vMatrix_sw = [vMatrix_sw;v_sw];
    iMatrix_sw = [iMatrix_sw;i_sw];
    vMatrix_tw = [vMatrix_tw;v_tw];
    iMatrix_tw = [iMatrix_tw;i_tw];
end
figure(1)
p = plot(vMatrix_sw,'EraseMode','background','MarkerSize',5);
axis([0 100 -2.5 2.5])
for kk=1:LenT
    set(p,'XData',[1:LenT],'YData',vMatrix_sw)
    drawnow
    pause(0.1);
end

% for kk=1:LenT
%     plot(vMatrix_sw(kk,:),'b')
%     plot(iMatrix_sw(kk,:),'r')
%     axis([0 100 -2.5 2.5]);
%     hold on
%     pause(0.05);
% end
% figure(2)
% for kkk=1:LenT
%     plot(vMatrix_tw(kkk,:),'b')
%     %plot(iMatrix_tw(kkk,:),'r')
%     axis([0 100 -1.5 1.5]);
%     hold on
%     pause(0.05);
% end