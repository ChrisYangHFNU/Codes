clear all
close all

C = 3e8;
f = 10e9;
lambda = C/f;
l = lambda;
omega = 2*pi*f;
beta = 2*pi/lambda;
d = linspace(0,2*l,100);
t = linspace(0,2*pi/omega,50);
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
savepath = 'D:\01.yangjing\Codes\GitHub\image\'; % 保存路径
for kk=1:LenT
    h = figure; % 将figure保存在h中
    set(h,'Visible','off'); %将figure设为不显示，防止一直弹框
    plot(d,vMatrix_sw(kk,:),d,iMatrix_sw(kk,:)) % 绘图
    axis([0 max(d) -2.5 2.5]);
    legend('电压驻波','电流驻波')
    saveas(h,[savepath 'standingwave-' num2str(kk-1)],'jpg');
    clf;
end

inputpath = 'D:\01.yangjing\Codes\GitHub\image\'; % 图片输入路径
format = '.jpg'; % 图片格式
pic = dir([inputpath,'*.jpg']); % 返回路径内容：文件名
WriteObj = VideoWriter('D:\01.yangjing\Codes\GitHub\movie\standingwave.avi'); % 合成视频目标文件路径
WriteObj.FrameRate = 5; % 调整帧率，用来调整视频长短
open(WriteObj); %打开视频
for ii=1:(length(pic))
    frame = imread(strcat(inputpath,'standingwave-',num2str(ii-1),format)); % 读取图片，放在变量frame中
    writeVideo(WriteObj,frame); % frame存到变量WriteObj中
    % ...这里是你想要对每张图片的操作，比如反色、二值化之类
end
close(WriteObj); %关闭视频
    
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