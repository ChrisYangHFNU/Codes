clear 
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
savepath = 'D:\7.Codes\GitHub\image\'; % ����·��
for kk=1:LenT
    h = figure; % ��figure������h��
    set(h,'Visible','off'); %��figure��Ϊ����ʾ����ֹһֱ����
    plot(d,vMatrix_tw(kk,:),d,iMatrix_tw(kk,:)) % ��ͼ
    axis([0 max(d) -2.5 2.5]);
    legend('��ѹ�в�','�����в�')
    saveas(h,[savepath 'travelingwave-' num2str(kk-1)],'jpg');
    clf;
end

inputpath = 'D:\7.Codes\GitHub\image\'; % ͼƬ����·��
format = '.jpg'; % ͼƬ��ʽ
pic = dir([inputpath,'*.jpg']); % ����·�����ݣ��ļ���
WriteObj = VideoWriter('D:\7.Codes\GitHub\movie\travelingwave.avi'); % �ϳ���ƵĿ���ļ�·��
WriteObj.FrameRate = 5; % ����֡�ʣ�����������Ƶ����
open(WriteObj); %����Ƶ
for ii=1:(length(pic))
    frame = imread(strcat(inputpath,'travelingwave-',num2str(ii-1),format)); % ��ȡͼƬ�����ڱ���frame��
    writeVideo(WriteObj,frame); % frame�浽����WriteObj��
    % ...����������Ҫ��ÿ��ͼƬ�Ĳ��������練ɫ����ֵ��֮��
end
close(WriteObj); %�ر���Ƶ
    
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