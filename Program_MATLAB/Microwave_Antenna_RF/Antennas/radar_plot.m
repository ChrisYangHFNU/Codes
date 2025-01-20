clc;clear;close all;
% 参数位置
% 雷达位置
radarLon = 113.3; % longtitude
radarLat = 28.1;  % latitude
radarHgt = 1000;  % height(m)

% 读取目标轨迹
load 'targetWGS84.txt';
targetLon = targetWGS84(:,2).'; % 目标经度（degree）
targetLat = targetWGS84(:,3).'; % 目标纬度
targetHgt = targetWGS84(:,4).'; % 目标高度(m)

% 绘制目标相对于雷达运动的轨迹
figure
plot3(radarLon,radarLat,radarHgt,'o','MarkerFaceColor','red','MarkerSize',6);
hold on;
for n=1:size(targetWGS84,1)
    plot3(targetLon,targetLat,targetHgt,'r-*'); % 使用plot3绘制三维图
    lonlon = [radarLon,targetLon(n)]; % [雷达经度，目标点经度]
    latlat = [radarLat,targetLat(n)]; % [雷达纬度，目标点纬度]
    hgthgt = [radarHgt,targetHgt(n)]; % [雷达高度，目标点高度]
    plot3(lonlon,latlat,hgthgt,'g');
end
hold off;
xlabel('经度(°)');
ylabel('纬度(°)');
zlabel('高度(m)');
title('目标轨迹');