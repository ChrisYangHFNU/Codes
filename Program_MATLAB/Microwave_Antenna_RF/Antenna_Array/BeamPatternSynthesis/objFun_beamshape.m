% 定义适应度函数
function y = objFun_beamshape(x)

thetaRange = -90:1:90;
eleNum = 20; % 直线阵列阵元数
d = 0.5; % 归一化半波长
% Amp = ones(1,eleNum); % 初始化幅度值
% x = zeros(1,eleNum); % 初始化相位

D_theta = objPattern(thetaRange); % 目标方向图 余割平方
fm = zeros(1,length(thetaRange));
for ii=1:eleNum
    fm = fm+x(ii+eleNum)*exp(1j*(ii-1)*2*pi*d*sind(thetaRange)+1j*x(ii)); % 将相位值放在行向量x前半部分，幅度值放在行向量后半部分
end
fm = abs(fm./max(abs(fm))); % 归一化实现的方向图
% % ---------------约束主瓣3dB波束宽度-------------------
% index = find(abs(fm)==max(abs(fm))); % 找到方向图中最大值所在位置
% % % 向前寻找第一个最大值的0.707倍的值，即为波束宽度的下限
% target = 0.707;
% idx1 = 0;
% idx2 = 0;
% for ii=index-1:-1:1
%     if abs(fm(ii)-target)<=0.05
%         idx1 = ii;
%         break;
%     end
% end
% 
% % 向后寻找第一个最大值0.707倍的值，即为波束宽度的上限
% for jj=index+1:1:length(fm)
%     if abs(fm(jj)-target)<=0.05
%         idx2 = jj;
%         break;
%     end
% end
% bw = (idx2-idx1);

% --------------- 约束副瓣 --------------------
TF = islocalmax(fm); % 找出所有极大值点
TF = TF(find(TF<90 | TF>150)); % 排除主瓣内的所有极大值点
sll_m = max(fm(TF)); % 找出主瓣区外最大副瓣值
% fm_tmp = fm(TF);
% idx_local_max = find(fm_tmp<max(fm_tmp)) ; % 刨掉最大值点
% sll_m = max(fm_tmp(idx_local_max));  % 找到剩下的极大值点中最大值，即原数组中的第二大值点
% ---------------------------------------------

w = ones(1,length(thetaRange));
w(91:148) = 3; 
% y = sum(w.*abs(fm-D_theta));
if sll_m<=0.1
    y = sum(w.*abs(fm-D_theta))+0*abs(sll_m-0.1);
else
    y = sum(w.*abs(fm-D_theta))+6*abs(sll_m-0.1);
end

% y(1) = sum(w.*abs(fm-D_theta));
% y(2) = abs(sll_m-0.1);
% y(3) = abs(bw-35);

end