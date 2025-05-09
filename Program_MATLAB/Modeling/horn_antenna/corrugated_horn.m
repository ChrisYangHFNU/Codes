clear;
close all;

% 设计参数 (示例值)
throat_diam = 10e-3; % 喉口直径 (m)
aperture_diam = 50e-3; % 口面直径 (m)
horn_length = 80e-3; % 喇叭总长度 (m)
corr_start_z = 20e-3; % 波纹开始位置 (从喉口开始)
corr_end_z = horn_length; % 波纹结束位置 (到口面结束)
corr_depth = 5e-3; % 波纹深度
slot_width = 2e-3; % 波纹槽宽度
tooth_width = 2e-3; % 波纹齿宽度
point_density = 1000; % 每米生成点数 (用于控制轮廓平滑度)

throat_radius = throat_diam / 2;
aperture_radius = aperture_diam / 2;

% 定义直线型廓形函数 R(z)
% R(z) = m * z + c
m = (aperture_radius - throat_radius) / horn_length;
c = throat_radius;
R_outer = @(z) m * z + c;

% 生成光滑段轮廓点 (从 z=0 到波纹开始位置)
z_smooth = linspace(0, corr_start_z, round(corr_start_z * point_density));
r_smooth = R_outer(z_smooth);

% 生成波纹段轮廓点
z_corr_range = corr_end_z - corr_start_z;
num_periods = floor(z_corr_range / (slot_width + tooth_width));
z_corr_points = [];
r_corr_points = [];

for i = 0:num_periods-1
    z_start_period = corr_start_z + i * (slot_width + tooth_width);

    % 波纹齿
    z_tooth_start = z_start_period;
    z_tooth_end = z_start_period + tooth_width;
    z_tooth = linspace(z_tooth_start, z_tooth_end, max(2, round(tooth_width * point_density)));
    r_tooth = R_outer(z_tooth);

    % 波纹槽 (底部和侧壁)
    z_slot_start = z_tooth_end;
    z_slot_end = z_slot_start + slot_width;

    % 槽内壁点 (底部)
    z_slot_bottom = linspace(z_slot_start, z_slot_end, max(2, round(slot_width * point_density)));
    r_slot_bottom = R_outer(z_slot_bottom) - corr_depth; % 槽底部半径

    % 槽侧壁点 (连接齿和槽底)
    r_outer_at_slot_start = R_outer(z_slot_start);
    r_bottom_at_slot_start = r_outer_at_slot_start - corr_depth;
    r_outer_at_slot_end = R_outer(z_slot_end);
    r_bottom_at_slot_end = r_outer_at_slot_end - corr_depth;

    z_side_points = [z_slot_start, z_slot_start, z_slot_end, z_slot_end];
    r_side_points = [r_outer_at_slot_start, r_bottom_at_slot_start, r_bottom_at_slot_end, r_outer_at_slot_end];


    % 合并当前周期内的点
    % 为了形成连续的轮廓，我们需要按顺序连接点
    if i == 0
         % 第一个周期从光滑段结束开始
         z_corr_points = [z_tooth, z_side_points(2), z_side_points(3), z_side_points(4)];
         r_corr_points = [r_tooth, r_side_points(2), r_side_points(3), r_side_points(4)];
    else
        % 后续周期连接上一个周期的结束点
        z_corr_points = [z_corr_points, z_tooth(2:end), z_side_points(2), z_side_points(3), z_side_points(4)];
        r_corr_points = [r_corr_points, r_tooth(2:end), r_side_points(2), r_side_points(3), r_side_points(4)];
    end

end

% 合并光滑段和波纹段点
% 确保点是唯一的且按 z 排序
z_profile = [z_smooth, z_corr_points];
r_profile = [r_smooth, r_corr_points];

% 移除重复点并按 Z 坐标排序
profile_points = unique([z_profile', r_profile'], 'rows');
profile_points = sortrows(profile_points, 1);

% 添加喉口中心点 (z=0, r=0) 和口面中心点 (z=horn_length, r=0) 以闭合轮廓，方便在AEDT中创建面
profile_points = [[0, 0]; profile_points; [horn_length, 0]];


% 提取最终的 z 和 r 坐标
final_z = profile_points(:, 1);
final_r = profile_points(:, 2);

% 可以选择绘制轮廓以检查
figure;
plot(final_z, final_r, '-o');
xlabel('Z (m)');
ylabel('R (m)');
title('Horn R-Z Profile');
axis equal;
grid on;

% 将轮廓点导出到文本文件
filename = 'horn_profile.csv';
% 创建一个包含 Z 和 R 列的表格
T = table(final_z, final_r, zeros(size(final_z)), 'VariableNames', {'X', 'Y', 'Z'}); % 导出为 X, Y, Z 格式，Z轴为喇叭轴线，所以 R 在这里对应 Y
writetable(T, filename, 'Delimiter', ',');

disp(['Profile points exported to ', filename]);