function [ xy_array, r_phi_array ] = get_circular_r_phi(a, n_points)
%GET_CIRCULAR_R_PHI Summary of this function goes here
%   Detailed explanation goes here
    x_arr = linspace(-a,a,n_points);
    y_arr = linspace(-a,a,n_points);
    xy_array = zeros(length(x_arr)*length(y_arr), 2);
    i = 1;
    for ix = 1:length(x_arr)
        for iy = 1:length(y_arr)
            xy_array(i, 1) = x_arr(ix);
            xy_array(i, 2) = y_arr(iy);
            i = i+1;
        end
    end

    i_outoff_wg = find(sqrt(xy_array(:, 1).^2 + xy_array(:, 2).^2) > a);
    xy_array(i_outoff_wg, :) = [];

    r_phi_array = zeros(size(xy_array));
    r_phi_array(:,1) = sqrt(xy_array(:,1).^2+xy_array(:,2).^2);
    r_phi_array(:,2) = atan2(xy_array(:,2),xy_array(:,1));

end

