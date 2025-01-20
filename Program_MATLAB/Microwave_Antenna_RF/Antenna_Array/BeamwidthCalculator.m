%% This Code Plots Beamwidth vs. Frequency and Scan Angle 
% Arik D. Brown

%% Input Parameters
BW.k = 0.886; % Beamwidth Factor (radians)
BW.f_vec = [1 5 10 15]; %Frequency in GHz
BW.lambda_vec = 0.3./BW.f_vec; %meters
BW.L = 1; %Aperture Length in meters
BW.theta_vec = 0:5:60; %Degrees

%% Calculate Beamwidths
[BW.lambda_mat, BW.theta_mat] = meshgrid(BW.lambda_vec,BW.theta_vec);
BW.mat_rad = BW.k*BW.lambda_mat./(BW.L*cosd(BW.theta_mat));
BW.mat_deg = rad2deg(BW.mat_rad);

%% Plot
figure(1),clf
plot(BW.theta_mat,BW.mat_deg,'linewidth',2)
grid
set(gca,'fontsize',16,'fontweight','b')
xlabel('Scan Angle (Degrees)','fontsize',16,'fontweight','b')
legend('1 GHz','5 GHz','10 GHz','15 GHz')