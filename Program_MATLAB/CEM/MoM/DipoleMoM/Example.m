close all; clear; clc;
%%
L       =   0.47;   % Length per wavelength
a       =   0.005;  % Radius per wavelength
N       =   10;     % Number of segments 
%% Call the MM solution
[I,Zin,theta_E,E_theta_E,phi_A,E_theta_A]=DipoleMM(L,a,N);
%%
