close all; clear; clc;
%% Definitions
phi_i       =   120;
% FileName    =   'Circle.dat';
% FileName    =   'Triangle.dat';
% FileName    =   'Square.dat';
% FileName    =   'Plane.dat';
FileName    =   'TwoCircles.dat';
%% Call the main function
[In_TE,sigma_TE,phi_TE,Data_TE]=Scat2D(FileName,phi_i,"TE");
fprintf('\n\n');
[In_TM,sigma_TM,phi_TM,Data_TM]=Scat2D(FileName,phi_i,"TM");
%% Plot totat fields approximation
fprintf('\n');
range       =   5;
PlotTE(phi_i,In_TE,Data_TE,range);
PlotTM(phi_i,In_TM,Data_TM,range);
%%
