% Chapter4_yeeder3d_demo.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% DEFINE INPUT ARGUMENTS FOR YEEDER2D
NS  = [2 2 2];
RES = [0.3 0.2 0.1];
BC  = [0 0 0];

% CALL YEEDER3D
[DEX,DEY,DEZ,DHX,DHY,DHZ] = yeeder3d(NS,RES,BC);

% SHOW DERIVATIVE MATRICES
disp('DEX = ');
disp(full(DEX));

disp('DEY = ');
disp(full(DEY));

disp('DEY = ');
disp(full(DEZ));

disp('DHX = ');
disp(full(DHX));

disp('DHY = ');
disp(full(DHY));

disp('DHZ = ');
disp(full(DHZ));







































