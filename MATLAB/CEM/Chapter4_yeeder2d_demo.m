% Chapter4_yeeder2d_demo.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% DEFINE INPUT ARGUMENTS FOR YEEDER2D
NS  = [3 4];
RES = [0.2 0.1];
BC  = [0 0];

% CALL YEEDER2D
[DEX,DEY,DHX,DHY] = yeeder2d(NS,RES,BC);

% SHOW DERIVATIVE MATRICES
disp('DEX = ');
disp(full(DEX))

disp('DEY = ');
disp(full(DEY))

disp('DHX = ');
disp(full(DHX))

disp('DHY = ');
disp(full(DHY))







































