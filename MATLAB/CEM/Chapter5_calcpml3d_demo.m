% Chapter5_calcpml3d_demo.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% DEFINE INPUT ARGUMENTS
NGRID = [6 6 6];
NPML  = [2 2 2 2 2 2];

% CALL ADDUPML2D
[sx,sy,sz] = calcpml3d(NGRID,NPML);

% DISPLAY THE RESULTS
disp('sx =');
disp(sx);

disp('sy =');
disp(sy);

disp('sz =');
disp(sz);


















