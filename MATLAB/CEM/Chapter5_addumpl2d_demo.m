% Chapter5_addumpl2d_demo.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% DEFINE GRID
Nx = 5;
Ny = 5;

% DEFINE PML
NPML = [2 2 2 2];

% BUILD ER2 AND UR2 ARRAYS
ER2 = ones(2*Nx,2*Ny);
UR2 = ones(2*Nx,2*Ny);

% CALL ADDUPML2D
[ERxx,ERyy,ERzz,URxx,URyy,URzz] ...
              = addupml2d(ER2,UR2,NPML);

% DISPLAY THE RESULTS
disp('ERxx =');
disp(ERxx);

disp('ERyy =');
disp(ERyy);

disp('ERzz =');
disp(ERzz);

disp('URxx =');
disp(URxx);

disp('URyy =');
disp(URyy);

disp('URzz =');
disp(URzz);

















