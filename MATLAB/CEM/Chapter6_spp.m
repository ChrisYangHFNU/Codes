% Chapter6_spp.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS 
micrometers = 1;
nanometers  = 1e-3 * micrometers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WAVELENGTH AND MODE
lam0 = 500 * nanometers;

% SURFACE PROPERTIES
erd = 2.31;
erm = -9.98 - 0.26i;
b1  = 2*lam0;
b2  = 1*lam0;

% GRID PARAMETERS
nmax   = max(real([sqrt(erd) sqrt(erm)]));
NRES   = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIRST GUESS AT GRID RESOLUION
dx = lam0/nmax/NRES;

% CALCULATE TOTAL GRID SIZE
Sx = b1 + b2;
Nx = ceil(Sx/dx);
Sx = Nx*dx;

% GRID AXIS
xa = [0:Nx-1]*dx;
xa = xa - b1;

% 2X GRID
Nx2 = 2*Nx;
dx2 = dx/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD MEDIUMS ONTO THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO FREE SPACE
ER2 = ones(Nx2,1);
UR2 = ones(Nx2,1);

% BUILD MATERIAL INTERFACE
nx = 2*round(b1/dx2/2);
ER2(1:nx)     = erd;
ER2(nx+1:Nx2) = erm;

% EXTRACT YEE GRID ARRAYS
ERxx = ER2(2:2:Nx2);
ERzz = ER2(1:2:Nx2);
URyy = UR2(2:2:Nx2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FINITE-DIFFERENCE ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FORM DIAGONAL MATERIALS ARRAYS
ERxx = diag(sparse(ERxx(:)));
ERzz = diag(sparse(ERzz(:)));
URyy = diag(sparse(URyy(:)));

% BUILD DERIVATIVE MATRICES
k0  = 2*pi/lam0;
NS  = [Nx 1];
RES = [dx 1];
BC  = [0 0];
[DEX,DEY,DHX,DHY] = yeeder2d(NS,k0*RES,BC);

% BUILD EIGEN-VALUE PROBLEM (H MODE)
A = -(DEX/ERzz*DHX + URyy);
B = inv(ERxx);

% SOLVE EIGEN-VALUE PROBLEM
ev      = -erd*erm/(erd + erm);
[Fy,D2] = eigs(A,B,1,ev);
D       = sqrt(D2);
NEFF    = -1i*diag(D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUALIZE THE SURFACE WAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PREPARE THE FIGURE WINDOW
clf;
hold on;

% DRAW MEDIUMS
x = 1.2*[0 1 1 0 0] - 0.1;
y = [0 0 -b1 -b1 0];
fill(x,y,0.9*[1 1 1]);
y = [0 0 b2 b2 0];
fill(x,y,0.5*[1 1 1]);

% DRAW THE MODES
Fy = real(Fy)/max(abs(Fy));
line(Fy,xa,'LineWidth',3,'Color','w');
line(Fy,xa,'LineWidth',1,'Color','b');

% SET VIEW
hold off;
axis equal tight off;
set(gca,'YDir','reverse');
title(['$n_{\rm{eff}} = ' num2str(NEFF,'%4.2f') '$'],...
      'Interpreter','LaTex');