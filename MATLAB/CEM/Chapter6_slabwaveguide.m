% Chapter6_slabwaveguide.m

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
lam0 = 1.55 * micrometers;
MODE = 'H';

% SLAB WAVEGUIDE
a  = 1500 * nanometers;
n1 = 1.0;
n2 = 2.0;
n3 = 1.5;

% GRID PARAMETERS
nmax = max([n1 n2 n3]);
NRES = 20;
b    = 3*lam0;

% NUMBER OF MODES TO CALCULATE
NMODES = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIRST GUESS AT GRID RESOLUION
dx = lam0/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(a/dx);
dx = a/nx;

% CALCULATE TOTAL GRID SIZE
Sx = b + a + b;
Nx = ceil(Sx/dx);
Sx = Nx*dx;

% GRID AXIS
xa = [1:Nx]*dx;
xa = xa - mean(xa);

% 2X GRID
Nx2 = 2*Nx;
dx2 = dx/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD SLAB WAVEGUIDE ONTO GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO FREE SPACE
ER2 = ones(Nx2,1);
UR2 = ones(Nx2,1);

% CALCULATE START AND STOP ARRAY INDICES
nx1 =   1 + ceil(b/dx2);
nx2 = nx1 + round(a/dx2) - 1;

% BUILD SLAB WAVEGUIDE
ER2(1:nx1-1)   = n1^2;
ER2(nx1:nx2)   = n2^2;
ER2(nx2+1:Nx2) = n3^2;

% EXTRACT YEE GRID ARRAYS
ERxx = ER2(2:2:Nx2);
ERyy = ER2(1:2:Nx2);
ERzz = ER2(1:2:Nx2);
URxx = UR2(1:2:Nx2);
URyy = UR2(2:2:Nx2);
URzz = UR2(2:2:Nx2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FINITE-DIFFERENCE ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FORM DIAGONAL MATERIALS ARRAYS
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));
URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

% BUILD DERIVATIVE MATRICES
k0  = 2*pi/lam0;
NS  = [Nx 1];
RES = [dx 1];
BC  = [0 0];
[DEX,DEY,DHX,DHY] = yeeder2d(NS,k0*RES,BC);

% BUILD EIGEN-VALUE PROBLEM
if MODE=='E'
    A = -(DHX/URzz*DEX + ERyy);
    B = inv(URxx);
else
    A = -(DEX/ERzz*DHX + URyy);
    B = inv(ERxx);
end

% SOLVE EIGEN-VALUE PROBLEM
ev      = -n2^2;
[Fy,D2] = eigs(A,B,NMODES,ev);
D       = sqrt(D2);
NEFF    = -1i*diag(D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUALIZE THE GUIDED MODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PREPARE THE FIGURE WINDOW
clf;
hold on;

% DRAW CORE
x = NMODES*[0 1 1 0 0] - 0.5;
y = 0.5*a*[-1 -1 +1 +1 -1];
fill(x,y,0.8*[1 1 1]);

% DRAW THE MODES
Fy = Fy/max(abs(Fy(:)));
for m = 1 : NMODES
    x0 = m - 1;
    line(x0+0.4*Fy(:,m),xa,'LineWidth',3,'Color','w');
    line(x0+0.4*Fy(:,m),xa,'LineWidth',1,'Color','b');
    T = ['$n_{' num2str(m) ',\textrm{eff}} = ' ...
         num2str(NEFF(m),'%4.2f') '$'];
    text(x0,-b/2,T,'Interpreter','LaTex',...
         'Rotation',-90,...
         'VerticalAlignment','bottom');
end

% SET VIEW
hold off;
axis equal tight off;