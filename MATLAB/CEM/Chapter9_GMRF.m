% Chapter9_GMRF.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
micrometers = 1;
nanometers  = 1e-3 * micrometers;
degrees     = pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE PLANE WAVE SOURCE PARAMETERS
lam1      = 1.4 * micrometers;
lam2      = 1.6 * micrometers;
NLAM      = 1000;
LAMBDA    = linspace(lam1,lam2,NLAM);
SRC.theta = 0*degrees;
SRC.MODE  = 'E';

% DEFINE GMRF PARAMETERS
f  = 0.5;
n1 = 1.0;
n2 = 1.4;
nL = 1.9;
nH = 2.1;
d1 = 265 * nanometers;
d2 = 380 * nanometers;
L  = 870 * nanometers;

% DEFINE FDFD PARAMETERS
NRES     = 40;
SPACER   = max(LAMBDA)*[0.5 0.5];
DEV.NPML = [20 20];
nmax     = max([n1 n2 nL nH]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID RESOLUTION
dx = min(LAMBDA)/nmax/NRES;
dy = min(LAMBDA)/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(L/dx);
dx = L/nx;
ny = ceil(d2/dy);
dy = d2/ny;

% GRID SIZE
Sx = L;
Nx = ceil(Sx/dx);
Sx = Nx*dx;

Sy = SPACER(1) + d1 + d2 + d1 + SPACER(2);
Ny = DEV.NPML(1) + ceil(Sy/dy) + DEV.NPML(2);
Sy = Ny*dy;

% 2X GRID
Nx2 = 2*Nx;             dx2 = dx/2;
Ny2 = 2*Ny;             dy2 = dy/2;

% CALCULATE AXIS VECTORS
xa  = [1:Nx]*dx;
ya  = [1:Ny]*dy;
xa2 = [1:Nx2]*dx2;
ya2 = [1:Ny2]*dy2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD GUIDED-MODE RESONANCE FILTER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO AIR
DEV.ER2 = n1^2*ones(Nx2,Ny2);
DEV.UR2 = ones(Nx2,Ny2);
 
% COMPUTE ARRAY INDICES
nx  = round(f*L/dx2);
nx1 = 1 + floor((Nx2 - nx)/2);
nx2 = nx1 + nx - 1;

ny1 = 2*DEV.NPML(1) + round(SPACER(1)/dy2) + 1;
ny2 = ny1 + round(d1/dy2) - 1;
ny3 = ny2 + 1;
ny4 = ny3 + round(d2/dy2) - 1;
ny5 = ny4 + 1;
ny6 = ny5 + round(d1/dy2) - 1;

% BUILD DEVICE
DEV.ER2(:,ny1:ny2)       = n2^2;
DEV.ER2(:,ny3:ny4)       = nL^2;
DEV.ER2(nx1:nx2,ny3:ny4) = nH^2;
DEV.ER2(:,ny5:ny6)       = n2^2;
 
% SHOW DIFFRACTION GRATING
subplot(151);
imagesc(xa2,ya2,DEV.ER2.');
axis equal tight off;
colorbar;
title('ER2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATE DEVICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FINISH INPUT ARGUMENTS
DEV.RES  = [dx dy];

% INITIALIZE RECORDS
REF = zeros(1,NLAM);
TRN = zeros(1,NLAM);
CON = zeros(1,NLAM);

%
% WAVELENGTH SWEEP
%
for nlam = 1 : NLAM
    
    % Set Wavelength
    SRC.lam0 = LAMBDA(nlam);
    
    % Call fdfd2d()
    DAT = fdfd2d(DEV,SRC);

    % Record Results
    REF(nlam) = DAT.REF;
    TRN(nlam) = DAT.TRN;
    CON(nlam) = DAT.REF + DAT.TRN;
    
    % Show Field
    subplot(152);
    imagesc(xa,ya,real(DAT.f).');
    axis equal tight;
    colorbar;
    title('Field');
    
    % Show Spectra
    subplot(1,5,[3:5]);
    plot(LAMBDA(1:nlam)/micrometers,CON(1:nlam),':k');
    hold on;
    plot(LAMBDA(1:nlam)/micrometers,TRN(1:nlam),'--b');
    plot(LAMBDA(1:nlam)/micrometers,REF(1:nlam),'-r');
    hold off;
    xlim([LAMBDA(1) LAMBDA(NLAM)]/micrometers);
    ylim([-0.05 1.05]);
    drawnow;
end