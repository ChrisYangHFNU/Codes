% Chapter9_polarizer_d.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
micrometers = 1;
nanometers  = 1e-3 * micrometers;
millimeters = 1e3 * micrometers;
meters      = 1e6 * micrometers;
degrees     = pi/180;

seconds   = 1;
hertz     = 1/seconds;
kilohertz = 1e3 * hertz;
megahertz = 1e6 * hertz;
gigahertz = 1e9 * hertz;
terahertz = 1e12 * hertz;

% CONSTANTS
c0 = 299279458 * meters/seconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE PLANE WAVE SOURCE PARAMETERS
f0        = 2.0 * terahertz;
SRC.lam0  = c0/f0;
SRC.theta = 0*degrees;

% DEFINE POLARIZER PARAMETERS
a = 4000 * nanometers;
t = 200 * nanometers;
w = 2000 * nanometers;

ersi = 11.8;
erau = 1e6;

d1    = 500 * nanometers;
d2    = 10 * micrometers;
NDAT  = 50;
d_DAT = linspace(d1,d2,NDAT);

% DEFINE FDFD PARAMETERS
NRES     = 200;
NDIM     = 1;
SPACER   = SRC.lam0*[0.1 0.1];
DEV.NPML = [10 10];
nmax     = sqrt(ersi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID RESOLUTION
dx0 = SRC.lam0/nmax/NRES;
dy0 = SRC.lam0/nmax/NRES;

% RESOLVE MINIMUM DIMENIONS
nx = w/dx0;
if nx<NDIM
    dx0 = w/NDIM;
end
    
ny = t/dy0;
if ny<NDIM
    dy0 = t/NDIM;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM PARAMETER SWEEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE RECORDS
REF_H = zeros(1,NDAT);
REF_E = zeros(1,NDAT);
TRN_H = zeros(1,NDAT);
TRN_E = zeros(1,NDAT);
CON_H = zeros(1,NDAT);
CON_E = zeros(1,NDAT);
ER    = zeros(1,NDAT);
ER_dB = zeros(1,NDAT);

%
% MAIN LOOP -- ITERATE OVER d
%
for ndat = 1 : NDAT
    
    % Get Next Value of d
    d = d_DAT(ndat);

    % Snap Grid to Critical Dimensions
    nx = ceil(a/dx0);
    dx = a/nx;
    ny = ceil(d/dy0);
    dy = d/ny;

    % Calculate Grid Size
    Sx = a;
    Nx = ceil(Sx/dx);
    Sx = Nx*dx;
    Sy = SPACER(1) + t + d + t + SPACER(2);
    Ny = DEV.NPML(1) + ceil(Sy/dy) + DEV.NPML(2);
    Sy = Ny*dy;

    % Calculate 2x Grid Parameters
    Nx2 = 2*Nx;             dx2 = dx/2;
    Ny2 = 2*Ny;             dy2 = dy/2;

    % Calculate Axis Vectors
    xa  = [1:Nx]*dx;
    ya  = [1:Ny]*dy;
    xa2 = [1:Nx2]*dx2;
    ya2 = [1:Ny2]*dy2;

    % Initialize to Air
    DEV.ER2 = ones(Nx2,Ny2);
    DEV.UR2 = ones(Nx2,Ny2);

    % Calculate Array Indices
    nx  = 2*round(w/dx2/2);
    nx1 = 1 + 2*floor((Nx2 - nx)/4);
    nx2 = nx1 + nx - 1;

    ny1 = 2*DEV.NPML(1) + 2*round(SPACER(1)/dy2/2) + 1;
    ny2 = ny1 + 2*round(t/dy2/2);
    ny3 = ny1 + round((t + d)/dy2);
    ny4 = ny3 + 2*round(t/dy2/2);

    % Build Polarizer on 2x Grid
    DEV.ER2(:,ny3:ny4) = erau;
    DEV.ER2(:,ny4+1:Ny2) = ersi;
    DEV.ER2(nx1:nx2,ny1:ny2) = erau;
    DEV.ER2(nx1:nx2,ny2+1:Ny2) = ersi;

    % Finish DEV
    DEV.RES  = [dx dy];

    % Simulate E Mode
    SRC.MODE    = 'E';
    DAT_E       = fdfd2d(DEV,SRC);
    REF_E(ndat) = DAT_E.REF;
    TRN_E(ndat) = DAT_E.TRN;
    CON_E(ndat) = DAT_E.REF + DAT_E.TRN;

    % Simulate H Mode
    SRC.MODE    = 'H';
    DAT_H       = fdfd2d(DEV,SRC);
    REF_H(ndat) = DAT_H.REF;
    TRN_H(ndat) = DAT_H.TRN;
    CON_H(ndat) = DAT_H.REF + DAT_H.TRN;
    
    % Calculate and Record Extinction Ratio
    ER(ndat)    = TRN_H(ndat)/TRN_E(ndat);
    ER_dB(ndat) = 10*log10(ER(ndat));
    
    % Show E Mode 
    subplot(141);
    imagesc(xa,ya,real(DAT_E.f).');
    axis equal tight;
    colorbar;
    title('Field');

    % Show H Mode 
    subplot(142);
    imagesc(xa,ya,real(DAT_H.f).');
    axis equal tight;
    colorbar;
    title('Field');

    subplot(1,4,[3 4]);
    plot(d_DAT(1:ndat)/micrometers,ER_dB(1:ndat),'-k');
    xlim([d_DAT(1) d_DAT(NDAT)]/micrometers);
    ylim([50 70]);

    drawnow;
end