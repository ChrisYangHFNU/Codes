% Chapter9_polarizer_freq.m

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
f1        = 0.6 * terahertz;
f2        = 3.0 * terahertz;
NFREQ     = 40;
FREQ      = linspace(f1,f2,NFREQ);
SRC.theta = 0*degrees;

% DEFINE POLARIZER PARAMETERS
a = 4000 * nanometers;
t = 200 * nanometers;
w = 2000 * nanometers;
d = 1500 * nanometers;

ersi = 11.8;
erau = 1e6;

% DEFINE FDFD PARAMETERS
NRES     = 100;
NDIM     = 1;
lam0     = c0/mean(FREQ);
SPACER   = lam0*[0.1 0.1];
DEV.NPML = [10 10];
nmax     = sqrt(ersi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID RESOLUTION
lam_min = c0/max(FREQ);
dx      = lam_min/nmax/NRES;
dy      = lam_min/nmax/NRES;

% RESOLVE MINIMUM DIMENIONS
nx = w/dx;
if nx<NDIM
    dx = w/NDIM;
end
    
ny = t/dy;
if ny<NDIM
    dy = t/NDIM;
end
    
% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(a/dx);
dx = a/nx;
ny = ceil(d/dy);
dy = d/ny;

% GRID SIZE
Sx = a;
Nx = ceil(Sx/dx);
Sx = Nx*dx;

Sy = SPACER(1) + t + d + t + SPACER(2);
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
%% BUILD POLARIZER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO AIR
DEV.ER2 = ones(Nx2,Ny2);
DEV.UR2 = ones(Nx2,Ny2);
 
% COMPUTE ARRAY INDICES
nx  = 2*round(w/dx2/2);
nx1 = 1 + 2*floor((Nx2 - nx)/4);
nx2 = nx1 + nx - 1;

ny1 = 2*DEV.NPML(1) + 2*round(SPACER(1)/dy2/2) + 1;
ny2 = ny1 + 2*round(t/dy2/2);
ny3 = ny1 + round((t + d)/dy2);
ny4 = ny3 + 2*round(t/dy2/2);

% BUILD DEVICE
DEV.ER2(:,ny3:ny4) = erau;
DEV.ER2(:,ny4+1:Ny2) = ersi;
DEV.ER2(nx1:nx2,ny1:ny2) = erau;
DEV.ER2(nx1:nx2,ny2+1:Ny2) = ersi;
 
% SHOW DIFFRACTION GRATING
subplot(151);
imagesc(xa2,ya2,real(DEV.ER2).');
axis equal tight off;
colorbar;
title('ER2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SIMULATE DEVICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FINISH INPUT ARGUMENTS
DEV.RES  = [dx dy];

% INITIALIZE RECORDS
REF_H = zeros(1,NFREQ);
REF_E = zeros(1,NFREQ);
TRN_H = zeros(1,NFREQ);
TRN_E = zeros(1,NFREQ);
CON_H = zeros(1,NFREQ);
CON_E = zeros(1,NFREQ);

%
% WAVELENGTH SWEEP
%
for nfreq = 1 : NFREQ
    
    % Set Wavelength
    SRC.lam0 = c0/FREQ(nfreq);
    
    % Simulate E Mode
    SRC.MODE     = 'E';
    DAT_E        = fdfd2d(DEV,SRC);
    REF_E(nfreq) = DAT_E.REF;
    TRN_E(nfreq) = DAT_E.TRN;
    CON_E(nfreq) = DAT_E.REF + DAT_E.TRN;
    
    % Simulate H Mode
    SRC.MODE     = 'H';
    DAT_H        = fdfd2d(DEV,SRC);
    REF_H(nfreq) = DAT_H.REF;
    TRN_H(nfreq) = DAT_H.TRN;
    CON_H(nfreq) = DAT_H.REF + DAT_H.TRN;
    
    % Show E Mode 
    subplot(231);
    imagesc(xa,ya,real(DAT_E.f).');
    axis equal tight;
    colorbar;
    title('Field');
    
    subplot(2,3,[2 3]);
    semilogy(FREQ(1:nfreq)/terahertz,CON_E(1:nfreq),':k');
    hold on;
    semilogy(FREQ(1:nfreq)/terahertz,TRN_E(1:nfreq),'--b');
    semilogy(FREQ(1:nfreq)/terahertz,REF_E(1:nfreq),'-r');
    hold off;
    xlim([FREQ(1) FREQ(NFREQ)]/terahertz);
    ylim([-0.05 1.05]);
    
    % Show H Mode 
    subplot(234);
    imagesc(xa,ya,real(DAT_H.f).');
    axis equal tight;
    colorbar;
    title('Field');
    
    subplot(2,3,[5 6]);
    semilogy(FREQ(1:nfreq)/terahertz,CON_H(1:nfreq),':k');
    hold on;
    semilogy(FREQ(1:nfreq)/terahertz,TRN_H(1:nfreq),'--b');
    semilogy(FREQ(1:nfreq)/terahertz,REF_H(1:nfreq),'-r');
    hold off;
    xlim([FREQ(1) FREQ(NFREQ)]/terahertz);
    ylim([-0.05 1.05]);
    
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE AND PLOT EXTINCTION RATIO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE EXTINCTION RATIO
ER    = TRN_H./TRN_E;
ER_dB = 10*log10(ER);

% PLOT EXTINCTION RATIO
clf;
plot(FREQ/terahertz,ER_dB,'-k');
xlabel('Frequency (THz)');
ylabel('Extinction Ratio');