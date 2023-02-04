% Chapter10_GMRF.m

% INITIALIZE MATLAB
close all;
clc
clear all;

% DEFINE UNITS
meters      = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
micrometers = 1e-6 * meters;
nanometers  = 1e-9 * meters;
seconds     = 1;
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;
degrees     = pi/180;

% CONSTANTS
c0 = 299792458 * meters/seconds;
e0 = 8.854187812813e-12 * 1/meters;
u0 = 1.256637062121e-6 * 1/meters;
N0 = 376.7303136686;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE SOURCE PARAMETERS
f1        = 4.5 * gigahertz;
f2        = 5.5 * gigahertz;
NFREQ     = 500;
FREQ      = linspace(f1,f2,NFREQ);
SRC.theta = 0 * degrees;
SRC.phi   = 0 * degrees;
SRC.pte   = 1/sqrt(2);
SRC.ptm   = 1i/sqrt(2);

% DEFINE GMRF PARAMETERS
er1 = 1.0;
er2 = 2.0;
erL = 3.0;
erH = 5.0;
a   = 3.7 * centimeters;
t   = 1.5 * centimeters;
f   = 0.5;

% DEFINE GRID PARAMETERS
NRES     = 20;
DEV.NPML = [10 10];
ermax    = max([er1 er2 erL erH]);
nmax     = sqrt(ermax);
lam_max  = c0/min(FREQ);
SPACER   = lam_max * [1 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE PRELIMINARY RESOLUTION
lam_min = c0/max(FREQ);
dx      = lam_min/nmax/NRES;
dy      = lam_min/nmax/NRES;
dz      = lam_min/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(a/dx);
dx = a/nx;
ny = ceil(a/dy);
dy = a/ny;
nz = ceil(t/dz);
dz = t/nz;

% CALCULATE GRID SIZE
Nx = ceil(a/dx);
Sx = Nx*dx;
Ny = ceil(a/dy);
Sy = Ny*dy;
Sz = SPACER(1) + t + SPACER(2);
Nz = DEV.NPML(1) + ceil(Sz/dz) + DEV.NPML(2);
Sz = Nz*dz;

% 2x GRID
Nx2 = 2*Nx;                 dx2 = dx/2;
Ny2 = 2*Ny;                 dy2 = dy/2;
Nz2 = 2*Nz;                 dz2 = dz/2;

% GRID AXES
xa = [1:Nx]*dx;             xa = xa - mean(xa);
ya = [1:Ny]*dy;             ya = ya - mean(ya);
za = [0.5:Nz-0.5]*dz;

xa2 = [1:Nx2]*dx2;          xa2 = xa2 - mean(xa2);
ya2 = [1:Ny2]*dy2;          ya2 = ya2 - mean(ya2);
za2 = [0.5:Nz2-0.5]*dz2;

% MESHGRIDS
[Y,X,Z]    = meshgrid(ya,xa,za);
[Y2,X2,Z2] = meshgrid(ya2,xa2,za2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD MATERIALS ARRAYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO FREE SPACE
DEV.ER2 = ones(Nx2,Ny2,Nz2);
DEV.UR2 = ones(Nx2,Ny2,Nz2);

% BUILD GRATING LAYER
r = a*sqrt(f/pi);
DEV.ER2 = (X2.^2 + Y2.^2) <= r^2;
DEV.ER2 = erL + (erH - erL)*DEV.ER2;

% ADD SUPERSTRATE AND SUBSTRATE
nz1 = 2*DEV.NPML(1) + round(SPACER(1)/dz2) + 1;
nz2 = nz1 + round(t/dz2) - 1;
DEV.ER2(:,:,1:nz1-1)  = er1;
DEV.ER2(:,:,nz2+1:Nz2) = er2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FREQUENCY SWEEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEV
DEV.RES = [dx dy dz];

% INITIALIZE RESPONSE
REF = zeros(1,NFREQ);
TRN = zeros(1,NFREQ);
CON = zeros(1,NFREQ);

%
% MAIN LOOP -- ITERATE OVER FREQUENCY
%
for nfreq = 1 : NFREQ
    
    % Get Nxt Frequency
    f0       = FREQ(nfreq);
    SRC.lam0 = c0/f0;
    
    % Call FDFD3D()
    DAT = fdfd3d(DEV,SRC);
        
    % Calculate Overall Response
    REF(nfreq) = DAT.REF;
    TRN(nfreq) = DAT.TRN;
    CON(nfreq) = DAT.CON;
    
    % Show Field
    subplot(131);
    slice(Y/centimeters,X/centimeters,Z/centimeters,...
          real(DAT.fy),0,0,0);
    shading interp;
    axis equal tight;
    view(250,20);
    colormap('gray');
    caxis(1.5*[-1 1]);
    
    % Show Response
    subplot(1,3,[2 3]);
    plot(FREQ(1:nfreq)/gigahertz,CON(1:nfreq),'--k');
    hold on;
    plot(FREQ(1:nfreq)/gigahertz,REF(1:nfreq),'-r');
    plot(FREQ(1:nfreq)/gigahertz,TRN(1:nfreq),'-b');
    hold off;
    xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
    ylim([-0.05 1.05]);
    title('Frequency Response');
    xlabel('Frequency (GHz)');
    drawnow;
    
end