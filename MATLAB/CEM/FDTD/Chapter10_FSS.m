% Chapter10_FSS.m

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
f1    = 8 * gigahertz;
f2    = 12 * gigahertz;
NFREQ = 50;
FREQ  = linspace(f1,f2,NFREQ);
SRC.theta = 30 * degrees;
SRC.phi   = 60 * degrees;
SRC.pte   = 1;
SRC.ptm   = 0;

% DEFINE FSS PARAMETERS
a  = 9.95 * millimeters;   % period of FSS
s1 = 1.00 * millimeters;
s2 = 0.90 * millimeters;
w  = 4.50 * millimeters;
d  = 2.32 * millimeters;
h  = 3.16 * millimeters;

erd = 2.5;                  % substrate dielectric
erm = 1e6;                  % FSS metal

% DEFINE GRID PARAMETERS
NRESxy   = 100;
NRESz    = NRESxy;
DEV.NPML = [10 10];
nmax     = sqrt(erd);
lam_max  = c0/min(FREQ);
SPACER   = 0.2*lam_max * [1 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE PRELIMINARY RESOLUTION
lam_min = c0/max(FREQ);
dx      = lam_min/nmax/NRESxy;
dy      = lam_min/nmax/NRESxy;
dz      = lam_min/nmax/NRESz;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(a/dx);
dx = a/nx;

ny = ceil(a/dy);
dy = a/ny;

nz = ceil(h/dz);
dz = h/nz;

% CALCULATE GRID SIZE
Nx = ceil(a/dx);
Sx = Nx*dx;

Ny = ceil(a/dy);
Sy = Ny*dy;

Sz = SPACER(1) + h + SPACER(2);
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

% ADD SUBSTRATE
nz1 = 2*DEV.NPML(1) + 2*round(SPACER(1)/dz2/2) + 1;
nz2 = nz1 + round(h/dz2) - 1;
DEV.ER2(:,:,nz1:nz2) = erd;

% BUILD FSS UNIT CELL
L  = 2*s2 + 2*d + s1;
nx = round(L/dx);
x1 = (a - L)/2;
x2 = x1 + s2;
x3 = (a - w)/2;
x4 = (a - s1)/2;
x5 = (a + s1)/2;
x6 = (a + w)/2;
x8 = (a + L)/2;
x7 = x8 - s2;

nx1 = round(x1/dx2);
nx2 = round(x2/dx2);
nx3 = round(x3/dx2);
nx4 = round(x4/dx2);
nx5 = round(x5/dx2);
nx6 = round(x6/dx2);
nx7 = round(x7/dx2);
nx8 = round(x8/dx2);

ER2 = zeros(Nx2,Ny2);
ER2(nx1:nx8-1,nx4:nx5-1) = 1;
ER2(nx4:nx5-1,nx1:nx8-1) = 1;
ER2(nx1:nx2-1,nx3:nx6-1) = 1;
ER2(nx7:nx8-1,nx3:nx6-1) = 1;
ER2(nx3:nx6-1,nx1:nx2-1) = 1;
ER2(nx3:nx6-1,nx7:nx8-1) = 1;

ER2 = 1 + (erm - 1)*ER2;

% ADD UNIT CELL TO GRID
DEV.ER2(:,:,nz1) = ER2;

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
    
    % Show Response
    subplot(121);
    plot(FREQ(1:nfreq)/gigahertz,CON(1:nfreq),'--k');
    hold on;
    plot(FREQ(1:nfreq)/gigahertz,REF(1:nfreq),'-r');
    plot(FREQ(1:nfreq)/gigahertz,TRN(1:nfreq),'-b');
    hold off;
    xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
    ylim([-0.05 1.05]);
    title('Frequency Response');
    xlabel('Frequency (GHz)');
    
    subplot(122);
    TdB = 10*log10(TRN);
    plot(FREQ(1:nfreq)/gigahertz,TdB(1:nfreq),'-k');
    xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
    ylim([-40 4]);
    title('Frequency Response');
    xlabel('Frequency (GHz)');

    drawnow;

end