% Chapter10_metamaterial.m

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
f1        = 9 * gigahertz;
f2        = 14 * gigahertz;
NFREQ     = 200;
FREQ      = linspace(f1,f2,NFREQ);
SRC.theta = 0 * degrees;
SRC.phi   = 0 * degrees;
SRC.pte   = 0;
SRC.ptm   = 1;

% DEFINE METAMATERIAL PARAMETERS
a  = 2.500 * millimeters;   % period of metamaterial
h  = 0.250 * millimeters;   % substrate thickness
w1 = 0.140 * millimeters;   % width of wire
w2 = 0.200 * millimeters;   % linewidth of ring
g  = 0.300 * millimeters;   % gap in ring
L2 = 2.200 * millimeters;   % length of outer ring
L1 = 1.500 * millimeters;   % length of inner ring

er   = 4.4;
tand = 0.02;
erd  = er*(1 - 1i*tand);

f0    = 10 * gigahertz;
sigma = 5.8e7;
erm   = 1 + sigma/(1i*2*pi*f0*e0);

% DEFINE GRID PARAMETERS
NRES     = 120;
DEV.NPML = [10 10];
nmax     = sqrt(erd);
lam0     = c0/mean(FREQ);
SPACER   = 0.1*lam0 * [1 1];

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
nz = ceil(a/dz);
dz = a/nz;

% CALCULATE GRID SIZE
Nx = ceil(a/dx);
Sx = Nx*dx;

Ny = ceil(a/dy);
Sy = Ny*dy;

Sz = SPACER(1) + a + SPACER(2);
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
ny  = 2*round(h/dy2/2);
ny1 = 1 + 2*floor((Ny2 - ny)/4);
ny2 = ny1 + ny - 1;
nz  = round(a/dz2);
nz1 = 1 + 2*floor((Nz2 - nz)/2/2);
nz2 = nz1 + nz - 1;
DEV.ER2(:,ny1:ny2,nz1:nz2) = erd;

% ADD WIRE
nza = nz1 + 2*round((a-w1)/2/dz2/2);
nzb = nza + round(w1/dz2) - 1;
DEV.ER2(:,ny1,nza:nzb) = erm;

% ADD OUTER RING
ny  = ny2 + 1;
nz  = 2*round(L2/dz2/2);
nz1 = 1 + 2*floor((Nz2 - nz)/4);
nz2 = nz1 + nz - 1;
nx  = round(L2/dz2);
nx1 = 1 + 2*floor((Nx2 - nx)/4);
nx2 = nx1 + nx - 1;
DEV.ER2(nx1:nx2,ny,nz1:nz2) = erm;
nxa = nx1 + round(w2/dz2) - 1;
nxb = nx2 - round(w2/dz2) + 1;
nz1 = nz1 + round(w2/dz2) - 1;
nz2 = nz2 - round(w2/dz2) + 1;
DEV.ER2(nxa:nxb,ny,nz1:nz2) = 1;
nz  = 2*round(g/dz2/2);
nz1 = 2*floor((Nz2 - nz)/4);
nz2 = nz1 + nz - 1;
DEV.ER2(nx1:nxb,ny,nz1:nz2) = 1;

% ADD INNER RING
ny  = ny2 + 1;
nz  = 2*round(L1/dz2/2);
nz1 = 1 + 2*floor((Nz2 - nz)/4);
nz2 = nz1 + nz - 1;
nx  = round(L1/dz2);
nx1 = 1 + 2*floor((Nx2 - nx)/4);
nx2 = nx1 + nx - 1;
DEV.ER2(nx1:nx2,ny,nz1:nz2) = erm;
nxa = nx1 + round(w2/dz2) - 1;
nxb = nx2 - round(w2/dz2) + 1;
nz1 = nz1 + round(w2/dz2) - 1;
nz2 = nz2 - round(w2/dz2) + 1;
DEV.ER2(nxa:nxb,ny,nz1:nz2) = 1;
nz  = 2*round(g/dz2/2);
nz1 = 2*floor((Nz2 - nz)/4);
nz2 = nz1 + nz - 1;
DEV.ER2(nxa:nx2,ny,nz1:nz2) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FREQUENCY SWEEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEV
DEV.RES = [dx dy dz];

% INITIALIZE RESPONSE
S11 = zeros(1,NFREQ);
S21 = zeros(1,NFREQ);
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
    k0       = 2*pi/SRC.lam0;
    
    % Call FDFD3D()
    DAT = fdfd3d(DEV,SRC);
        
    % Record Response
    dsrc = SPACER(1) - 3*dz;
    dref = SPACER(1) - 1*dz;
    dtrn = SPACER(2) - 1*dz;
    nxc  = 1 + floor(Nx/2);
    nyc  = 1 + floor(Nx/2);
    S11(nfreq) = DAT.s11(nxc,nyc) ./ exp(-1i*k0*(dsrc + dref));
    S21(nfreq) = DAT.s21(nxc,nyc) ./ exp(-1i*k0*(dsrc + dtrn));
    REF(nfreq) = DAT.REF;
    TRN(nfreq) = DAT.TRN;
    CON(nfreq) = DAT.CON;
    
    % Show Response
    subplot(2,3,[2 3]);
    plot(FREQ(1:nfreq)/gigahertz,CON(1:nfreq),'--k');
    hold on;
    plot(FREQ(1:nfreq)/gigahertz,REF(1:nfreq),'-r');
    plot(FREQ(1:nfreq)/gigahertz,TRN(1:nfreq),'-b');
    hold off;
    xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
    ylim([-0.05 1.05]);
    title('Power Response');
    xlabel('Frequency (GHz)');
    
    subplot(2,3,[5 6]);
    plot(FREQ(1:nfreq)/gigahertz,angle(S11(1:nfreq)),'-r');
    hold on;
    plot(FREQ(1:nfreq)/gigahertz,angle(S21(1:nfreq)),'-b');
    hold off;
    xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
    title('Phase Response');
    xlabel('Frequency (GHz)');
    
    drawnow;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM PARAMETER RETRIEVAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CHANGE SIGN CONVENTION
S11 = conj(S11);
S21 = conj(S21);

% UNWRAP PHASE
phi = angle(S11);
for n = 1 : NFREQ-1
    dphi = phi(n+1) - phi(n);
    if abs(dphi) > 0.5*pi
        phi(n+1:NFREQ) = phi(n+1:NFREQ) - dphi;
    end
end
S11 = abs(S11).*exp(1i*phi);
phi11 = phi;

phi = angle(S21);
for n = 1 : NFREQ-1
    dphi = phi(n+1) - phi(n);
    if abs(dphi) > 0.5*pi
        phi(n+1:NFREQ) = phi(n+1:NFREQ) - dphi;
    end
end
S21 = abs(S21).*exp(1i*phi);
phi21 = phi;

% RETRIEVE IMPEDANCE AND REFRACTIVE INDEX
X = (1 - S21.^2 + S11.^2)./(2*S11);

s      = ones(1,NFREQ);
r      = X + s.*sqrt(X.^2 - 1);
ind    = find(abs(r)>1);
s(ind) = -1;
r      = X + s.*sqrt(X.^2 - 1);

t = (S11 + S21 - r)./(1 - (S11 + S21).*r);

k0   = 2*pi*FREQ/c0;
eta  = N0*(1 + r)./(1 - r);
neff = log(t)./(1i*k0*a);

% RETRIEVE PERMITTIVITY AND PERMEABILITY
er = neff*N0./eta;
ur = neff.*eta./N0;

% SHOW RESPONSE
subplot(321);
plot(FREQ/gigahertz,CON,'--k','LineWidth',2);
hold on;
plot(FREQ/gigahertz,REF,':k','LineWidth',2);
plot(FREQ/gigahertz,TRN,'-k','LineWidth',2);
hold off;
xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
ylim([-0.05 1.05]);
xlabel('Frequency (GHz)');
title('Reflectance & Transmittance');

subplot(322);
plot(FREQ/gigahertz,phi11,':k','LineWidth',2);
hold on;
plot(FREQ/gigahertz,phi21,'-k','LineWidth',2);
hold off;
xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
ylim([-2.25 2.25]);
xlabel('Frequency (GHz)');
title('Phase Response');

% SHOW REFRACTIVE INDEX AND IMPEDANCE
subplot(323);
plot(FREQ/gigahertz,real(neff),'-k','LineWidth',2);
hold on;
plot(FREQ/gigahertz,imag(neff),'--k','LineWidth',2);
hold off;
xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
ylim([-5 +7]);
xlabel('Frequency (GHz)');
title('Refractive Index $n$','Interpreter','LaTex');

subplot(324);
plot(FREQ/gigahertz,real(eta)/N0,'-k','LineWidth',2);
hold on;
plot(FREQ/gigahertz,imag(eta)/N0,'--k','LineWidth',2);
hold off;
xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
ylim([-3 +3]);
xlabel('Frequency (GHz)');
title('Impedance $\eta/\eta_0$','Interpreter','LaTex');

% SHOW RELATIVE PERMITTIVITY AND PERMEABILITY
subplot(325);
plot(FREQ/gigahertz,real(er),'-k','LineWidth',2);
hold on;
plot(FREQ/gigahertz,imag(er),'--k','LineWidth',2);
hold off;
xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
ylim([-5 2]);
xlabel('Frequency (GHz)');
title('Relative Permittivity $\varepsilon_r$',...
      'Interpreter','LaTex');

subplot(326);
plot(FREQ/gigahertz,real(ur),'-k','LineWidth',2);
hold on;
plot(FREQ/gigahertz,imag(ur),'--k','LineWidth',2);
hold off;
xlim([FREQ(1) FREQ(NFREQ)]/gigahertz);
ylim([-7 +14]);
xlabel('Frequency (GHz)');
title('Relative Permeability $\mu_r$','Interpreter','LaTex');