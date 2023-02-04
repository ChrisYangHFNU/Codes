% Chapter8_sawtooth.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
degrees     = pi/180;
meters      = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
micrometers = 1e-6 * meters;
nanometers  = 1e-9 * meters;
inches      = 2.54 * centimeters;
feet        = 12 * inches;
seconds     = 1;
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;
terahertz   = 1e12 * hertz;
petahertz   = 1e15 * hertz;

% CONSTANTS
e0 = 8.85418782e-12 * 1/meters;
u0 = 1.25663706e-6 * 1/meters;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE PLANE WAVE SOURCE PARAMETERS
f0    = 30.0 * gigahertz;
theta = 30*degrees;
lam0  = c0/f0;
k0    = 2*pi/lam0;
MODE  = 'E';

% DEFINE DIFFRACTION GRATING PARAMETERS
L   = 1.5 * centimeters;
d   = 1.0 * centimeters;
er1 = 1.0;
er2 = 9.0;

% DEFINE FDFD PARAMETERS
NRES   = 100;
SPACER = lam0*[1 1];
NPML   = [20 20];
ermax  = max([er1 er2]);
nmax   = sqrt(ermax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(L/dx);
dx = L/nx;

ny = ceil(d/dy);
dy = d/ny;

% GRID SIZE
Sx = L;
Nx = ceil(Sx/dx);
Sx = Nx*dx;

Sy = SPACER(1) + d + SPACER(2);
Ny = NPML(1) + ceil(Sy/dy) + NPML(1);
Sy = Ny*dy;

% 2X GRID
Nx2 = 2*Nx;             dx2 = dx/2;
Ny2 = 2*Ny;             dy2 = dy/2;

% CALCULATE AXIS VECTORS
xa = [1:Nx]*dx;
ya = [1:Ny]*dy;
[Y,X] = meshgrid(ya,xa);
xa2 = [1:Nx2]*dx2;
ya2 = [1:Ny2]*dy2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DIFFRACTION GRATING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO AIR
ER2 = er1*ones(Nx2,Ny2);
UR2 = ones(Nx2,Ny2);
 
% BUILD SAWTOOTH
ny1 = 2*NPML(1) + round(SPACER(1)/dy2) + 1;
ny2 = ny1 + round(d/dy2) - 1;
for ny = ny1 : ny2
    f = (ny - ny1 + 1)/(ny2 - ny1 + 1);
    nx = round(f*L/dx2);
    ER2(Nx2-nx+1:Nx2,ny) = er2;
end
 
% ADD SUBSTRATE
ER2(:,ny2+1:Ny2) = er2;
 
% SHOW DIFFRACTION GRATING
subplot(121);
imagesc(xa2,ya2,ER2.');
axis equal tight off;
colorbar;
title('ER2');

% INCORPORATE PML
[ERxx,ERyy,ERzz,URxx,URyy,URzz] ...
                    = addupml2d(ER2,UR2,[0 0 NPML]);

% DIAGONALIZE MATERIAL TENSORS
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));
URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FDFD ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INCIDENT WAVE VECTOR
nsrc  = sqrt(UR2(1,1)*ER2(1,1));
kxinc = k0*nsrc*sin(theta);
kyinc = k0*nsrc*cos(theta);
kinc  = [ kxinc ; kyinc ];

% BUILD DERIVATIVE MATRICES
NS  = [Nx Ny];
RES = [dx dy];
BC  = [1 0];
[DEX,DEY,DHX,DHY] = yeeder2d(NS,k0*RES,BC,kinc/k0);
        
% BUILD WAVE MATRIX
if MODE == 'E'
    A = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;
else
    A = DEX/ERyy*DHX + DEY/ERxx*DHY + URzz;
end

% CALCULATE SOURCE FIELD
fsrc  = exp(-1i*(kxinc*X + kyinc*Y));

% CALCULATE SCATTERED-FIELD MASKING MATRIX
ny = NPML(1) + 2;
Q  = zeros(Nx,Ny);
Q(:,1:ny) = 1;
Q  = diag(sparse(Q(:)));

% CALCULATE SOURCE VECTOR
b = (Q*A - A*Q)*fsrc(:);

% SOLVE FOR FIELD
f = A\b;
f = reshape(f,Nx,Ny);

% SHOW FIELD
subplot(122);
pcolor(xa,ya,real(f).');
shading interp;
axis equal tight;
set(gca,'YDir','reverse');
caxis([-1 1]);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYZE REFLECTION AND TRANSMISSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXTRACT MATERIAL PROPERTIES IN EXTERNAL REGIONS
urref = UR2(1,1);
urtrn = UR2(1,Ny2);
erref = ER2(1,1);
ertrn = ER2(1,Ny2);
nref  = sqrt(urref*erref);
ntrn  = sqrt(urtrn*ertrn);

% CALCULATE WAVE VECTOR COMPONENTS
m     = [-floor(Nx/2):+floor((Nx-1)/2)].';
kx    = kxinc - m*2*pi/Sx;
kyref = sqrt((k0*nref)^2 - kx.^2);
kytrn = sqrt((k0*ntrn)^2 - kx.^2);

% EXTRACT REFELECTED AND TRANSMITTED WAVES
fref = f(:,NPML(1)+1)./fsrc(:,1);
ftrn = f(:,Ny-NPML(2))./fsrc(:,1);

% CALCULATE AMPLITUDES OF SPATIAL HARMONICS
aref = fftshift(fft(fref))/Nx;
atrn = fftshift(fft(ftrn))/Nx;

% CALCULATE DIFFRACTION EFFICENCIES
RDE = abs(aref).^2.*real(kyref/kyinc);
if MODE == 'E'
    TDE = abs(atrn).^2.*real(urref/urtrn*kytrn/kyinc);
else
    TDE = abs(atrn).^2.*real(erref/ertrn*kytrn/kyinc);
end

% CALCULATE OVERALL REFLECTANCE & TRANSMITTANCE
REF = sum(RDE(:))
TRN = sum(TDE(:))
CON = REF + TRN

% SHOW RESULTS
ind = find((RDE(:)+TDE(:))>1e-4);
[ m(ind) 100*RDE(ind) 100*TDE(ind) ]