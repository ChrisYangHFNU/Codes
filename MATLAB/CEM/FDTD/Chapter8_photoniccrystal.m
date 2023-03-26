% Chapter8_photoniccrystal.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
degrees = pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE
lam0  = 1;
k0    = 2*pi/lam0;
bw    = 1*lam0;
theta = 20*degrees;

% PHOTONIC CRYSTAL PARAMETERS
wn  = 0.217;
a   = wn*lam0;
r   = 0.48*a;
er1 = 1.0;
er2 = 12.0;
Lx  = 60*a;
Ly  = 40*a;

% FDFD PARAMETERS
NRES   = 20;
SPACER = lam0*[6 6 0 0];
NPML   = [20 20 20 20];
ermax  = max([er1 er2]);
nmax   = sqrt(ermax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
b  = a*sqrt(2);
nx = ceil(b/dx);
dx = b/nx;

ny = ceil(b/dy);
dy = b/ny;

% GRID SIZE
Sx = SPACER(1) + Lx + SPACER(2);
Nx = NPML(1) + ceil(Sx/dx) + NPML(2);
Sx = Nx*dx;

Sy = SPACER(3) + Ly + SPACER(4);
Ny = NPML(3) + ceil(Sy/dy) + NPML(4);
Sy = Ny*dy;

% 2X GRID
Nx2 = 2*Nx;             dx2 = dx/2;
Ny2 = 2*Ny;             dy2 = dy/2;

% CALCULATE AXIS VECTORS
xa = [1:Nx]*dx;         
xa = xa - NPML(1)*dx - SPACER(1);
ya = [1:Ny]*dy;         
ya = ya - mean(ya);

xa2 = [1:Nx2]*dx2;      
xa2 = xa2 - 2*NPML(1)*dx2 - SPACER(1);
ya2 = [1:Ny2]*dy2;      
ya2 = ya2 - mean(ya2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD LATTICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UNIT CELL GRID
Nxuc = round(b/dx2);
Nyuc = round(b/dy2);
xauc = [1:Nxuc]*dx2;      xauc = xauc - mean(xauc);
yauc = [1:Nyuc]*dy2;      yauc = yauc - mean(yauc);
[Y,X] = meshgrid(yauc,xauc);

% BUILD UNIT CELL
ERuc = (X.^2 + Y.^2) <= r^2;
ERuc = ERuc | ((X - b/2).^2 + (Y - b/2).^2) <= r^2;
ERuc = ERuc | ((X + b/2).^2 + (Y - b/2).^2) <= r^2;
ERuc = ERuc | ((X - b/2).^2 + (Y + b/2).^2) <= r^2;
ERuc = ERuc | ((X + b/2).^2 + (Y + b/2).^2) <= r^2;
ERuc = er2 + (er1 - er2)*ERuc;
URuc = ones(Nxuc,Nyuc);

% INITIALIZE MATERIALS ARRAYS
ER2 = ones(Nx2,Ny2);
UR2 = ones(Nx2,Ny2);

% STACK UNIT CELL TO BUILD LATTICE
for nyuc = 1 : round(Ly/b)
    ny1 = 2*NPML(3) + (nyuc - 1)*Nyuc ...
        + round(SPACER(3)/dy2) + 1;
    ny2 = ny1 + Nyuc - 1;
    for nxuc = 1 : round(Lx/b)
        nx1 = 2*NPML(1) + (nxuc - 1)*Nxuc ...
            + round(SPACER(1)/dx2) + 1;
        nx2 = nx1 + Nxuc - 1;
        ER2(nx1:nx2,ny1:ny2) = ERuc;
        UR2(nx1:nx2,ny1:ny2) = URuc;
    end
end

% SHOW UNIT CELL
subplot(131);
imagesc(xa2,ya2,ER2.');
axis equal tight;
colorbar;
title('Lattice');

% INCORPORATE PML
[ERxx,ERyy,ERzz,URxx,URyy,URzz] = addupml2d(ER2,UR2,NPML);

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

% BUILD DERIVATIVE MATRICES
NS  = [Nx Ny];
RES = [dx dy];
BC  = [0 0];
[DEX,DEY,DHX,DHY] = yeeder2d(NS,k0*RES,BC);
        
% BUILD WAVE MATRIX
A = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;

% CALCULATE SOURCE FIELD
nsrc      = sqrt(UR2(1,1)*ER2(1,1));
[Y,X]     = meshgrid(ya,xa);
[TH,R]    = cart2pol(X,Y);
[XR,YR]   = pol2cart(TH-theta,R);
fsrc      = exp(-(YR/bw).^2).*exp(-1i*k0*nsrc*XR);
fsrc(X>0) = 0;

% CALCULATE SCATTERED-FIELD MASKING MATRIX
nx1 = NPML(1) + 2;
nx2 = Nx - NPML(2) - 1;
ny1 = NPML(3) + 2;
ny2 = Ny - NPML(4) - 1;
Q   = ones(Nx,Ny);
Q(nx1:nx2,ny1:ny2) = 0;
Q  = diag(sparse(Q(:)));

% CALCULATE SOURCE VECTOR
b = (Q*A - A*Q)*fsrc(:);

% SOLVE FOR FIELD
f = A\b;
f = reshape(f,Nx,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DRAW PRETTY SOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PREPARE GRAPHICS WINDOW
clf;
hold on;

% DRAW FIELD
pcolor(xa,ya,real(f).');
shading interp;

% SUPERIMPOSE LATTICE
h = contour(xa2,ya2,ER2.',10,'Color','w','LineWidth',0.5);

% SET GRAPHICS VIEW
hold off;
axis equal tight;
set(gca,'YDir','reverse');
caxis([-1 1]);
colorbar;