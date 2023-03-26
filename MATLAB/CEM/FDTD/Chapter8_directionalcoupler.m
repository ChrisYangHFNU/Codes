% Chapter8_directionalcoupler.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
micrometers = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE
lam0  = 1.55 * micrometers;
k0    = 2*pi/lam0;

% RIB WAVEGUIDE PARAMETERS
rib_n1 = 1.0;
rib_n2 = 3.5;
rib_n3 = 1.5;
rib_h  = 0.30 * micrometers;
rib_t  = 0.10 * micrometers;
rib_w  = 0.40 * micrometers;

% DIRECTIONAL COUPLER PARAMETERS
L1 = 0.5*lam0;
L2 = 5.0*lam0;
g  = 0.1*lam0;

% FDFD PARAMETERS
NRES   = 20;
SPACER = lam0*[2 2 1 1];
NPML   = [20 20 20 20];
nmax   = max([rib_n1 rib_n2 rib_n3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(rib_w/dx);
dx = rib_w/nx;
ny = ceil(rib_w/dy);
dy = rib_w/ny;

% GRID SIZE
Sx = SPACER(1) + L1 + L2 + SPACER(2);
Nx = NPML(1) + ceil(Sx/dx) + NPML(2);
Sx = Nx*dx;

Sy = SPACER(3) + rib_w + g + rib_w + SPACER(4);
Ny = NPML(3) + ceil(Sy/dy) + NPML(4);
Sy = Ny*dy;

% 2X GRID
Nx2 = 2*Nx;             dx2 = dx/2;
Ny2 = 2*Ny;             dy2 = dy/2;

% CALCULATE AXIS VECTORS
xa = [1:Nx]*dx;
ya = [1:Ny]*dy;
xa2 = [1:Nx2]*dx2;
ya2 = [1:Ny2]*dy2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE EFFECTIVE INDICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ANALYZE OFF WAVEGUIDE

    % Define Waveguide
    a = rib_t;
    b = 5*lam0;
    
    % Calculate 1D Grid for Slab
    dz  = lam0/nmax/NRES;
    nz  = ceil(a/dz);
    dz  = a/nz;
    Sz  = b + a + b;
    Nz  = ceil(Sz/dz);
    Sz  = Nz*dz;
    za  = [1:Nz]*dz;
    za  = za - mean(za);
    Nz2 = 2*Nz;
    dz2 = dz/2;
    
    % Build Slab
    ER2            = ones(Nz2,1);
    UR2            = ones(Nz2,1);
    nz1            = 1 + ceil(b/dz2);
    nz2            = nz1 + round(a/dz2) - 1;
    ER2(1:nz1-1)   = rib_n1^2;
    ER2(nz1:nz2)   = rib_n2^2;
    ER2(nz2+1:Nz2) = rib_n3^2;
    
    % Build Eigen-Value Problem
    ERxx = diag(sparse(ER2(2:2:Nz2)));
    ERzz = diag(sparse(ER2(1:2:Nz2)));
    URyy = diag(sparse(UR2(2:2:Nz2)));
    
    NS  = [1 Nz];
    RES = [1 dz];
    BC  = [0 0];
    [~,DEZ,~,DHZ] = yeeder2d(NS,k0*RES,BC);

    A = -(DEZ/ERxx*DHZ + URyy);
    B = inv(ERzz);

    % Solve Eigen-Value Problem
    ev         = -rib_n2^2;
    [Hy,D2]    = eigs(A,B,4,ev);
    D          = sqrt(D2);
    NEFF       = -1i*diag(D);
    [NEFF,ind] = sort(NEFF,'descend');
    Hy         = Hy(:,ind);
    
    % Get Effective Index of Fundamental Mode
    neff1 = NEFF(1);
    
% ANALYZE ON WAVEGUIDE

    % Define Waveguide
    a = rib_t + rib_h;
    b = 5*lam0;
    
    % Calculate 1D Grid for Slab
    dz  = lam0/nmax/NRES;
    nz  = ceil(a/dz);
    dz  = a/nz;
    Sz  = b + a + b;
    Nz  = ceil(Sz/dz);
    Sz  = Nz*dz;
    za  = [1:Nz]*dz;
    za  = za - mean(za);
    Nz2 = 2*Nz;
    dz2 = dz/2;
    
    % Build Slab
    ER2            = ones(Nz2,1);
    UR2            = ones(Nz2,1);
    nz1            = 1 + ceil(b/dz2);
    nz2            = nz1 + round(a/dz2) - 1;
    ER2(1:nz1-1)   = rib_n1^2;
    ER2(nz1:nz2)   = rib_n2^2;
    ER2(nz2+1:Nz2) = rib_n3^2;
    
    % Build Eigen-Value Problem
    ERxx = diag(sparse(ER2(2:2:Nz2)));
    ERzz = diag(sparse(ER2(1:2:Nz2)));
    URyy = diag(sparse(UR2(2:2:Nz2)));
    
    NS  = [1 Nz];
    RES = [1 dz];
    BC  = [0 0];
    [~,DEZ,~,DHZ] = yeeder2d(NS,k0*RES,BC);

    A = -(DEZ/ERxx*DHZ + URyy);
    B = inv(ERzz);

    % Solve Eigen-Value Problem
    ev         = -rib_n2^2;
    [Hy,D2]    = eigs(A,B,4,ev);
    D          = sqrt(D2);
    NEFF       = -1i*diag(D);
    [NEFF,ind] = sort(NEFF,'descend');
    Hy         = Hy(:,ind);
    
    % Get Effective Index of Fundamental Mode
    neff2 = NEFF(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD OPTICAL INTEGRATED CIRCUIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATERIALS ARRAYS
ER2 = neff1^2*ones(Nx2,Ny2);
UR2 = ones(Nx2,Ny2);

% INPUT WAVEGUIDE
ny1 = 2*NPML(3) + round(SPACER(3)/dy2) + 1;
ny2 = ny1 + round(rib_w/dy2) - 1;
ER2(:,ny1:ny2) = neff2^2;

% OUTPUT WAVEGUIDE
nx = 2*NPML(1) + round(SPACER(1)/dx2) + round(L1/dx2) + 1;
ny1 = ny2 + round(g/dy2) - 1;
ny2 = ny1 + round(rib_w/dy2) - 1;
ER2(nx:Nx2,ny1:ny2) = neff2^2;

% SHOW UNIT CELL
imagesc(xa2,ya2,ER2.');
axis equal tight;
colorbar;
title('Lattice');
drawnow;

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
%% ANALYZE INPUT WAVEGUIDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXTRACT SLAB WAVEGUIDE FROM GRID
nx   = 2*(NPML(1) + 1);
erzz = diag(sparse(ER2(nx,1:2:Ny2)));
urxx = diag(sparse(UR2(nx,2:2:Ny2)));
uryy = diag(sparse(UR2(nx+1,1:2:Ny2)));
    
% BUILD DERIVATIVE MATRICES
NS  = [1 Ny];
RES = [1 dy];
BC  = [0 0];
[~,DEY,~,DHY] = yeeder2d(NS,k0*RES,BC);

% BUILD EIGEN-VALUE PROBLEM
A = -(DHY/urxx*DEY + erzz);
B = inv(uryy);

% SOLVE FULL EIGEN-VALUE PROBLEM
[Ez,D2]    = eig(full(A),full(B));
D          = sqrt(D2);
NEFF       = -1i*diag(D);
[NEFF,ind] = sort(real(NEFF),'descend');
Ez         = Ez(:,ind);

% GET SOURCE MODE
neff  = NEFF(1);
Ezsrc = Ez(:,1);

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
fsrc = zeros(Nx,Ny);
for nx = 1 : Nx
    fsrc(nx,:) = Ezsrc*exp(-1i*k0*neff*nx*dx);
end

% CALCULATE SCATTERED-FIELD MASKING MATRIX
nx        = NPML(1) + 2;
Q         = zeros(Nx,Ny);
Q(1:nx,:) = 1;
Q         = diag(sparse(Q(:)));

% CALCULATE SOURCE VECTOR
b = (Q*A - A*Q)*fsrc(:);

% SOLVE FOR FIELD
f = A\b;
f = reshape(f,Nx,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYZE TRANSMITTED AND REFLECTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXTRACT FIELDS
nx   = NPML(1) + 1;
fref = f(nx,:);
nx   = Nx - NPML(2);
ftrn = f(nx,:);

% CALCULATE MODE AMPLITUDES
aref = Ez\fref(:);
atrn = Ez\ftrn(:);

% CALCULATE REFLECTANCE AND TRANSMITTANCE
R = 100*abs(aref(1))^2
T = 100*abs(atrn(1))^2

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
er0 = (neff1^2 + neff2^2)/2;
contour(xa2,ya2,ER2.',er0,'Color','w','LineWidth',0.5);

% SET GRAPHICS VIEW
hold off;
axis equal tight;
set(gca,'YDir','reverse');
caxis(0.25*[-1 1]);
colorbar;
colormap('gray');