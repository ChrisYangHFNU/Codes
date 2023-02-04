% Chapter6_ribwaveguide.m

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

% WAVELENGTH
lam0 = 1.55 * micrometers;

% SOI RIB WAVEGUIDE PARAMETERS
rib_n1 = 1.0;
rib_n2 = 3.5;
rib_n3 = 1.5;
rib_h  = 0.60 * micrometers;
rib_t  = 0.60 * micrometers;
rib_w  = 0.80 * micrometers;

% GRID PARAMETERS
nmax   = max([rib_n1 rib_n2 rib_n3]);
NRES   = 20;
SPACER = lam0*[1.5 1.5 0.5 0.5];

% NUMBER OF MODES TO CALCULATE
NMODES = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIRST GUESS AT GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(rib_w/dx);
dx = rib_w/nx;
ny = ceil(rib_t/dy);
dy = rib_t/ny;

% GRID SIZE
Sx = SPACER(1) + rib_w + SPACER(2);
Nx = ceil(Sx/dx);
Sx = Nx*dx;

Sy = SPACER(3) + rib_h + rib_t + SPACER(4);
Ny = ceil(Sy/dy);
Sy = Ny*dy;

% 2X GRID
Nx2 = 2*Nx;         dx2 = dx/2;
Ny2 = 2*Ny;         dy2 = dy/2;

% GRID AXES
xa = [1:Nx]*dx;     xa = xa - mean(xa);
ya = [1:Ny]*dy;     ya = ya - mean(ya);
xa2 = [1:Nx2]*dx2;     xa2 = xa2 - mean(xa2);
ya2 = [1:Ny2]*dy2;     ya2 = ya2 - mean(ya2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO AIR
ER2 = rib_n1^2*ones(Nx2,Ny2);
UR2 = ones(Nx2,Ny2);

% CALCULATE POSITION INDICES
nx1 = 1 + round(SPACER(1)/dx2);
nx2 = nx1 + round(rib_w/dx2) - 1;

ny1 = 1 + round(SPACER(3)/dy2);
ny2 = ny1 + round(rib_h/dy2) - 1;
ny3 = ny2 + 1;
ny4 = ny3 + round(rib_t/dy2) - 1;
ny5 = ny4 + 1;

% BUILD RIB WAVEGUIDE
ER2(nx1:nx2,ny1:ny2) = rib_n2^2;
ER2(:,ny3:ny4)       = rib_n2^2;
ER2(:,ny5:Ny2)       = rib_n3^2;

% EXTRACT MATERIALS TENSORS
ERxx = ER2(2:2:Nx2,1:2:Ny2);
ERyy = ER2(1:2:Nx2,2:2:Ny2);
ERzz = ER2(1:2:Nx2,1:2:Ny2);

URxx = UR2(1:2:Nx2,2:2:Ny2);
URyy = UR2(2:2:Nx2,1:2:Ny2);
URzz = UR2(2:2:Nx2,2:2:Ny2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FINITE-DIFFERENCE ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DIAGONALIZE MATERIAL ARRAYS
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));
URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

% BUILD DERIVATE MATRICES
k0  = 2*pi/lam0;
NS  = [Nx Ny];
RES = [dx dy];
BC  = [0 0];
[DEX,DEY,DHX,DHY] = yeeder2d(NS,k0*RES,BC);
 
% BUILD P AND Q MATRICES
P = [ DEX/ERzz*DHY , -(DEX/ERzz*DHX + URyy) ; ...
      DEY/ERzz*DHY + URxx , -DEY/ERzz*DHX ];
Q = [ DHX/URzz*DEY , -(DHX/URzz*DEX + ERyy) ; ...
      DHY/URzz*DEY + ERxx , -DHY/URzz*DEX ];

% SOLVE EIGEN-VALUE PROBLEM
ev       = -rib_n2^2;
[Exy,D2] = eigs(P*Q,NMODES,ev);
D        = sqrt(D2);
NEFF     = -1i*diag(D);

% CALCULATE OTHER FIELD COMPONENTS
M  = Nx*Ny;
Ex = Exy(1:M,:);
Ey = Exy(M+1:2*M,:);

Hxy = (Q*Exy)/D;
Hx  = Hxy(1:M,:);
Hy  = Hxy(M+1:2*M,:);

Ez = ERzz\(DHX*Hy - DHY*Hx);
Hz = URzz\(DEX*Ey - DEY*Ex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUALIZE THE MODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VISUALIZE CALCULATED MODES
for m = 1 : NMODES
    % Extract Fields
    ex = reshape(Ex(:,m),Nx,Ny);
    ey = reshape(Ey(:,m),Nx,Ny);
    
    % Scale Fields
    fmax = max(abs([ ex(:) ; ey(:) ]));
    ex   = ex/fmax;
    ey   = ey/fmax;
    
    % Plot Modes
    subplot(NMODES,2,(m - 1)*2 + 1);
    imagesc(xa,ya,abs(ex).');
    axis equal tight off;
    colorbar;
    caxis([0 1]);
    title(['Ex, n_{eff} = ' num2str(NEFF(m))]);
    
    subplot(NMODES,2,(m - 1)*2 + 2);
    imagesc(xa,ya,abs(ey).');
    axis equal tight off;
    colorbar;
    caxis([0 1]);
    title('Ey');
            
end