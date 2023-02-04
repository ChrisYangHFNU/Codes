% Chapter6_microstrip.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
meters      = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
micrometers = 1e-6 * meters;
seconds     = 1;
hertz       = 1/seconds;
kilohertz   = 1e3*hertz;
megahertz   = 1e6*hertz;
gigahertz   = 1e9*hertz;

% CONSTANTS
c0 = 299792458 * meters/seconds;
N0 = 376.73031346177;
e0 = 8.8541878176e-12 * 1/meters;
u0 = 1.2566370614e-6 * 1/meters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE
freq = 1.0 * gigahertz;
lam0 = c0/freq;

% MATERIAL PROPERTIES
sigma = 5.8e7;
erm   = 1 + sigma/(1i*2*pi*freq*e0);

er    = 4.4;
tand  = 0.03;
erd   = er*(1 - 1i*tand);

era   = 1.0;

% MICROSTRIP PARAMETERS
w = 1 * millimeters;
h = 0.5 * millimeters;
t = 35 * micrometers;

% GRID PARAMETERS
nmax   = sqrt(real(erd));
NRES   = 40;
NDIM   = 10;
SPACER = w*[3 3 2 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COMPUTE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID RESOLUTION
dx = lam0/nmax/NRES;
dy = lam0/nmax/NRES;

% RESOLVE MINIMUM DIMENSIONS
nx = round(w/dx);
if nx<NDIM
    nx = NDIM;
    dx = w/nx;
end

ny = round(t/dy);
if ny<NDIM
    ny = NDIM;
    dy = t/ny;
end

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(w/dx);
dx = w/nx;

ny = ceil(h/dy);
dy = h/ny;

% COMPUTE GRID SIZE
Sx = SPACER(1) + w + SPACER(2);
Nx = ceil(Sx/dx);
Sx = Nx*dx;

Sy = SPACER(3) + t + h + SPACER(4);
Ny = ceil(Sy/dy) + 1;
Sy = Ny*dy;

% GRID AXES
xa = [1:Nx]*dx;         xa = xa - mean(xa);
ya = [0:Ny-1]*dy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEVICE ON GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE MATERIALS TO FREE SPACE
URxx = ones(Nx,Ny);
URyy = ones(Nx,Ny);
URzz = ones(Nx,Ny);
ERxx = era*ones(Nx,Ny);
ERyy = era*ones(Nx,Ny);
ERzz = era*ones(Nx,Ny);

% COMPUTE POSITION INDICES
nx  = round(w/dx);
nx1 = 1 + floor((Nx - nx)/2);
nx2 = nx1 + nx - 1;

ny6 = Ny;
ny5 = ny6 - 2;
ny4 = ny5 - 1;
ny3 = ny4 - round(h/dy) + 1;
ny2 = ny3 - 1;
ny1 = ny2 - round(t/dy) + 1;

% ADD DIELECTRIC
ERxx(:,ny3:Ny) = erd;
ERyy(:,ny3:Ny) = erd;
ERzz(:,ny3:Ny) = erd;

% ADD METAL
ERxx(nx1:nx2-1,ny1:ny2  ) = erm;
ERyy(nx1:nx2  ,ny1:ny2-1) = erm;
ERzz(nx1:nx2  ,ny1:ny2  ) = erm;
ERxx(:,ny5:ny6) = erm;
ERyy(:,ny5:ny6) = erm;
ERzz(:,ny5:ny6) = erm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FINITE-DIFFERENCE ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FORM DIAGONAL MATERIALS MATRICES
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));

URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

% BUILD DERIVATIVE MATRICES
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

% SOLVING EIGEN-VALUE PROBLEM
ereff    = (erd+1)/2 + (erd-1)/(2*sqrt(1+12*(h/w)));
[Exy,D2] = eigs(P*Q,1,-real(ereff));
gamman   = sqrt(D2);

% CALCULATE E AND H FIELD COMPONENTS
M   = Nx*Ny;
Ex  = Exy(1:M,:);
Ey  = Exy(M+1:2*M,:);

Hxy = (Q*Exy)/(-1i*N0*gamman);
Hx  = Hxy(1:M,:);
Hy  = Hxy(M+1:2*M,:);

% RESHAPE FIELDS BACK TO 2D GRID
Ex = reshape(Ex,Nx,Ny);
Ey = reshape(Ey,Nx,Ny);
Hx = reshape(Hx,Nx,Ny);
Hy = reshape(Hy,Nx,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE TRANSMISSION LINE PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE LINE VOLTAGE
nxc = round(Nx/2);
V0  = sum(Ey(nxc,ny2:ny4))*dy;

% CALCULATE LINE CURRENT
s  = 1;
I1 = sum(Hx(nx1-s:nx2+s,ny1-s))*dx;
I2 = sum(Hy(nx2+s        ,ny1-s+1:ny2+s))*dy;
I3 = sum(Hx(nx1-s:nx2+s,ny2+s))*dx;
I4 = sum(Hy(nx1-s        ,ny1-s+1:ny2+s))*dy;
I0 = I1 + I2 - I3 - I4;

% CALCULATE Z0 and gamma
Z0    = V0/I0;
gamma = k0*gamman;
neff  = imag(gamman);

% CALCUALTE DISTRIBUTED PARAMETERS
R = real(gamma*Z0);
L = imag(gamma*Z0)/(2*pi*freq);
G = real(gamma/Z0);
C = imag(gamma/Z0)/(2*pi*freq);

% REPORT RESULTS
disp(['gamma = ' num2str(gamma) ' 1/m']);
disp(['Z0    = ' num2str(Z0) ' ohms']);
disp(['R     = ' num2str(R/1e-3) ' mOhm/m']);
disp(['L     = ' num2str(L/1e-9) ' nH/m']);
disp(['G     = ' num2str(G/1e-6) ' uS/m']);
disp(['C     = ' num2str(C/1e-12) ' pF/m']);
disp(['neff  = ' num2str(neff)]);
disp(['ereff = ' num2str(neff^2)]);

% CALCULATE ELECTRIC FIELD MAGNITUDE
E = sqrt(abs(Ex).^2+abs(Ey).^2);
E = E/max(E(:));

% PLOT THE FIELD
pcolor(xa/millimeters,ya/millimeters,E.');
shading interp;
axis equal tight;
title(['Magnitude of Electric Field ' ...
      '$\left| \vec{E} \right|$'],...
      'Interpreter','LaTex');
set(gca,'YDir','reverse');
colormap(flipud(colormap('gray')));
xlabel('$x$ (mm)','Interpreter','LaTex');
ylabel('$y$ (mm)','Interpreter','LaTex',...
       'HorizontalAlignment','right');