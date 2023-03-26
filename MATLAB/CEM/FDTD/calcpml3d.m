function [sx,sy,sz] = calcpml3d(NGRID,NPML)
% CALCPML3D     Calculate PML Parameters
%
% [sx,sy,sz] = calcpml3d(NGRID,NPML);
%
% INPUT ARGUMENTS
% ================
% NGRID     [Nx Ny Nz] Number of Cells on Grid
% NPML      [NXLO NXHI NYLO NYHI NZLO NZHI] Size of PML

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HANDLE INPUT AND OUTPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE PML PARAMETERS
amax = 4;
cmax = 1;
p    = 3;

% EXTRACT GRID SIZE
Nx = NGRID(1);
Ny = NGRID(2);
Nz = NGRID(3);

% EXTRACT PML SIZE
NXLO = NPML(1);
NXHI = NPML(2);
NYLO = NPML(3);
NYHI = NPML(4);
NZLO = NPML(5);
NZHI = NPML(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE PML PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE PML PARAMETERS TO PROBLEM SPACE
sx = ones(Nx,Ny,Nz);
sy = ones(Nx,Ny,Nz);
sz = ones(Nx,Ny,Nz);

% CALCULATE XLO PML
for nx = 1 : NXLO
    ax = 1 + (amax - 1)*(nx/NXLO)^p;
    cx = cmax*sin(0.5*pi*nx/NXLO)^2;
    sx(NXLO - nx + 1,:,:) = ax*(1 - 1i*60*cx);
end

% CALCULATE XHI PML
for nx = 1 : NXHI
    ax = 1 + (amax - 1)*(nx/NXHI)^p;
    cx = cmax*sin(0.5*pi*nx/NXHI)^2;
    sx(Nx - NXHI + nx,:,:) = ax*(1 - 1i*60*cx);
end

% CALCULATE YLO PML
for ny = 1 : NYLO
    ay = 1 + (amax - 1)*(ny/NYLO)^p;
    cy = cmax*sin(0.5*pi*ny/NYLO)^2;
    sy(:,NYLO - ny + 1,:) = ay*(1 - 1i*60*cy);
end

% CALCULATE YHI PML
for ny = 1 : NYHI
    ay = 1 + (amax - 1)*(ny/NYHI)^p;
    cy = cmax*sin(0.5*pi*ny/NYHI)^2;
    sy(:,Ny - NYHI + ny,:) = ay*(1 - 1i*60*cy);
end
% CALCULATE ZLO PML
for nz = 1 : NZLO
    az = 1 + (amax - 1)*(nz/NZLO)^p;
    cz = cmax*sin(0.5*pi*nz/NZLO)^2;
    sz(:,:,NZLO - nz + 1) = az*(1 - 1i*60*cz);
end

% CALCULATE ZHI PML
for nz = 1 : NZHI
    az = 1 + (amax - 1)*(nz/NZHI)^p;
    cz = cmax*sin(0.5*pi*nz/NZHI)^2;
    sz(:,:,Nz - NZHI + nz) = az*(1 - 1i*60*cz);
end