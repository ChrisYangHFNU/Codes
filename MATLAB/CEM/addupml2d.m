function [ERxx,ERyy,ERzz,URxx,URyy,URzz] = addupml2d(ER2,UR2,NPML)
% ADDUPML2D     Add UPML to a 2D Yee Grid
%
% [ERxx,ERyy,ERzz,URxx,URyy,URzz] = addupml2d(ER2,UR2,NPML);
%
% INPUT ARGUMENTS
% ================
% ER2       Relative Permittivity on 2x Grid
% UR2       Relative Permeability on 2x Grid
% NPML      [NXLO NXHI NYLO NYHI] Size of UPML on 1x Grid
%
% OUTPUT ARGUMENTS
% ================
% ERxx      xx Tensor Element for Relative Permittivity
% ERyy      yy Tensor Element for Relative Permittivity
% ERzz      zz Tensor Element for Relative Permittivity
% URxx      xx Tensor Element for Relative Permeability
% URyy      yy Tensor Element for Relative Permeability
% URzz      zz Tensor Element for Relative Permeability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE PML PARAMETERS
amax = 4;
cmax = 1;
p    = 3;

% EXTRACT GRID PARAMETERS
[Nx2,Ny2] = size(ER2);

% EXTRACT PML PARAMETERS
NXLO = 2*NPML(1);
NXHI = 2*NPML(2);
NYLO = 2*NPML(3);
NYHI = 2*NPML(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE PML PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE PML PARAMETERS TO PROBLEM SPACE
sx = ones(Nx2,Ny2);
sy = ones(Nx2,Ny2);

% ADD XLO PML
for nx = 1 : NXLO
    ax = 1 + (amax - 1)*(nx/NXLO)^p;
    cx = cmax*sin(0.5*pi*nx/NXLO)^2;
    sx(NXLO - nx + 1,:) = ax*(1 - 1i*60*cx);
end

% ADD XHI PML
for nx = 1 : NXHI
    ax = 1 + (amax - 1)*(nx/NXHI)^p;
    cx = cmax*sin(0.5*pi*nx/NXHI)^2;
    sx(Nx2 - NXHI + nx,:) = ax*(1 - 1i*60*cx);
end

% ADD YLO PML
for ny = 1 : NYLO
    ay = 1 + (amax - 1)*(ny/NYLO)^p;
    cy = cmax*sin(0.5*pi*ny/NYLO)^2;
    sy(:,NYLO - ny + 1) = ay*(1 - 1i*60*cy);
end

% ADD YHI PML
for ny = 1 : NYHI
    ay = 1 + (amax - 1)*(ny/NYHI)^p;
    cy = cmax*sin(0.5*pi*ny/NYHI)^2;
    sy(:,Ny2 - NYHI + ny) = ay*(1 - 1i*60*cy);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INCORPORATE PML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE TENSOR ELEMENTS WITH UPML
ERxx = ER2./sx.*sy;
ERyy = ER2.*sx./sy;
ERzz = ER2.*sx.*sy;

URxx = UR2./sx.*sy;
URyy = UR2.*sx./sy;
URzz = UR2.*sx.*sy;

% EXTRACT TENSOR ELEMENTS ON YEE GRID
ERxx = ERxx(2:2:Nx2,1:2:Ny2);
ERyy = ERyy(1:2:Nx2,2:2:Ny2);
ERzz = ERzz(1:2:Nx2,1:2:Ny2);

URxx = URxx(1:2:Nx2,2:2:Ny2);
URyy = URyy(2:2:Nx2,1:2:Ny2);
URzz = URzz(2:2:Nx2,2:2:Ny2);

