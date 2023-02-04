function [DEX,DEY,DHX,DHY] = yeeder2d(NS,RES,BC,kinc)
% YEEDER2D      Derivative Matrices on a 2D Yee Grid
%
% [DEX,DEY,DHX,DHY] = yeeder2d(NS,RES,BC,kinc);
%
% INPUT ARGUMENTS
% =================
% NS    [Nx Ny] Grid Size
% RES   [dx dy] Grid Resolution
% BC    [xbc ybc] Boundary Conditions
%         0: Dirichlet boundary conditions
%         1: Periodic boundary conditions
% kinc  [kx ky] Incident Wave Vector
%       This argument is only needed for PBCs.
%
% Note: For normalized grids use k0*RES and kinc/k0
%
% OUTPUT ARGUMENTS
% =================
% DEX   Derivative Matrix wrt x for Electric Fields
% DEY   Derivative Matrix wrt to y for Electric Fields
% DHX   Derivative Matrix wrt to x for Magnetic Fields
% DHY   Derivative Matrix wrt to y for Magnetic Fields

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HANDLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXTRACT GRID PARAMETERS
Nx = NS(1);     dx = RES(1);
Ny = NS(2);     dy = RES(2);

% DEFAULT KINC
if ~exist('kinc')
    kinc = [0 0];
end

% DETERMINE MATRIX SIZE
M = Nx*Ny;

% ZERO MATRIX
Z = sparse(M,M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HANDLE IF SIZE IS 1 CELL
if Nx==1
    DEX = -1i*kinc(1)*speye(M,M);
    
% HANDLE ALL OTHER CASES
else
    
    % Center Diagonal
    d0 = -ones(M,1);
    
    % Upper Diagonal
    d1 = ones(M,1);
    d1(Nx+1:Nx:M) = 0;
    
    % Build Derivative Matrix with Dirichlet BCs
    DEX = spdiags([d0 d1]/dx,[0 1],Z);
    
    % Incorporate Periodic Boundary Conditions
    if BC(1)==1
        d1 = zeros(M,1);
        d1(1:Nx:M) = exp(-1i*kinc(1)*Nx*dx)/dx;
        DEX = spdiags(d1,1-Nx,DEX);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HANDLE IF SIZE IS 1 CELL
if Ny==1
    DEY = -1i*kinc(2)*speye(M,M);
    
% HANDLE ALL OTHER CASES
else
    
    % Center Diagonal
    d0 = -ones(M,1);
    
    % Upper Diagonal
    d1 = ones(M,1);
    
    % Build Derivative Matrix with Dirichlet BCs
    DEY = spdiags([d0 d1]/dy,[0 Nx],Z);
    
    % Incorporate Periodic Boundary Conditions
    if BC(2)==1
        d1 = (exp(-1i*kinc(2)*Ny*dy)/dy)*ones(M,1);
        DEY = spdiags(d1,Nx-M,DEY);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DHX AND DHY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DHX = -DEX';
DHY = -DEY';