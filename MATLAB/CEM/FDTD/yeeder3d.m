function [DEX,DEY,DEZ,DHX,DHY,DHZ] = yeeder3d(NS,RES,BC,kinc)
% YEEDER3D      Derivative Matrices on a 3D Yee Grid
%
% [DEX,DEY,DEZ,DHX,DHY,DHZ] = yeeder3d(NS,RES,BC,kinc);
%
% INPUT ARGUMENTS
% =================
% NS    [Nx Ny Nz] Grid Size
% RES   [dx dy dz] Grid Resolution
% BC    [xbc ybc zbc] Boundary Conditions
%         0: Dirichlet boundary conditions
%         1: Periodic boundary conditions
% kinc  [kx ky kz] Incident Wave Vector
%       This argument is only needed for PBCs.
%
% Note: For normalized grids use k0*RES and kinc/k0
%
% OUTPUT ARGUMENTS
% =================
% DEX   Derivative Matrix wrt x for Electric Fields
% DEY   Derivative Matrix wrt y for Electric Fields
% DEZ   Derivative Matrix wrt z for Electric Fields
% DHX   Derivative Matrix wrt x for Magnetic Fields
% DHY   Derivative Matrix wrt y for Magnetic Fields
% DHZ   Derivative Matrix wrt z for Magnetic Fields

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HANDLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXTRACT GRID PARAMETERS
Nx = NS(1);     dx = RES(1);
Ny = NS(2);     dy = RES(2);
Nz = NS(3);     dz = RES(3);

% DEFAULT KINC
if ~exist('kinc')
    kinc = [0 0 0];
end

% DETERMINE MATRIX SIZE
M = Nx*Ny*Nz;

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
    d1 = [ ones((Ny-1)*Nx,1) ; zeros(Nx,1) ];
    d1 = repmat(d1,[Nz-1 1]);
    d1 = [ zeros(Nx,1) ; d1 ; ones((Ny-1)*Nx,1) ];
    
    % Build Derivative Matrix with Dirichlet BCs
    DEY = spdiags([d0 d1]/dy,[0 Nx],Z);
    
    % Incorporate Periodic Boundary Conditions
    if BC(2)==1
        ph = exp(-1i*kinc(2)*Ny*dy)/dy;
                
        d1 = [ ones(Nx,1) ; zeros((Ny-1)*Nx,1) ];
        d1 = repmat(d1,[Nz-1 1]);
        d1 = [ d1 ; ones((Ny-1)*Nx,1) ; zeros(Nx,1) ];

        DEY = spdiags(ph*d1,-Nx*(Ny-1),DEY);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DEZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HANDLE IF SIZE IS 1 CELL
if Nz==1
    DEZ = -1i*kinc(3)*speye(M,M);
    
% HANDLE ALL OTHER CASES
else
    
    % Center Diagonal
    d0 = ones(M,1);
        
    % Build Derivative Matrix with Dirichlet BCs
    DEZ = spdiags([-d0 +d0]/dz,[0 Nx*Ny],Z);
    
    % Incorporate Periodic Boundary Conditions
    if BC(3)==1
        d0  = (exp(-1i*kinc(3)*Nz*dz)/dz)*ones(M,1);
        DEZ = spdiags(d0,-Nx*Ny*(Nz-1),DEZ);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD DHX, DHY AND DHZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DHX = -DEX';
DHY = -DEY';
DHZ = -DEZ';