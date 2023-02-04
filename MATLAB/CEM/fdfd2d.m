function DAT = fdfd2d(DEV,SRC)
% FDFD2D    Two-Dimensional FDFD for Periodic Structures
%
% DAT = fdfd2d(DEV,SRC);
%
% INPUT ARGUMENTS
% ================
% DEV       Device Parameters
%   .ER2    relative permittivity on 2x grid
%   .UR2    relative permeability on 2x grid
%   .RES    [dx dy] grid resolution on Yee grid
%   .NPML   [NYLO NYHI] size of UPML on Yee grid
%
% SRC       Source Parameters
%   .lam0   free space wavelength
%   .theta  angle of incidence
%   .MODE   'E' or 'H'
%
% OUTPUT ARGUMENTS
% ================
% DAT       Output Data
%   .f      simulated field
%   .m      diffraction order indices
%   .RDE    diffraction efficiencies of reflected waves
%   .REF    overall reflectance
%   .TDE    diffraction efficiencies of transmitted waves
%   .TRN    overall transmittance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HANDLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GET SIZE OF GRID
[Nx2,Ny2] = size(DEV.ER2);

% CALCULATE GRID
Nx = Nx2/2;
Ny = Ny2/2;
dx = DEV.RES(1);
dy = DEV.RES(2);
Sx = Nx*dx;
Sy = Ny*dy;
xa = [1:Nx]*dx;
ya = [1:Ny]*dy;
[Y,X] = meshgrid(ya,xa);

% WAVE NUMBER
k0 = 2*pi/SRC.lam0;

% EXTRACT MATERIAL PROPERTIES IN EXTERNAL REGIONS
urref = DEV.UR2(1,1);
urtrn = DEV.UR2(1,Ny2);
erref = DEV.ER2(1,1);
ertrn = DEV.ER2(1,Ny2);
nref  = sqrt(urref*erref);
ntrn  = sqrt(urtrn*ertrn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FDFD ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INCORPORATE UPML
[ERxx,ERyy,ERzz,URxx,URyy,URzz] ...
    = addupml2d(DEV.ER2,DEV.UR2,[0 0 DEV.NPML]);

% DIAGONALIZE MATERIAL TENSORS
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));
URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

% INCIDENT WAVE VECTOR
kxinc = k0*nref*sin(SRC.theta);
kyinc = k0*nref*cos(SRC.theta);
kinc  = [ kxinc ; kyinc ];

% BUILD DERIVATIVE MATRICES
NS  = [Nx Ny];
BC  = [1 0];
[DEX,DEY,DHX,DHY] = yeeder2d(NS,k0*DEV.RES,BC,kinc/k0);
        
% BUILD WAVE MATRIX
if SRC.MODE == 'E'
    A = DHX/URyy*DEX + DHY/URxx*DEY + ERzz;
else
    A = DEX/ERyy*DHX + DEY/ERxx*DHY + URzz;
end

% CALCULATE SOURCE FIELD
fsrc  = exp(-1i*(kxinc*X + kyinc*Y));

% CALCULATE SCATTERED-FIELD MASKING MATRIX
ny = DEV.NPML(1) + 2;
Q  = zeros(Nx,Ny);
Q(:,1:ny) = 1;
Q  = diag(sparse(Q(:)));

% CALCULATE SOURCE VECTOR
b = (Q*A - A*Q)*fsrc(:);

% SOLVE FOR FIELD
DAT.f = A\b;
DAT.f = reshape(DAT.f,Nx,Ny);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYZE REFLECTION AND TRANSMISSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE WAVE VECTOR COMPONENTS
DAT.m = [-floor(Nx/2):+floor((Nx-1)/2)].';
kx    = kxinc - DAT.m*2*pi/Sx;
kyref = sqrt((k0*nref)^2 - kx.^2);
kytrn = sqrt((k0*ntrn)^2 - kx.^2);

% EXTRACT REFELECTED AND TRANSMITTED WAVES
fref = DAT.f(:,DEV.NPML(1)+1)./fsrc(:,1);
ftrn = DAT.f(:,Ny-DEV.NPML(2))./fsrc(:,1);

% CALCULATE AMPLITUDES OF SPATIAL HARMONICS
aref = fftshift(fft(fref))/Nx;
atrn = fftshift(fft(ftrn))/Nx;

% CALCULATE DIFFRACTION EFFICENCIES
DAT.RDE = abs(aref).^2.*real(kyref/kyinc);
if SRC.MODE == 'E'
    DAT.TDE = abs(atrn).^2.*real(urref/urtrn*kytrn/kyinc);
else
    DAT.TDE = abs(atrn).^2.*real(erref/ertrn*kytrn/kyinc);
end

% CALCULATE OVERALL REFLECTANCE & TRANSMITTANCE
DAT.REF = sum(DAT.RDE(:));
DAT.TRN = sum(DAT.TDE(:));