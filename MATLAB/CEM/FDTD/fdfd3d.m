function DAT = fdfd3d(DEV,SRC)
% FDFD3D    Three-Dimensional FDFD for Periodic Structures
%
% DAT = fdfd3d(DEV,SRC);
%
% INPUT ARGUMENTS
% ================
% DEV       Device Parameters
%
%   .ER2 or
%   .ER2xx   .ER2xy   .ER2xz
%   .ER2yx   .ER2yy   .ER2yz     relative perm. on 2x grid
%   .ER2zx   .ER2zy   .ER2zz
%
%   .UR2 or
%   .UR2xx   .UR2xy   .UR2xz
%   .UR2yx   .UR2yy   .UR2yz     relative perm. on 2x grid
%   .UR2zx   .UR2zy   .UR2zz
%
%   .RES    [dx dy] grid resolution on Yee grid
%   .NPML   [NZLO NZHI] size of SCPML on Yee grid
%
% SRC       Source Parameters
%
%   .lam0   free space wavelength
%   .theta  angle of incidence
%   .phi    angle of incidence
%   .pte    complex amplitude of the TE polarization
%   .ptm    complex amplitude of the TM polarization
%
% OUTPUT ARGUMENTS
% ================
% DAT       Output Data
%
%   .fx   .fy   .fz       simulated field
%   .m    .n              diffraction order indices
%
%   .RDE    diffraction efficiencies of reflected waves
%   .REF    overall reflectance
%   .TDE    diffraction efficiencies of transmitted waves
%   .TRN    overall transmittance
%
%   .s11    complex ref. coefficients of diffraction orders
%   .s21    complex tran. coefficients of diffraction orders

% DEFINE SOLVER PARAMETERS
tol   = 1e-6;
maxit = 15000;

% ANONYMOUS FUNCTIONS
diagonalize = @(x) diag(sparse(x(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HANDLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GET SIZE OF GRID
if isfield(DEV,'ER2')
    [Nx2,Ny2,Nz2] = size(DEV.ER2);
else
    [Nx2,Ny2,Nz2] = size(DEV.ER2xx);
end

% CALCULATE GRID
Nx = Nx2/2;
Ny = Ny2/2;
Nz = Nz2/2;
dx = DEV.RES(1);
dy = DEV.RES(2);
dz = DEV.RES(3);
Sx = Nx*dx;
Sy = Ny*dy;
Sz = Nz*dz;

% GRID AXES
xa2 = [1:Nx2]*dx/2;          xa2 = xa2 - mean(xa2);
ya2 = [1:Ny2]*dy/2;          ya2 = ya2 - mean(ya2);
za2 = [0.5:Nz2-0.5]*dz/2;
[Y2,X2,Z2] = meshgrid(ya2,xa2,za2);

% WAVE NUMBER
k0 = 2*pi/SRC.lam0;

% SPECIAL MATRICES
M  = Nx*Ny*Nz;
ZZ = sparse(M,M);
I  = speye(M,M);

% PERMITTIVITY TENSOR ELEMENTS
if isfield(DEV,'ER2')
    ERxx = diagonalize(DEV.ER2(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    ERxy = ZZ;
    ERxz = ZZ;
    ERyx = ZZ;
    ERyy = diagonalize(DEV.ER2(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    ERyz = ZZ;
    ERzx = ZZ;
    ERzy = ZZ;
    ERzz = diagonalize(DEV.ER2(1:2:Nx2,1:2:Ny2,2:2:Nz2));
else
    ERxx = diagonalize(DEV.ER2xx(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    ERxy = diagonalize(DEV.ER2xy(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    ERxz = diagonalize(DEV.ER2xz(1:2:Nx2,1:2:Ny2,2:2:Nz2));
    ERyx = diagonalize(DEV.ER2yx(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    ERyy = diagonalize(DEV.ER2yy(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    ERyz = diagonalize(DEV.ER2yz(1:2:Nx2,1:2:Ny2,2:2:Nz2));
    ERzx = diagonalize(DEV.ER2zx(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    ERzy = diagonalize(DEV.ER2zy(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    ERzz = diagonalize(DEV.ER2zz(1:2:Nx2,1:2:Ny2,2:2:Nz2));
end

% PERMEABILITY TENSOR ELEMENTS
if isfield(DEV,'UR2')
    URxx = diagonalize(DEV.UR2(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    URxy = ZZ;
    URxz = ZZ;
    URyx = ZZ;
    URyy = diagonalize(DEV.UR2(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    URyz = ZZ;
    URzx = ZZ;
    URzy = ZZ;
    URzz = diagonalize(DEV.UR2(1:2:Nx2,1:2:Ny2,2:2:Nz2));
else
    URxx = diagonalize(DEV.UR2xx(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    URxy = diagonalize(DEV.UR2xy(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    URxz = diagonalize(DEV.UR2xz(1:2:Nx2,1:2:Ny2,2:2:Nz2));
    URyx = diagonalize(DEV.UR2yx(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    URyy = diagonalize(DEV.UR2yy(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    URyz = diagonalize(DEV.UR2yz(1:2:Nx2,1:2:Ny2,2:2:Nz2));
    URzx = diagonalize(DEV.UR2zx(2:2:Nx2,1:2:Ny2,1:2:Nz2));
    URzy = diagonalize(DEV.UR2zy(1:2:Nx2,2:2:Ny2,1:2:Nz2));
    URzz = diagonalize(DEV.UR2zz(1:2:Nx2,1:2:Ny2,2:2:Nz2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FDFD ANALYS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXTRACT MATERIAL PROPERTIES IN EXTERNAL REGIONS
ur_ref = full(URxx(1,1));
ur_trn = full(URxx(M,M));
er_ref = full(ERxx(1,1));
er_trn = full(ERxx(M,M));
n_ref  = sqrt(ur_ref*er_ref);
n_trn  = sqrt(ur_trn*er_trn);

% CALCULATE PML PARAMETERS
[sx2,sy2,sz2] = calcpml3d([Nx2 Ny2 Nz2],2*[0 0 0 0 DEV.NPML]);

sx_ey = diagonalize(1./sx2(1:2:Nx2,2:2:Ny2,1:2:Nz2));
sx_ez = diagonalize(1./sx2(1:2:Nx2,1:2:Ny2,2:2:Nz2));
sy_ex = diagonalize(1./sy2(2:2:Nx2,1:2:Ny2,1:2:Nz2));
sy_ez = diagonalize(1./sy2(1:2:Nx2,1:2:Ny2,2:2:Nz2));
sz_ex = diagonalize(1./sz2(2:2:Nx2,1:2:Ny2,1:2:Nz2));
sz_ey = diagonalize(1./sz2(1:2:Nx2,2:2:Ny2,1:2:Nz2));

sx_hy = diagonalize(1./sx2(2:2:Nx2,1:2:Ny2,2:2:Nz2));
sx_hz = diagonalize(1./sx2(2:2:Nx2,2:2:Ny2,1:2:Nz2));
sy_hx = diagonalize(1./sy2(1:2:Nx2,2:2:Ny2,2:2:Nz2));
sy_hz = diagonalize(1./sy2(2:2:Nx2,2:2:Ny2,1:2:Nz2));
sz_hx = diagonalize(1./sz2(1:2:Nx2,2:2:Ny2,2:2:Nz2));
sz_hy = diagonalize(1./sz2(2:2:Nx2,1:2:Ny2,2:2:Nz2));

% CALCULATE WAVE VECTORS
kinc   = k0*n_ref*[ sin(SRC.theta)*cos(SRC.phi) ; ...
                    sin(SRC.theta)*sin(SRC.phi)  ; ...
                    cos(SRC.theta) ];
m      = [-floor(Nx/2):+floor((Nx-1)/2)].';
n      = [-floor(Ny/2):+floor((Ny-1)/2)].';
kx     = kinc(1) - m*2*pi/Sx;
ky     = kinc(2) - n*2*pi/Sy;
kz_ref = sqrt((k0*n_ref)^2 - kx.^2 - ky.^2);
kz_trn = sqrt((k0*n_trn)^2 - kx.^2 - ky.^2);

% BUILD DERIVATIVE MATRICES
NS  = [Nx Ny Nz];
RES = [dx dy dz];
BC  = [1 1 0];
[DEX,DEY,DEZ,DHX,DHY,DHZ] = yeeder3d(NS,k0*RES,BC,kinc/k0);

% CALCULATE INTERPOLATION MATRICES
RX = (0.5*k0*dx)*DEX;
RY = (0.5*k0*dy)*DEY;
RZ = (0.5*k0*dz)*DEZ;

% FORM MATERIALS TENSORS
ER = [ ERxx        RX*RY'*ERxy RX*RZ'*ERxz ; ...
       RY*RX'*ERyx ERyy        RY*RZ'*ERyz ; ...
       RZ*RX'*ERzx RZ*RY'*ERzy ERzz ];
UR = [ URxx        RX*RY'*URxy RX*RZ'*URxz ; ...
       RY*RX'*URyx URyy        RY*RZ'*URyz ; ...
       RZ*RX'*URzx RZ*RY'*URzy URzz ];

% BUILD THE WAVE MATRIX A
CE = [ ZZ , -sz_hx*DEZ , sy_hx*DEY ; ...
       sz_hy*DEZ , ZZ , -sx_hy*DEX ; ...
       -sy_hz*DEY , sx_hz*DEX , ZZ ];
CH = [ ZZ , -sz_ex*DHZ , sy_ex*DHY ; ...
       sz_ey*DHZ , ZZ , -sx_ey*DHX ; ...
       -sy_ez*DHY , sx_ez*DHX , ZZ ];
A  = CH/UR*CE - ER;

% CALCULATE POLARIZATION VECTOR P
az = [0;0;1];
if abs(SRC.theta)<1e-6
    ate = [0;1;0];
else
    ate = cross(az,kinc);
    ate = ate/norm(ate);
end
atm = cross(kinc,ate);
atm = atm/norm(atm);
P   = SRC.pte*ate + SRC.ptm*atm;

% CALCULATE SOURCE FIELD fsrc
fphase = exp(-1i*(kinc(1)*X2 + kinc(2)*Y2 + kinc(3)*Z2));
fx     = P(1)*fphase(2:2:Nx2,1:2:Ny2,1:2:Nz2);
fy     = P(2)*fphase(1:2:Nx2,2:2:Ny2,1:2:Nz2);
fz     = P(3)*fphase(1:2:Nx2,1:2:Ny2,2:2:Nz2);
fsrc   = [ fx(:) ; fy(:) ; fz(:) ];

% CALCULATE SCATERED FIELD MASKING MATRIX
nz = DEV.NPML(1) + 3;
Q  = zeros(Nx,Ny,Nz);
Q(:,:,1:nz) = 1;
Q = diagonalize(Q);
Q = [ Q ZZ ZZ ; ZZ Q ZZ ; ZZ ZZ Q ];

% CALCULATE SOURCE VECTOR B
b = (Q*A - A*Q)*fsrc;

% SOLVE FOR FIELD
f = zeros(3*M,1);
A = gpuArray(A);
b = gpuArray(b);
f = bicg(A,b,tol,maxit,[],[],f);
f = gather(f);

% EXTRACT FIELD COMPONENTS
DAT.fx = reshape(f(    1:  M),Nx,Ny,Nz);
DAT.fy = reshape(f(  M+1:2*M),Nx,Ny,Nz);
DAT.fz = reshape(f(2*M+1:3*M),Nx,Ny,Nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANALYZE REFLECTION AND TRANSMISSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% EXTRACT FIELD FROM RECORD PLANES
nz_ref = DEV.NPML(1) + 2;
fx_ref = DAT.fx(:,:,nz_ref);
fy_ref = DAT.fy(:,:,nz_ref);
fz_ref = DAT.fz(:,:,nz_ref-1:nz_ref);

nz_trn = Nz - DEV.NPML(1);
fx_trn = DAT.fx(:,:,nz_trn);
fy_trn = DAT.fy(:,:,nz_trn);
fz_trn = DAT.fz(:,:,nz_trn-1:nz_trn);

% INTERPOLATE FIELDS TO ORIGIN OF YEE CELLS
Ex_ref = fx_ref;
Ex_ref(1,:)    = (fx_ref(Nx,:)     + fx_ref(1,:))/2;
Ex_ref(2:Nx,:) = (fx_ref(1:Nx-1,:) + fx_ref(2:Nx,:))/2;

Ey_ref = fy_ref;
Ey_ref(:,1)    = (fy_ref(:,Ny)     + fy_ref(:,1))/2;
Ey_ref(:,2:Ny) = (fy_ref(:,1:Ny-1) + fy_ref(:,2:Ny))/2;

Ez_ref = (fz_ref(:,:,1) + fz_ref(:,:,2))/2;

Ex_trn = fx_trn;
Ex_trn(1,:)    = (fx_trn(Nx,:)     + fx_trn(1,:))/2;
Ex_trn(2:Nx,:) = (fx_trn(1:Nx-1,:) + fx_trn(2:Nx,:))/2;

Ey_trn = fy_trn;
Ey_trn(:,1,:)    = (fy_trn(:,Ny,:)     + fy_trn(:,1,:))/2;
Ey_trn(:,2:Ny) = (fy_trn(:,1:Ny-1) + fy_trn(:,2:Ny))/2;

Ez_trn = (fz_trn(:,:,1) + fz_trn(:,:,2))/2;

% REMOVE PHASE TILT
Ex_ref = Ex_ref./fphase(2:2:Nx2,1:2:Ny2,1);
Ey_ref = Ey_ref./fphase(1:2:Nx2,2:2:Ny2,1);
Ez_ref = Ez_ref./fphase(1:2:Nx2,1:2:Ny2,2);

Ex_trn = Ex_trn./fphase(2:2:Nx2,1:2:Ny2,1);
Ey_trn = Ey_trn./fphase(1:2:Nx2,2:2:Ny2,1);
Ez_trn = Ez_trn./fphase(1:2:Nx2,1:2:Ny2,2);

% CALCULATE AMPLITUDES OF DIFFRACTION ORDERS
ax_ref = fftshift(fft2(Ex_ref))/(Nx*Ny);
ay_ref = fftshift(fft2(Ey_ref))/(Nx*Ny);
az_ref = fftshift(fft2(Ez_ref))/(Nx*Ny);
a_ref  = abs(ax_ref).^2 + abs(ay_ref).^2 + abs(az_ref).^2;
DAT.s11 = sqrt(ax_ref.^2 + ay_ref.^2 + az_ref.^2);

ax_trn = fftshift(fft2(Ex_trn))/(Nx*Ny);
ay_trn = fftshift(fft2(Ey_trn))/(Nx*Ny);
az_trn = fftshift(fft2(Ez_trn))/(Nx*Ny);
a_trn  = abs(ax_trn).^2 + abs(ay_trn).^2 + abs(az_trn).^2;
DAT.s21 = sqrt(ax_trn.^2 + ay_trn.^2 + az_trn.^2);

% CALCULATE DIFFRACTION EFFICIENCIES
DAT.RDE = a_ref.*real(kz_ref/kinc(3));
DAT.TDE = a_trn.*real(ur_ref/ur_trn*kz_trn/kinc(3));

% Calculate Overall Response
DAT.REF = sum(DAT.RDE(:));
DAT.TRN = sum(DAT.TDE(:));
DAT.CON = DAT.REF + DAT.TRN;
