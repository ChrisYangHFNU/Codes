% Chapter10_cloak.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS
degrees = pi/180;

% ANONYMOUS FUNCTIONS
diagonalize = @(x) diag(sparse(x(:)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOURCE
lam0  = 1;
theta = 30*degrees;
P     = [0;0;1];

% CLOAK
R1 = 1.0*lam0;
R2 = 3.5*lam0;

% GRID
Sx   = 10*lam0;
Sy   = Sx;
NRES = 60;
NPML = [20 20 20 20];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RESOLUTION
dx = lam0/NRES;
dy = lam0/NRES;

% GRID SIZE
Nx = ceil(Sx/dx);
dx = Sx/Nx;
Ny = ceil(Sy/dy);
dy = Sy/Ny;

% 2X GRID
Nx2 = 2*Nx;         dx2 = dx/2;
Ny2 = 2*Ny;         dy2 = dy/2;

% GRID AXES
xa = [1:Nx]*dx;     xa = xa - mean(xa);
ya = [1:Ny]*dy;     ya = ya - mean(ya);
xa2 = [1:Nx2]*dx2;  xa2 = xa2 - mean(xa2);
ya2 = [1:Ny2]*dy2;  ya2 = ya2 - mean(ya2);

% MESHGRID
[Y2,X2] = meshgrid(ya2,xa2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD CLOAK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO FREE SPACE
ERxx = ones(Nx2,Ny2);
ERxy = zeros(Nx2,Ny2);
ERxz = zeros(Nx2,Ny2);
ERyx = zeros(Nx2,Ny2);
ERyy = ones(Nx2,Ny2);
ERyz = zeros(Nx2,Ny2);
ERzx = zeros(Nx2,Ny2);
ERzy = zeros(Nx2,Ny2);
ERzz = ones(Nx2,Ny2);

% CALCULATE PERMITTIVITY TENSORS
for ny = 1 : Ny2
    for nx = 1 : Nx2
        % Calculate Polar Coordinates
        r   = sqrt(xa2(nx)^2 + ya2(ny)^2);
        phi = atan2(ya2(ny),xa2(nx));
        
        % Populate Tensor
        if r>=R1 && r<=R2
            % Calculate Polar Tensor
            ERr = (r - R1)/r;
            ERp = r/(r - R1);
            ERz = (R2/(R2-R1))^2*ERr;
            ER  = [ ERr 0 0 ; 0 ERp 0 ; 0 0 ERz ];
            % Calculate Cartesian Tensor
            R  = [ cos(phi) -sin(phi) 0 ; ...
                   sin(phi)  cos(phi) 0 ; ...
                   0 0 1 ];
            ER = R*ER/R;
            % Populate Grid
            ERxx(nx,ny) = ER(1,1);
            ERxy(nx,ny) = ER(1,2);
            ERxz(nx,ny) = ER(1,3);
            ERyx(nx,ny) = ER(2,1);
            ERyy(nx,ny) = ER(2,2);
            ERyz(nx,ny) = ER(2,3);
            ERzx(nx,ny) = ER(3,1);
            ERzy(nx,ny) = ER(3,2);
            ERzz(nx,ny) = ER(3,3);
        end
    end
end

% CALCULATE IMPERMEABILITY
D = ERxx.*ERyy.*ERzz - ERxx.*ERyz.*ERzy ...
  - ERxy.*ERyx.*ERzz + ERxy.*ERyz.*ERzx ...
  + ERxz.*ERyx.*ERzy - ERxz.*ERyy.*ERzx;
YRxx = (ERyy.*ERzz - ERyz.*ERzy)./D;
YRxy = (ERxz.*ERzy - ERxy.*ERzz)./D;
YRxz = (ERxy.*ERyz - ERxz.*ERyy)./D;
YRyx = (ERyz.*ERzx - ERyx.*ERzz)./D;
YRyy = (ERxx.*ERzz - ERxz.*ERzx)./D;
YRyz = (ERxz.*ERyx - ERxx.*ERyz)./D;
YRzx = (ERyx.*ERzy - ERyy.*ERzx)./D;
YRzy = (ERxy.*ERzx - ERxx.*ERzy)./D;
YRzz = (ERxx.*ERyy - ERxy.*ERyx)./D;

% FORM DIAGONAL MATERIALS MATRICES
ERxx = diagonalize(ERxx(2:2:Nx2,1:2:Ny2));
ERxy = diagonalize(ERxy(1:2:Nx2,2:2:Ny2));
ERxz = diagonalize(ERxz(1:2:Nx2,1:2:Ny2));
ERyx = diagonalize(ERyx(2:2:Nx2,1:2:Ny2));
ERyy = diagonalize(ERyy(1:2:Nx2,2:2:Ny2));
ERyz = diagonalize(ERyz(1:2:Nx2,1:2:Ny2));
ERzx = diagonalize(ERzx(2:2:Nx2,1:2:Ny2));
ERzy = diagonalize(ERzy(1:2:Nx2,2:2:Ny2));
ERzz = diagonalize(ERzz(1:2:Nx2,1:2:Ny2));

YRxx = diagonalize(YRxx(1:2:Nx2,2:2:Ny2));
YRxy = diagonalize(YRxy(2:2:Nx2,1:2:Ny2));
YRxz = diagonalize(YRxz(2:2:Nx2,2:2:Ny2));
YRyx = diagonalize(YRyx(1:2:Nx2,2:2:Ny2));
YRyy = diagonalize(YRyy(2:2:Nx2,1:2:Ny2));
YRyz = diagonalize(YRyz(2:2:Nx2,2:2:Ny2));
YRzx = diagonalize(YRzx(1:2:Nx2,2:2:Ny2));
YRzy = diagonalize(YRzy(2:2:Nx2,1:2:Ny2));
YRzz = diagonalize(YRzz(2:2:Nx2,2:2:Ny2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FDFD SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULATE SOURCE
k0    = 2*pi/lam0;
ersrc = ERxx(1,1);
nsrc  = sqrt(ersrc);
kinc  = k0*nsrc*[sin(theta);cos(theta);0];

E     = exp(-1i*(kinc(1)*X2 + kinc(2)*Y2));
Ex    = P(1)*E(2:2:Nx2,1:2:Ny2);
Ey    = P(2)*E(1:2:Nx2,2:2:Ny2);
Ez    = P(3)*E(1:2:Nx2,1:2:Ny2);
fsrc  = [ Ex(:) ; Ey(:) ; Ez(:) ];

% BUILD DERIVATIVE MATRICES
NS  = [Nx Ny 1];
RES = [dx dy 1];
BC  = [0 0 0];
[DEX,DEY,DEZ,DHX,DHY,DHZ] = yeeder3d(NS,k0*RES,BC,kinc/k0);

% COMPUTE INTEPOLATION MATRICES
RX = (k0*dx)*abs(DEX)/2;
RY = (k0*dy)*abs(DEY)/2;
RZ = (k0*1)*abs(DEZ)/2;

% CALCULATE SCPML MATRICES
[sx2,sy2,sz2] = calcpml3d([Nx2 Ny2 2],2*[NPML 0 0]);

sx_ey = diagonalize(1./sx2(1:2:Nx2,2:2:Ny2,1));
sx_ez = diagonalize(1./sx2(1:2:Nx2,1:2:Ny2,2));
sy_ex = diagonalize(1./sy2(2:2:Nx2,1:2:Ny2,1));
sy_ez = diagonalize(1./sy2(1:2:Nx2,1:2:Ny2,2));
sz_ex = diagonalize(1./sz2(2:2:Nx2,1:2:Ny2,1));
sz_ey = diagonalize(1./sz2(1:2:Nx2,2:2:Ny2,1));

sx_hy = diagonalize(1./sx2(2:2:Nx2,1:2:Ny2,2));
sx_hz = diagonalize(1./sx2(2:2:Nx2,2:2:Ny2,1));
sy_hx = diagonalize(1./sy2(1:2:Nx2,2:2:Ny2,2));
sy_hz = diagonalize(1./sy2(2:2:Nx2,2:2:Ny2,1));
sz_hx = diagonalize(1./sz2(1:2:Nx2,2:2:Ny2,2));
sz_hy = diagonalize(1./sz2(2:2:Nx2,1:2:Ny2,2));

% BUILD MATERIAL TENSORS
ER = [ ERxx , RX*RY'*ERxy , RX*RZ'*ERxz ; ...
       RY*RX'*ERyx , ERyy , RY*RZ'*ERyz ; ...
       RZ*RX'*ERzx , RZ*RY'*ERzy , ERzz ];
YR = [ YRxx , RX'*RY*YRxy , RX'*RZ*YRxz ; ...
       RY'*RX*YRyx , YRyy , RY'*RZ*YRyz ; ...
       RZ'*RX*YRzx , RZ'*RY*YRzy , YRzz ];
   
% CALCULATE CURL MATRICES WITH SCPML
M  = Nx*Ny;
ZZ = sparse(M,M);
CE = [ ZZ , -sz_hx*DEZ , sy_hx*DEY ; ...
       sz_hy*DEZ , ZZ , -sx_hy*DEX ; ...
       -sy_hz*DEY , sx_hz*DEX , ZZ ];
CH = [ ZZ , -sz_ex*DHZ , sy_ex*DHY ; ...
       sz_ey*DHZ , ZZ , -sx_ey*DHX ; ...
       -sy_ez*DHY , sx_ez*DHX , ZZ ];

% BUILD SCATTERED-FIELD MASKING MATRIX
nx1 = NPML(1) + 2;
nx2 = Nx - NPML(2) - 1;
ny1 = NPML(3) + 2;
ny2 = Ny - NPML(4) - 1;
Q   = ones(Nx,Ny);
Q(nx1:nx2,ny1:ny2) = 0;
Q   = diagonalize(Q);
Q   = [ Q ZZ ZZ ; ZZ Q ZZ ; ZZ ZZ Q ];

% CALCULATE WAVE MATRIX
A = CH*YR*CE - ER;

% CALCULATE SOURCE VECTOR
b = (Q*A - A*Q)*fsrc;

% Solve for Field
disp('Solving...'); drawnow;
f = A\b;

% EXTRACT FIELD COMPONENTS
fx = f(1:M);
fy = f(M+1:2*M);
fz = f(2*M+1:3*M);

% RESHAPE BACK TO TWO-DIMENSIONAL GRID
fx = reshape(fx,Nx,Ny);
fy = reshape(fy,Nx,Ny);
fz = reshape(fz,Nx,Ny);

% VISUALIZE RESULT
clf;
pcolor(xa,ya,real(fz).');
shading interp;

hold on;
phi = linspace(0,2*pi,100);
x = cos(phi);
y = sin(phi);
line(R1*x,R1*y,'Color','w');
line(R2*x,R2*y,'Color','w');
hold off;

axis equal tight;
colorbar;
set(gca,'YDir','reverse');