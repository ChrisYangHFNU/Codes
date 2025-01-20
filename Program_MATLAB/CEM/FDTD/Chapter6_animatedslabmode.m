% Chapter6_animatedslabmode.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% UNITS 
micrometers = 1;
nanometers  = 1e-3 * micrometers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WAVELENGTH AND MODE
lam0 = 1.55 * micrometers;
MODE = 'H';

% SLAB WAVEGUIDE
a  = 1500 * nanometers;
n1 = 1.0;
n2 = 2.0;
n3 = 1.5;

% GRID PARAMETERS
nmax = max([n1 n2 n3]);
NRES = 20;
b    = 3*lam0;

% NUMBER OF MODES TO CALCULATE
NMODES = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIRST GUESS AT GRID RESOLUION
dx = lam0/nmax/NRES;

% SNAP GRID TO CRITICAL DIMENSIONS
nx = ceil(a/dx);
dx = a/nx;

% CALCULATE TOTAL GRID SIZE
Sx = b + a + b;
Nx = ceil(Sx/dx);
Sx = Nx*dx;

% GRID AXIS
xa = [1:Nx]*dx;
xa = xa - mean(xa);

% 2X GRID
Nx2 = 2*Nx;
dx2 = dx/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD SLAB WAVEGUIDE ONTO GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INITIALIZE TO FREE SPACE
ER2 = ones(Nx2,1);
UR2 = ones(Nx2,1);

% CALCULATE START AND STOP ARRAY INDICES
nx1 =   1 + ceil(b/dx2);
nx2 = nx1 + round(a/dx2) - 1;

% BUILD SLAB WAVEGUIDE
ER2(1:nx1-1)   = n1^2;
ER2(nx1:nx2)   = n2^2;
ER2(nx2+1:Nx2) = n3^2;

% EXTRACT YEE GRID ARRAYS
ERxx = ER2(2:2:Nx2);
ERyy = ER2(1:2:Nx2);
ERzz = ER2(1:2:Nx2);
URxx = UR2(1:2:Nx2);
URyy = UR2(2:2:Nx2);
URzz = UR2(2:2:Nx2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FINITE-DIFFERENCE ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FORM DIAGONAL MATERIALS ARRAYS
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));
URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

% BUILD DERIVATIVE MATRICES
k0  = 2*pi/lam0;
NS  = [Nx 1];
RES = [dx 1];
BC  = [0 0];
[DEX,DEY,DHX,DHY] = yeeder2d(NS,k0*RES,BC);

% BUILD EIGEN-VALUE PROBLEM
if MODE=='E'
    A = -(DHX/URzz*DEX + ERyy);
    B = inv(URxx);
else
    A = -(DEX/ERzz*DHX + URyy);
    B = inv(ERxx);
end

% SOLVE EIGEN-VALUE PROBLEM
ev      = -n2^2;
[Fy,D2] = eigs(A,B,NMODES,ev);
D       = sqrt(D2);
NEFF    = -1i*diag(D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ANIMATE A SLAB WAVEGUIDE MODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE GIF ANIMATION
NFRAMES  = 40;
gif_name = 'slab_animation.gif';

% PICK A MODE
m = 2;

% CALCULATE A 2D GRID
Sy = 2*Sx;
dy = dx;
Ny = ceil(Sy/dy);
Sy = Ny*dy;
ya = [0:Ny-1]*dy;
[Y,X] = meshgrid(ya,xa);

% PROPAGATE MODE ACROSS GRID
F = zeros(Nx,Ny);
for ny = 1 : Ny
    SM(:,ny) = Fy(:,m)*exp(-1i*k0*NEFF(m)*ny*dy);
end

% VISUALIZE USING IMAGESC
for nframe = 1 : NFRAMES
    
    % Add Phase
    phase = 2*pi*nframe/NFRAMES;
    fi = SM*exp(1i*phase);
    
    % Prepare Graphics
    clf;
    hold on;
    
    % Draw Field
    pcolor(ya,xa,real(fi));
    shading interp;
    
    % Calculate Transparencies
    fmin = 0.05;
    fmax = 0.5;
    nmin = min([n1 n2 n3]);
    nmax = max([n1 n2 n3]);
    fa1  = fmin + (fmax - fmin)*(n1 - nmin)/(nmax - nmin);
    fa2  = fmin + (fmax - fmin)*(n2 - nmin)/(nmax - nmin);
    fa3  = fmin + (fmax - fmin)*(n3 - nmin)/(nmax - nmin);
    
    % Draw Slab Waveguide
    x = [ ya(1) ya(Ny) ya(Ny) ya(1) ya(1) ];
    y = [ a/2 , a/2 , b+a/2 , b+a/2 , a/2 ];
    fill(x,-y,'w','FaceAlpha',fa1);
    fill(x,+y,'w','FaceAlpha',fa2);
    y = [ -a/2 , -a/2 , a/2 , a/2 , -a/2 ];
    fill(x,-y,'w','FaceAlpha',fa3);
    
    % Set Graphics View
    hold off;
    axis equal tight off;
    set(gca,'YDir','reverse');
    drawnow;
    
    % Capture Frame
    F = getframe(gca);
    F = frame2im(F);
    [ind,cmap] = rgb2ind(F,256,'nodither');
    
    % Add Frame to GIF
    if nframe == 1
        imwrite(ind,cmap,gif_name,'gif',...
                'DelayTime',0,'Loopcount',inf);
    else
        imwrite(ind,cmap,gif_name,'gif',...
                'DelayTime',0,'WriteMode','append');
    end
end