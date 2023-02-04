% Chapter7_IFCs.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PHOTONIC CRYSTAL PARAMETERS
a      = 1;
r      = 0.48*a;
erhole = 1.0;
erfill = 12.0;

% FDFD PARAMETERS
Nx     = 60;
Ny     = Nx;
NBx    = 41;
NBy    = NBx;
NBANDS = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCULATE OPTIMIZED GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% YEE GRID
dx = a/Nx;
dy = a/Ny;

% 2X GRID
Nx2 = 2*Nx;             dx2 = dx/2;
Ny2 = 2*Ny;             dy2 = dy/2;

% CALCULATE 2X MESHGRID
xa2 = [1:Nx2]*dx2;      xa2 = xa2 - mean(xa2);
ya2 = [1:Ny2]*dy2;      ya2 = ya2 - mean(ya2);
[Y2,X2] = meshgrid(ya2,xa2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD UNIT CELL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BUILD UNIT CELL
ER2 = (X2.^2 + Y2.^2) <= r^2;
ER2 = erfill + (erhole - erfill)*ER2;
UR2 = ones(Nx2,Ny2);

% EXTRACT YEE GRID MATERIAL ARRAYS
ERxx = ER2(2:2:Nx2,1:2:Ny2);
ERyy = ER2(1:2:Nx2,2:2:Ny2);
ERzz = ER2(1:2:Nx2,1:2:Ny2);
URxx = UR2(1:2:Nx2,2:2:Ny2);
URyy = UR2(2:2:Nx2,1:2:Ny2);
URzz = UR2(2:2:Nx2,2:2:Ny2);

% SHOW UNIT CELL
subplot(131);
imagesc(xa2,ya2,ER2.');
axis equal tight;
colorbar;
title('Unit Cell');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALCUALTE LIST OF BLOCH WAVE VECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CREATE BETA MESHGRID
bx      = linspace(-pi/a,+pi/a,NBx);
by      = linspace(-pi/a,+pi/a,NBy);
[BY,BX] = meshgrid(by,bx);

% EXTRACT IBZ AREA
ind_ibz  = find(BX>=0 & BY>=0 & BX>=BY);
indx_ibz = 1 + mod(ind_ibz-1,NBx);
indy_ibz = 1 + floor((ind_ibz-1)/NBx);
NBETA    = length(ind_ibz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM FDFD ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DIAGONALIZE MATERIAL TENSORS
ERxx = diag(sparse(ERxx(:)));
ERyy = diag(sparse(ERyy(:)));
ERzz = diag(sparse(ERzz(:)));
URxx = diag(sparse(URxx(:)));
URyy = diag(sparse(URyy(:)));
URzz = diag(sparse(URzz(:)));

% INITIALIZE BAND DATA
WNTE = zeros(NBx,NBy,NBANDS);
WNTM = zeros(NBx,NBy,NBANDS);

%
% MAIN LOOP -- ITERATE OVER IBZ
%
for nbeta = 1 : NBETA
    
    % Get Next Block Wave Vector
    beta = [ BX(ind_ibz(nbeta)) ; BY(ind_ibz(nbeta)) ];
        
    % Build Derivative Matrices
    NS  = [Nx Ny];
    RES = [dx dy];
    BC  = [1 1];
    [DEX,DEY,DHX,DHY] = yeeder2d(NS,RES,BC,beta);
        
    % TM Mode Analysis
    A = - DHX/URyy*DEX - DHY/URxx*DEY;
    B = ERzz;
    D = eigs(A,B,NBANDS,0);
    D = real(sqrt(D));
    D = sort(D);
    WNTM(indx_ibz(nbeta),indy_ibz(nbeta),:) ...
                         = a/(2*pi)*D(1:NBANDS);

    % TE Mode Analysis
    A = - DEX/ERyy*DHX - DEY/ERxx*DHY;
    B = URzz;
    D = eigs(A,B,NBANDS,0);
    D = real(sqrt(D));
    D = sort(D);
    WNTE(indx_ibz(nbeta),indy_ibz(nbeta),:) ...
                         = a/(2*pi)*D(1:NBANDS);
        
end

% UNFOLD DATA
for nb = 1 : NBANDS
    wn = WNTE(:,:,nb);
    wn = wn + wn.';
    wn = wn - 0.5*diag(diag(wn));
    wn(1:floor(NBx/2),:) = wn(NBx:-1:2+floor(NBx/2),:);
    wn(:,1:floor(NBy/2)) = wn(:,NBy:-1:2+floor(NBy/2));
    WNTE(:,:,nb) = wn;
    
    wn = WNTM(:,:,nb);
    wn = wn + wn.';
    wn = wn - 0.5*diag(diag(wn));
    wn(1:floor(NBx/2),:) = wn(NBx:-1:2+floor(NBx/2),:);
    wn(:,1:floor(NBy/2)) = wn(:,NBy:-1:2+floor(NBy/2));
    WNTM(:,:,nb) = wn;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DRAW PRETTY ISO FREQUENCY CONTOURS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SHOW TE BAND 1
subplot(2,3,2);
wn = WNTE(:,:,1);
v = linspace(min(wn(:)),max(wn(:)),20);
contour(bx,by,wn.',v,'LineWidth',2);
axis equal tight;
colorbar;
xlabel('\beta_x');
ylabel('\beta_y');
title('TE Mode 1');
    
% SHOW TE BAND 2
subplot(2,3,3);
wn = WNTE(:,:,2);
v = linspace(min(wn(:)),max(wn(:)),20);
contour(bx,by,wn.',v,'LineWidth',2);
axis equal tight;
colorbar;
xlabel('\beta_x');
ylabel('\beta_y');
title('TE Mode 2');

% SHOW TM BAND 1
subplot(2,3,5);
wn = WNTM(:,:,1);
v = linspace(min(wn(:)),max(wn(:)),20);
contour(bx,by,wn.',v,'LineWidth',2);
axis equal tight;
colorbar;
xlabel('\beta_x');
ylabel('\beta_y');
title('TM Mode 1');
    
% SHOW TM BAND 2
subplot(2,3,6);
wn = WNTM(:,:,2);
v = linspace(min(wn(:)),max(wn(:)),20);
contour(bx,by,wn.',v,'LineWidth',2);
axis equal tight;
colorbar;
xlabel('\beta_x');
ylabel('\beta_y');
title('TM Mode 2');