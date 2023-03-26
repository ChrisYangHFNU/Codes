% Chapter7_photonicbanddiagram.m

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
Nx     = 40;
Ny     = Nx;
NBETA  = 100;
NBANDS = 5;
wnmax  = 0.6;

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

% RECIPROCAL LATTICE VECTORS
T1 = (2*pi/a) * [ 1 ; 0 ];
T2 = (2*pi/a) * [ 0 ; 1 ];

% KEY POINTS OF SYMMETRY
G = [ 0 ; 0 ];
X = 0.5*T1;
M = 0.5*T1 + 0.5*T2;

% CHOOSE PATH AROUND IBZ
KP = [ G X M G ];
KL = { '\Gamma' 'X' 'M' '\Gamma' };

% DETERMINE LENGTH OF IBZ PERIMETER
NKP  = length(KP);
LIBZ = 0;
for m = 1 : NKP-1
    LIBZ = LIBZ + norm(KP(:,m+1) - KP(:,m));
end

% GENERATE LIST OF POINTS AROUND IBZ
dibz  = LIBZ/NBETA;
BETA  = KP(:,1);
KT    = 1;
NBETA = 1;
for m = 1 : NKP-1
    dK      = KP(:,m+1) - KP(:,m);
    N       = ceil(norm(dK)/dibz);
    BETA    = [ BETA , KP(:,m) + dK*[1:N]/N ];
    NBETA   = NBETA + N;
    KT(m+1) = NBETA;
end

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
WNTE = zeros(NBANDS,NBETA);
WNTM = zeros(NBANDS,NBETA);

%
% MAIN LOOP -- ITERATE OVER IBZ
%
for nbeta = 1 : NBETA
    
    % Get Next Block Wave Vector
    beta = BETA(:,nbeta);
    
    % Build Derivative Matrices
    NS  = [Nx Ny];
    RES = [dx dy];
    BC  = [1 1];
    [DEX,DEY,DHX,DHY] = yeeder2d(NS,RES,BC,beta);
        
    % TM Mode Analysis
    A = - DHX/URyy*DEX - DHY/URxx*DEY;
    B = ERzz;
    D = eigs(A,B,NBANDS,0);
    D = sort(D);
    WNTM(:,nbeta) = D(1:NBANDS);

    % TE Mode Analysis
    A = - DEX/ERyy*DHX - DEY/ERxx*DHY;
    B = URzz;
    D = eigs(A,B,NBANDS,0);
    D = sort(D);
    WNTE(:,nbeta) = D(1:NBANDS);
    
end

% NORMALIZE THE FREQUENCIES
WNTE = a/(2*pi)*real(sqrt(WNTE));
WNTM = a/(2*pi)*real(sqrt(WNTM));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DRAW PROFESSIONAL BAND DIAGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PREPARE FIGURE WINDOW
subplot(1,3,2:3);

% DRAW BAND DATA
plot(1:NBETA,WNTM.','.b');
hold on;
plot(1:NBETA,WNTE.','.r');
hold off;

% SET VIEW
xlim([1 NBETA]);
ylim([0 wnmax]);
set(gca,'XTick',KT,'XTickLabel',KL);
xlabel('Bloch Wave Vector $\vec{\beta}$',...
       'Interpreter','LaTex');
ylabel('Frequency $\omega_{\textrm{n}} = a/\lambda_0$',...
       'Interpreter','LaTex');
title('Photonic Band Diagram');


