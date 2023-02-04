% Chapter3_mvfdm.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% GRID
a  = 0;
b  = 10;
Nx = 200;
dx = (b - a)/Nx;
xa = [0.5:Nx-0.5]*dx;

% BUILD DERIVATIVE MATRICES
d   = ones(Nx,1)/dx;

DFX = sparse(Nx,Nx);
DFX = spdiags(-d,0,DFX);
DFX = spdiags(+d,+1,DFX);

DGX = -DFX';

% BUILD MATRIX EQUATION
I = speye(Nx,Nx);
A = DGX*DFX + 6*I;
b = zeros(Nx,1);

% ADD BOUNDARY VALUES
A(1,:) = 0;
A(1,1) = 1;
b(1)   = 10;

A(Nx,:)  = 0;
A(Nx,Nx) = 1;
b(Nx)    = 1;

% SOLVE
f = A\b;
g = -(1/3)*DFX*f;

% PLOT RESULTS
plot(xa,f,'-k','LineWidth',2);
hold on;
plot(xa,g,'--k','LineWidth',2);
hold off;
xlabel('$x$','Interpreter','LaTex');
legend('$f(x)$','$g(x)$',...
       'Interpreter','LaTex',...
       'Location','NorthEastOutside');
set(gca,'FontSize',18,'LineWidth',2);