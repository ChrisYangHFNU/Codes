% Chapter3_fdm1d.m

% INITIALIZE MATLAB
close all;
clc;
clear all;

% DEFINE AND CALCULATE GRID
a  = 0;
b  = 10;
Nx = 100;
xa = linspace(a,b,Nx);
dx = xa(2) - xa(1);

% BUILD DERIVATIVE MATRIX DX
d  = ones(Nx,1)/(2*dx);
DX = spdiags(-d,-1,sparse(Nx,Nx));
DX = spdiags(+d,+1,DX);
DX(Nx,Nx-1:Nx) = [-1 +1]/dx;

% BUILD MATRIX EQUATION A*f=b
A = DX;
b = -(1/3)*ones(Nx,1);

% ADD BOUNDARY VALUE f(0)=1
A(1,:) = 0;
A(1,1) = 1;
b(1)   = 1;

% SOLVE PROBLEM
f = A\b;

% PLOT SOLUTION
plot(xa,f);
axis equal tight;