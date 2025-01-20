%  A simple finite-difference method to get potential distribution in a
%  rectangular box. Boundary conditions are V0=0 except for V0=10 at y=b.
%  The code was written by Wei Sha from HKU
%  Email: wsha@eee.hku.hk 

clc;clear

%  boundary value at y=b
V0=10;

% line and column number
Nx=80;
Ny=60;
N=Nx*Ny;

% setting for sparse matrix
a=ones(N,5);  % line 
b=ones(N,5);  % column
c=zeros(N,5); % value

% voltage vector
v=zeros(N,1); 

% matrix filling (five-point difference)
for n=1:Ny
    for m=1:Nx
        index=(n-1)*Nx+m;
        
        %  middle point
        a(index,1)=index;
        b(index,1)=index;
        c(index,1)=4;
        
        if (m-1>=1)
            a(index,2)=index;
            b(index,2)=index-1;  %  left
            c(index,2)=-1;
        end
        
        if (m+1<=Nx)
            a(index,3)=index;
            b(index,3)=index+1;  %  right
            c(index,3)=-1;
        end
        
        if (n+1<=Ny)
            a(index,4)=index;
            b(index,4)=index+Nx; %  down
            c(index,4)=-1;
        end
        
        if (n-1>=1)
            a(index,5)=index;
            b(index,5)=index-Nx; %  up
            c(index,5)=-1;
        end
        
    end
end

% voltage filling V0=10 at y=b
for m=1:Nx
    index=(Ny-1)*Nx+m;
    v(index)=V0;
end

% generate sparse matrix
X=sparse(a,b,c,N,N);

% get the solution
res=X\v;

% convert 1d solution to 2d solution
yy=zeros(Ny,Nx);

for n=1:Ny
    for m=1:Nx
        index=(n-1)*Nx+m;
        yy(n,m)=res(index);
    end
end

% get analytical solution by separation of variables
yy_ref=zeros(Ny,Nx);
for n=1:Ny
    for m=1:Nx
        for q=1:2:101
            yy_ref(n,m)=yy_ref(n,m)+4*10/pi*sin(q*pi*m/(Nx+1))*sinh(q*pi*n/(Nx+1))/(q*sinh(q*pi*(Ny+1)/(Nx+1)));
        end
    end
end

%  show results
figure(1)
contour(yy);
xlabel('x')
ylabel('y')
title('equipotential lines (numerical)')

figure(2)
contour(yy_ref);
xlabel('x')
ylabel('y')
title('equipotential lines (analytical)')

%  error between analytical solution and numerical solution
norm(yy_ref-yy,'fro')/norm(yy,'fro')



