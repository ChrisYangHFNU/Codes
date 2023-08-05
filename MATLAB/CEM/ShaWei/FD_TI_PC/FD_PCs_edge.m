% Ref: Menglin L.N. Chen, Li Jun Jiang, Zhihao Lan, and Wei E.I. Sha
% Pseudospin-Polarized Topological Line Defects in Dielectric Photonic Crystals
% IEEE Xplore, IEEE Transactions on Antennas and Propagation, 2019.
% https://doi.org/10.1109/TAP.2019.2934816

function FD_PCs;

c=3e8;

R=6e-3;
a=2.8;

a0=a*R;
Px=a0;             % lattice constant along x
Py=sqrt(3)*a0;     % lattice constant along y

% mesh
epr1=11.7;          % dielectric ****

Nx=100;  % grids along x
Ny=round(Nx*Py/Px); % grids along y

deltax=Px/Nx;
deltay=Py/Ny;

N=Nx*Ny;        % total grid numbers;

x_c1=R; y_c1=0; % x center y center
x_c11=a0-R;

x_c2=R/2; y_c2=sqrt(3)*R/2;
x_c22=a0-R/2;

x_c3=a0/2-R/2; y_c3=sqrt(3)*(a0-R)/2;
x_c33=a0/2+R/2;

x_c4=a0/2-R; y_c4=sqrt(3)*(a0)/2;
x_c44=a0/2+R;

x_c5=a0/2-R/2; y_c5=sqrt(3)*(a0+R)/2;
x_c55=a0/2+R/2;

x_c6=R/2; y_c6=sqrt(3)*(2*a0-R)/2;
x_c66=a0-R/2;

x_c7=R; y_c7=sqrt(3)*a0;
x_c77=a0-R;

eps_arr=ones(N,1);  % permittivity of air (initial condition)

% averaged scheme for permittivity
for m=1:Nx
    for n=1:Ny
        index=(n-1)*Nx+m;
        flag1=shape_f((m-1)*deltax+0.5*deltax,(n-1)*deltay+0.5*deltay,x_c1,x_c2,x_c3,x_c4,x_c5,x_c6,x_c7,x_c11,x_c22,x_c33...
            ,x_c44,x_c55,x_c66,x_c77,y_c1,y_c2,y_c3,y_c4,y_c5,y_c6,y_c7);
        flag2=shape_f((m-1)*deltax+0.5*deltax,(n-1)*deltay-0.5*deltay,x_c1,x_c2,x_c3,x_c4,x_c5,x_c6,x_c7,x_c11,x_c22,x_c33...
            ,x_c44,x_c55,x_c66,x_c77,y_c1,y_c2,y_c3,y_c4,y_c5,y_c6,y_c7);
        flag3=shape_f((m-1)*deltax-0.5*deltax,(n-1)*deltay+0.5*deltay,x_c1,x_c2,x_c3,x_c4,x_c5,x_c6,x_c7,x_c11,x_c22,x_c33...
            ,x_c44,x_c55,x_c66,x_c77,y_c1,y_c2,y_c3,y_c4,y_c5,y_c6,y_c7);
        flag4=shape_f((m-1)*deltax-0.5*deltax,(n-1)*deltay-0.5*deltay,x_c1,x_c2,x_c3,x_c4,x_c5,x_c6,x_c7,x_c11,x_c22,x_c33...
            ,x_c44,x_c55,x_c66,x_c77,y_c1,y_c2,y_c3,y_c4,y_c5,y_c6,y_c7);
        eps_arr1(index)=eps_arr(index)+(flag1+flag2+flag3+flag4)/4*(epr1-1);
    end
end

%  construct edege state
ss1=reshape(eps_arr1,Nx,Ny);

% averaged scheme for permittivity, keep four lines: 1 2 6 7
for m=1:Nx
    for n=1:Ny
        index=(n-1)*Nx+m;
        flag1=shape_f2((m-1)*deltax+0.5*deltax,(n-1)*deltay+0.5*deltay,x_c1,x_c2,x_c6,x_c7,x_c11,x_c22,x_c66...
            ,x_c77,y_c1,y_c2,y_c6,y_c7);
        flag2=shape_f2((m-1)*deltax+0.5*deltax,(n-1)*deltay-0.5*deltay,x_c1,x_c2,x_c6,x_c7,x_c11,x_c22,x_c66...
            ,x_c77,y_c1,y_c2,y_c6,y_c7);
        flag3=shape_f2((m-1)*deltax-0.5*deltax,(n-1)*deltay+0.5*deltay,x_c1,x_c2,x_c6,x_c7,x_c11,x_c22,x_c66...
            ,x_c77,y_c1,y_c2,y_c6,y_c7);
        flag4=shape_f2((m-1)*deltax-0.5*deltax,(n-1)*deltay-0.5*deltay,x_c1,x_c2,x_c6,x_c7,x_c11,x_c22,x_c66...
            ,x_c77,y_c1,y_c2,y_c6,y_c7);
        eps_arr2(index)=eps_arr(index)+(flag1+flag2+flag3+flag4)/4*(epr1-1);
    end
end

ss2=reshape(eps_arr2,Nx,Ny);
sss=[ss1 ss1 ss1 ss2 ss1 ss1 ss1];

% introduce gap
[Nx,Ny]=size(sss);
g=36;
sss=[sss(1:Nx,1:(Ny+1)/2-1) ones(Nx,g) sss(1:Nx,(Ny+1)/2+1:Ny)]; 

% redefine Nx Ny Px Py
[Nx,Ny]=size(sss);
Px=deltax*Nx;
Py=deltay*Ny;
N=Nx*Ny;

% Get new eps_array
eps_arr=reshape(sss,N,1);

% %check the mesh
figure
pcolor(reshape(eps_arr,Nx,Ny).')
axis equal
shading interp
colormap jet

N_sam=10;      % Brillouin sampling points ***
N_max=50;     % maximum band number ***

% sampling at the edge of Brillouin zone
kx=linspace(0,pi/Px,N_sam); ky=linspace(0,0,N_sam);

omega_f=zeros(N_max,length(kx));

%calculate eigenvalue
for m=1:length(kx)
    length(kx)-m
    [omega_f(:,m)]=eigs_PCs(Nx,Ny,Px,Py,eps_arr,kx(m),ky(m),N_max,0,0);
end

%draw band diagram
figure;
for m=1:N_max
    plot(omega_f(m,:)*c/(2*pi),'rx-')
    hold on
    plot(1:length(kx),7.941e9*ones(1,length(kx)),'k--')
    plot(1:length(kx),8.67e9*ones(1,length(kx)),'k--')
    axis([1 length(kx) 7.5e9 9e9])
    xlabel('Bloch Wavenumber')
    ylabel('Eigen-frequency')
end

% draw eigenmode at some k
NN=2;
omega_r=eigs_PCs(Nx,Ny,Px,Py,eps_arr,kx(NN),ky(NN),N_max,NN,1);


% solve the eigen-equation with a fixed kx and ky
% N_max is the maximum eigenvalues of interest
% Px Py are lattice constants along x and y directions
% flag=1 for drawing eigenmodes
function [omega]=eigs_PCs(Nx,Ny,Px,Py,eps_arr,kx,ky,N_max,NN,flag)

deltax=Px/Nx;
deltay=Py/Ny;

N=Nx*Ny;      % total grid numbers;

% sparse matrix structure (line number, column number, and values)
% five-point difference
line=zeros(5,N);
column=zeros(5,N);
value=zeros(5,N);

% center region
for m=2:Nx-1
    for n=2:Ny-1
        index=(n-1)*Nx+m;
        line(:,index)=[index,index,index,index,index];
        column(:,index)=[index,index+1,index-1,index+Nx,index-Nx];
        value(:,index)=[2*(1/deltax^2+1/deltay^2)/eps_arr(index),...
            -1/deltax^2/eps_arr(index),-1/deltax^2/eps_arr(index),...
            -1/deltay^2/eps_arr(index),-1/deltay^2/eps_arr(index)];
        
    end
end

% left
for m=1:1
    for n=2:Ny-1
        index=(n-1)*Nx+m;
        line(:,index)=[index,index,index,index,index];
        column(:,index)=[index,index+1,(n-1)*Nx+Nx,index+Nx,index-Nx];
        value(:,index)=[2*(1/deltax^2+1/deltay^2)/eps_arr(index),...
            -1/deltax^2/eps_arr(index),-exp(j*kx*Px)/deltax^2/eps_arr(index),...
            -1/deltay^2/eps_arr(index),-1/deltay^2/eps_arr(index)];
        
    end
end

% right
for m=Nx:Nx
    for n=2:Ny-1
        index=(n-1)*Nx+m;
        line(:,index)=[index,index,index,index,index];
        column(:,index)=[index,(n-1)*Nx+1,index-1,index+Nx,index-Nx];
        value(:,index)=[2*(1/deltax^2+1/deltay^2)/eps_arr(index),...
            -exp(-j*kx*Px)/deltax^2/eps_arr(index),-1/deltax^2/eps_arr(index),...
            -1/deltay^2/eps_arr(index),-1/deltay^2/eps_arr(index)];
        
    end
end

% top
for m=2:Nx-1
    for n=1:1
        index=(n-1)*Nx+m;
        line(:,index)=[index,index,index,index,index];
        column(:,index)=[index,index+1,index-1,index+Nx,(Ny-1)*Nx+m];
        value(:,index)=[2*(1/deltax^2+1/deltay^2)/eps_arr(index),...
            -1/deltax^2/eps_arr(index),-1/deltax^2/eps_arr(index),...
            -1/deltay^2/eps_arr(index),-exp(j*ky*Py)/deltay^2/eps_arr(index)];
        
    end
end

% bottom
for m=2:Nx-1
    for n=Ny:Ny
        index=(n-1)*Nx+m;
        line(:,index)=[index,index,index,index,index];
        column(:,index)=[index,index+1,index-1,(1-1)*Nx+m,index-Nx];
        value(:,index)=[2*(1/deltax^2+1/deltay^2)/eps_arr(index),...
            -1/deltax^2/eps_arr(index),-1/deltax^2/eps_arr(index),...
            -exp(-j*ky*Py)/deltay^2/eps_arr(index),-1/deltay^2/eps_arr(index)];
        
    end
end

% left-top
m=1; n=1;
index=(n-1)*Nx+m;
line(:,index)=[index,index,index,index,index];
column(:,index)=[index,index+1,(n-1)*Nx+Nx,index+Nx,(Ny-1)*Nx+m];
value(:,index)=[2*(1/deltax^2+1/deltay^2)/eps_arr(index),...
    -1/deltax^2/eps_arr(index),-exp(j*kx*Px)/deltax^2/eps_arr(index),...
    -1/deltay^2/eps_arr(index),-exp(j*ky*Py)/deltay^2/eps_arr(index)];

% left-bottom
m=1; n=Ny;
index=(n-1)*Nx+m;
line(:,index)=[index,index,index,index,index];
column(:,index)=[index,index+1,(n-1)*Nx+Nx,(1-1)*Nx+m,index-Nx];
value(:,index)=[2*(1/deltax^2+1/deltay^2)/eps_arr(index),...
    -1/deltax^2/eps_arr(index),-exp(j*kx*Px)/deltax^2/eps_arr(index),...
    -exp(-j*ky*Py)/deltay^2/eps_arr(index),-1/deltay^2/eps_arr(index)];

% right-top
m=Nx; n=1;
index=(n-1)*Nx+m;
line(:,index)=[index,index,index,index,index];
column(:,index)=[index,(n-1)*Nx+1,index-1,index+Nx,(Ny-1)*Nx+m];
value(:,index)=[2*(1/deltax^2+1/deltay^2)/eps_arr(index),...
    -exp(-j*kx*Px)/deltax^2/eps_arr(index),-1/deltax^2/eps_arr(index),...
    -1/deltay^2/eps_arr(index),-exp(j*ky*Py)/deltay^2/eps_arr(index)];

% right-bottom
m=Nx; n=Ny;
index=(n-1)*Nx+m;
line(:,index)=[index,index,index,index,index];
column(:,index)=[index,(n-1)*Nx+1,index-1,(1-1)*Nx+m,index-Nx];
value(:,index)=[2*(1/deltax^2+1/deltay^2)/eps_arr(index),...
    -exp(-j*kx*Px)/deltax^2/eps_arr(index),-1/deltax^2/eps_arr(index),...
    -exp(-j*ky*Py)/deltay^2/eps_arr(index),-1/deltay^2/eps_arr(index)];

% Solve eigen-equation
A=sparse(line,column,value,N,N);  % establish sparse matrix
[V,D] = eigs(A,N_max,'sm');       % find eigenmodes and eigen-frequency
[omega,order]=sort(sqrt(abs(real(diag(D)))));  %  from small to large eigenvalues

% draw eigenmodes
field_2d=zeros(Ny,Nx);
if (flag==1) % flag control
    
    for mm=40:42  %  choose which band to plot ****
        field_1d=V(:,order(mm));  % eigenmodes from low to high order
        for m=1:Ny
            for n=1:Nx
                index=(m-1)*Nx+n;
                field_2d(m,n)=field_1d(index);
            end
        end
        figure(mm)
        pcolor(real(field_2d));
        shading interp
        colormap jet
        colorbar
        axis([1 Nx 1 Ny])
        axis equal
    end
end

%  model the unit cell
function flag=shape_f(x,y,xc1,xc2,xc3,xc4,xc5,xc6,xc7,xc11,xc22,xc33,xc44,xc55,xc66,xc77,yc1,yc2,yc3,yc4,yc5,yc6,yc7)
flag=0;
radius=2e-3;  % radius of circle in unit cells ***
if (((x-xc1)^2+(y-yc1)^2)<radius^2 | ((x-xc2)^2+(y-yc2)^2)<radius^2 | ((x-xc3)^2+(y-yc3)^2)<radius^2 | ...
        ((x-xc4)^2+(y-yc4)^2)<radius^2 | ((x-xc5)^2+(y-yc5)^2)<radius^2 | ((x-xc6)^2+(y-yc6)^2)<radius^2 | ...
        ((x-xc7)^2+(y-yc7)^2)<radius^2 | ((x-xc11)^2+(y-yc1)^2)<radius^2 | ((x-xc22)^2+(y-yc2)^2)<radius^2 |...
        ((x-xc33)^2+(y-yc3)^2)<radius^2 | ((x-xc44)^2+(y-yc4)^2)<radius^2 | ((x-xc55)^2+(y-yc5)^2)<radius^2 |...
        ((x-xc66)^2+(y-yc6)^2)<radius^2 | ((x-xc77)^2+(y-yc7)^2)<radius^2)
    flag=1;
end

%  model the unit cell
function flag2=shape_f2(x,y,xc1,xc2,xc3,xc4,xc11,xc22,xc33,xc44,yc1,yc2,yc3,yc4)
flag2=0;
radius=2e-3;  % radius of circle in unit cells ***
if (((x-xc1)^2+(y-yc1)^2)<radius^2 | ((x-xc2)^2+(y-yc2)^2)<radius^2 | ((x-xc3)^2+(y-yc3)^2)<radius^2 | ...
        ((x-xc4)^2+(y-yc4)^2)<radius^2  | ((x-xc11)^2+(y-yc1)^2)<radius^2 | ((x-xc22)^2+(y-yc2)^2)<radius^2 |...
        ((x-xc33)^2+(y-yc3)^2)<radius^2 | ((x-xc44)^2+(y-yc4)^2)<radius^2)
    flag2=1;
end
