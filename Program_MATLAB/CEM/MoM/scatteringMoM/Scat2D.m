function[In,sigma,phi,Data]=Scat2D(FileName,phi_i,Type)
% Solves 2D scattering problem using MoM with Pieacewise Linear basis and 
% testing functions (Galerkin's method) and calculates RCS.
% Call [In,sigma,phi,Data]=Scat2D(FileName,phi_i,Type)
%
% Input:-
% FileName  : a 1D mesh of PEC scatter file with format .dat (see sample files) 
% phi_i     : incident field angle in degrees
% Type      : type of incident field "TE" or "TM"
%
% Output:-
% In        : Induced current on segments
% sigma     : RCS data
% phi       : scattered field angle phi data
% Data      : segments data with format: xn yn xn- yn- xn+ yn+
%% Main program calls functions
fprintf('Processing the mesh...\n');
tic; Data=MakeMesh(FileName); toc; 
if Type=="TE"
    fprintf('\nTE solution:\n');
    [In,phi,sigma]=TE(phi_i,Data);
end
if Type=="TM"
    fprintf('\nTM solution:\n');
    [In,phi,sigma]=TM(phi_i,Data);
end
fprintf('\nComplete!\n');
end
%% Process the mesh
function[Mesh]=MakeMesh(str)
%% Call read .dat file
DATA        =   ReadFile(str);
[Nseg,Nodes,Segments]=NodesR(DATA);
%% Create PWL segments
Segments    =   [Segments;[Segments(1,1) Segments(1,2)]];
Mesh        =   [];
for i=1:Nseg
    N1n        	=   Segments(i,1);
    N2n        	=   Segments(i,2);
    N1p        	=   Segments(i+1,1);
    N2p        	=   Segments(i+1,2);
    if N1n==N2p
        xn      =   Nodes(Segments(i,1),1);
        yn      =   Nodes(Segments(i,1),2);
        xm      =   Nodes(Segments(i+1,1),1);
        ym      =   Nodes(Segments(i+1,1),2);
        xp      =   Nodes(Segments(i,2),1);
        yp      =   Nodes(Segments(i,2),2);
        Mesh    =   [Mesh;[xn yn xm ym xp yp]];
    end
    if N2n==N1p
        xn      =   Nodes(Segments(i,2),1);
        yn      =   Nodes(Segments(i,2),2);
        xm      =   Nodes(Segments(i,1),1);
        ym      =   Nodes(Segments(i,1),2);
        xp      =   Nodes(Segments(i+1,2),1);
        yp      =   Nodes(Segments(i+1,2),2);
        Mesh    =   [Mesh;[xn yn xm ym xp yp]];
    end
end
%%
end
%% Extract segments from nodes/connectivity list
function[Nseg,Nodes,Segments]=NodesR(DATA)
Nnodes      =   DATA(1,1);
DATA        =   DATA(2:end,:);
Nodes       =   zeros(Nnodes,3);
for i=1:Nnodes
    Nodes(i,:)	=   [DATA(i,2) DATA(i,3) DATA(i,4)];
end
DATA_Con    =   DATA(Nnodes+1:end,:);
[Nseg,~]    =   size(DATA_Con);
Segments    =   zeros(Nseg,2);
for i=1:Nseg
    Segments(i,:)	=   [DATA_Con(i,3) DATA_Con(i,4)];
end   
%% Center the object at the origin
% x_min       =   min(Nodes(:,1));
% x_max       =   max(Nodes(:,1));
% y_min       =   min(Nodes(:,2));
% y_max       =   max(Nodes(:,2));
% Nodes(:,1)	=   Nodes(:,1)-0.5*(x_max+x_min);
% Nodes(:,2)	=   Nodes(:,2)-0.5*(y_max+y_min);
%% Check the direction of contour C
sum         =   0;
for i=1:Nseg
    x1      =   Nodes(Segments(i,1),1);
    x2      =   Nodes(Segments(i,2),1);
    y1      =   Nodes(Segments(i,1),2);
  	y2      =   Nodes(Segments(i,2),2);
    L       =   sqrt((x1-x2)^2+(y1-y2)^2);
    phi     =   atan((y2-y1)/(x2-x1));
    phi(isnan(phi))=0;
    if (x2-x1)<0
        phi     =   phi+pi;
    end
    Mag     =   sqrt((x1+0.5*L*cos(phi))^2+(y1+0.5*L*sin(phi))^2);
    rho_x   =   (x1+0.5*L*cos(phi))/Mag;
    rho_y   =   (y1+0.5*L*sin(phi))/Mag;
    lx      =   cos(phi);
    ly      =   sin(phi);
    sum     =   sum+rho_x*ly-rho_y*lx;
end
%% Correct the direction of contour C 
if sum<0
    Segments    =   [Segments(:,2) Segments(:,1)];
end
Plot()
%% Plot the nodes and segments
function[]=Plot()
    figure()
    plot(Nodes(:,1),Nodes(:,2),'.k','MarkerSize',15)
    hold on
    for ii=1:Nseg
        x1      =   Nodes(Segments(ii,1),1);
        x2      =   Nodes(Segments(ii,2),1);
        y1      =   Nodes(Segments(ii,1),2);
        y2      =   Nodes(Segments(ii,2),2);
        plot([x1 x2],[y1 y2],'-k','LineWidth',1)
        L       =   sqrt((x1-x2)^2+(y1-y2)^2);
        phi     =   atan((y2-y1)/(x2-x1));
        phi(isnan(phi))=0;
        if (x2-x1)<0
            phi     =   phi+pi;
        end
        axis equal
        xlabel('$x/\lambda$','Interpret','Latex','FontSize',14)
        ylabel('$y/\lambda$','Interpret','Latex','FontSize',14)
        title('Nodes and segments','Interpret','Latex','FontSize',14)
        set(gca,'TickLabel','Latex','FontSize',14)
    end
    hold off
end   
end
%% Read raw .dat file and prepare DATA matrix
function[DATA]=ReadFile(str)
fID         =	fopen(str,'r');
RDATA   	=	fscanf(fID,'%f');
[Nf,~]    	=   size(RDATA);
Nd          =   0.25*(Nf+2);
DATA        =   zeros(Nd,4);
count       =   1;
for i=1:Nd
    if i==1
        D1          =   RDATA(count,1);count=count+1;
        D2          =   RDATA(count,1);count=count+1;
        DATA(i,:)   =   [ D1 D2 0 0];
    else
        D1          =   RDATA(count,1);count=count+1;
        D2          =   RDATA(count,1);count=count+1;
        D3          =   RDATA(count,1);count=count+1;
        D4          =   RDATA(count,1);count=count+1;
        DATA(i,:)   =   [D1 D2 D3 D4];
    end
end
end
%% TE solution
function[In,phi,sigma]=TE(phi_i,Data)
%%
alpha   =   0.5;
%%
[Nss,~] =   size(Data);
%% Call EFIE and MFIE
eta     =   120*pi;
tic;
fprintf('Solving the EFIE...\n');
[Zmn_EFIE,Vm_EFIE]=MM_EFIE_TE(Data,Nss,phi_i);
toc;
tic;
fprintf('Solving the MFIE...\n');
[Zmn_MFIE,Vm_MFIE]=MM_MFIE_TE(Data,Nss,phi_i);
toc; 
%% Solve for currents
fprintf('Solving the CFIE...\n');
Zmn     =   alpha*Zmn_EFIE+(1-alpha)*eta*Zmn_MFIE;
Vm      =   alpha*Vm_EFIE+(1-alpha)*eta*Vm_MFIE;
In      =   Zmn\Vm;
%% Find RCS
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
%%
[Nss,~] =   size(Data);
phi     =   linspace(-pi,pi,1e3);
sigma   =   0;
for n=1:Nss
    xn      =   Data(n,1);
    yn      =   Data(n,2);
    xnm    	=   Data(n,3);
    ynm    	=   Data(n,4);
    xnp    	=   Data(n,5);
    ynp    	=   Data(n,6);
    Lxnm   	=   xn-xnm;
    Lynm   	=   yn-ynm;
    Lxnp   	=   xnp-xn;
    Lynp   	=   ynp-yn;
    kap_np  =   exp(j*k*(xnp*cos(phi)+ynp*sin(phi)));
    kap_nm  =   exp(j*k*(xnm*cos(phi)+ynm*sin(phi)));
    dot_p   =   -Lxnp*sin(phi)+Lynp*cos(phi);
    dot_m   =   -Lxnm*sin(phi)+Lynm*cos(phi);
    ap      =   Lxnp*cos(phi)+Lynp*sin(phi);
    am      =   Lxnm*cos(phi)+Lynm*sin(phi);
    Tp      =   (exp(-j*k*ap).*(1+j*k*ap)-1)./(k^2*ap.^2);
    Tm      =   (exp(+j*k*am).*(1-j*k*am)-1)./(k^2*am.^2); 
    Tp(abs(ap)<1e-5)=0.5;
    Tm(abs(am)<1e-5)=0.5;
    %% Standard
    Term    =   dot_p.*kap_np.*Tp+dot_m.*kap_nm.*Tm;
    sigma   =   sigma+(eta/2)*In(n,1)*Term;
end
sigma   =   2*pi*abs(sigma).^2;
%% Plot RCS
figure()
plot(phi*180/pi,10*log10(sigma),'-k','LineWidth',1)
xlabel('$\varphi$ [deg]','Interpret','Latex','FontSize',14)
ylabel('$\sigma/\lambda$ [dB]','Interpret','Latex','FontSize',14)
title('TE','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
xlim([-180 180])
%%
end
%% Solve EFIE for TE
function[Zmn,Vm]=MM_EFIE_TE(Data,Nss,phi_i)
%% Definitions
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
phi_i   =   deg2rad(phi_i);
%% Fill Zmn and Vm matrices
Zmn   	=   zeros(Nss,Nss);
Vm   	=   zeros(Nss,1);
for m=1:Nss
    for n=1:Nss
        %% Data m
        xm      =   Data(m,1);
        ym      =   Data(m,2);
        xmm    	=   Data(m,3);
        ymm    	=   Data(m,4);
        xmp    	=   Data(m,5);
        ymp    	=   Data(m,6);
        [lmm,lmp,phi_mm,phi_mp]=rn(xm,ym,xmm,ymm,xmp,ymp);
        %% Data n
       	xn      =   Data(n,1);
        yn      =   Data(n,2);
        xnm    	=   Data(n,3);
        ynm    	=   Data(n,4);
        xnp    	=   Data(n,5);
        ynp    	=   Data(n,6);
        [lnm,lnp,phi_nm,phi_np]=rn(xn,yn,xnm,ynm,xnp,ynp);
        %% 
        if m==n
            Zmn(m,n)    =   Term_singular(lmp)+Term_singular(lmm)+Term_pm+Term_mp;
        elseif (xnp==xm)&&(ynp==ym)
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+Term_singular_s(lmm); 
        elseif (xnm==xm)&&(ynm==ym)
            Zmn(m,n)    =   Term_pp+Term_mm+Term_mp+Term_singular_s(lmp);
        else
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+Term_mp;  
        end
    end
  	Xmp     =   @(l) (1-l/lmp).*exp(+j*k*l*cos(phi_i-phi_mp));
  	Xmm     =   @(l) (1-l/lmm).*exp(-j*k*l*cos(phi_i-phi_mm));
   	Vm(m,1)	=   exp(j*k*(xm*cos(phi_i)+ym*sin(phi_i)))*...
        (sin(phi_mp-phi_i)*Quad(Xmp,0,lmp)+sin(phi_mm-phi_i)*Quad(Xmm,0,lmm));
end
%% Term ++
function[Termpp]=Term_pp()
dot_pp  =   cos(phi_mp)*cos(phi_np)+...
            sin(phi_mp)*sin(phi_np);
R_pp    =   @(l,l_) Rmn_pp(l,l_,xm,ym,xn,yn,phi_mp,phi_np);
g_pp  	=   @(l,l_) besselh(0,2,k*R_pp(l,l_))/(4*j);
Z1_pp   =   @(l,l_) j*k*eta*dot_pp*g_pp(l,l_).*(1-l/lmp).*(1-l_/lnp);
Z2_pp   =   @(l,l_) -j*(eta/k)*(1/(lmp*lnp))*g_pp(l,l_);
Termpp 	=   Quad2(Z1_pp,0,lmp,0,lnp)+Quad2(Z2_pp,0,lmp,0,lnp);
end
%% Term --
function[Termmm]=Term_mm()
dot_mm  =   cos(phi_mm)*cos(phi_nm)+...
            sin(phi_mm)*sin(phi_nm);
R_mm    =   @(l,l_) Rmn_mm(l,l_,xm,ym,xn,yn,phi_mm,phi_nm);
g_mm  	=   @(l,l_) besselh(0,2,k*R_mm(l,l_))/(4*j);
Z1_mm   =   @(l,l_) j*k*eta*dot_mm*g_mm(l,l_).*(1-l/lmm).*(1-l_/lnm);
Z2_mm   =   @(l,l_) -j*(eta/k)*(1/(lmm*lnm))*g_mm(l,l_);
Termmm 	=   Quad2(Z1_mm,0,lmm,0,lnm)+Quad2(Z2_mm,0,lmm,0,lnm);
end
%% Term +-
function[Termpm]=Term_pm()
dot_pm  =   cos(phi_mp)*cos(phi_nm)+...
            sin(phi_mp)*sin(phi_nm);
R_pm    =   @(l,l_) Rmn_pm(l,l_,xm,ym,xn,yn,phi_mp,phi_nm);
g_pm   	=   @(l,l_) besselh(0,2,k*R_pm(l,l_))/(4*j);
Z1_pm   =   @(l,l_) j*k*eta*dot_pm*g_pm(l,l_).*(1-l/lmp).*(1-l_/lnm);
Z2_pm   =   @(l,l_) -j*(eta/k)*(1/(lmp*lnm))*g_pm(l,l_);
Termpm 	=   Quad2(Z1_pm,0,lmp,0,lnm)-Quad2(Z2_pm,0,lmp,0,lnm);
end
%% Term -+
function[Termmp]=Term_mp()
dot_mp  =   cos(phi_mm)*cos(phi_np)+...
            sin(phi_mm)*sin(phi_np);
R_mp    =   @(l,l_) Rmn_mp(l,l_,xm,ym,xn,yn,phi_mm,phi_np);
g_mp  	=   @(l,l_) besselh(0,2,k*R_mp(l,l_))/(4*j);
Z1_mp   =   @(l,l_) j*k*eta*dot_mp*g_mp(l,l_).*(1-l/lmm).*(1-l_/lnp);
Z2_mp   =   @(l,l_) -j*(eta/k)*(1/(lmm*lnp))*g_mp(l,l_);
Termmp 	=   Quad2(Z1_mp,0,lmm,0,lnp)-Quad2(Z2_mp,0,lmm,0,lnp);
end
%% R ++
function[R]=Rmn_pp(l,l_,xm,ym,xn,yn,phi_m,phi_n)
Dmnx    =   xm-xn+l*cos(phi_m)-l_*cos(phi_n);
Dmny    =   ym-yn+l*sin(phi_m)-l_*sin(phi_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2);
end
%% R --
function[R]=Rmn_mm(l,l_,xm,ym,xn,yn,phi_m,phi_n)
Dmnx    =   xm-xn-l*cos(phi_m)+l_*cos(phi_n);
Dmny    =   ym-yn-l*sin(phi_m)+l_*sin(phi_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2);
end
%% R +-
function[R]=Rmn_pm(l,l_,xm,ym,xn,yn,phi_m,phi_n)
Dmnx    =   xm-xn+l*cos(phi_m)+l_*cos(phi_n);
Dmny    =   ym-yn+l*sin(phi_m)+l_*sin(phi_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2);
end
%% R -+
function[R]=Rmn_mp(l,l_,xm,ym,xn,yn,phi_m,phi_n)
Dmnx    =   xm-xn-l*cos(phi_m)-l_*cos(phi_n);
Dmny    =   ym-yn-l*sin(phi_m)-l_*sin(phi_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2);
end
%% Singular term ++ or --
function[I]=Term_singular(L)
gamma   =   1.781072417990198;
func1  	=   @(l) ((1-l/L).^2).*l.*(log(k*0.5*gamma*l)-1)+...
            (1-l/L).*(0.25*(l.^2)/L).*(2*log(k*0.5*gamma*l)-1);
I1    	=   L^2/4-j*(2/pi)*2*Quad(func1,0,L);
func2  	=   @(l) l.*(log(k*0.5*gamma*l)-1);
I2    	=   1-j*(2/pi)*(2/L^2)*Quad(func2,0,L);
I       =   j*k*eta*I1/(j*4)-j*(eta/k)*I2/(j*4);
end
%% Singular term +- or -+ 
function[I]=Term_singular_s(L)
gamma   =   1.781072417990198;
func1  	=   @(l) (1-l/L).*((L^2-l.^2).*log(0.5*k*gamma*(L-l))+...
                (l.^2).*log(0.5*k*gamma*l)-L*(L+2*l)/2);        
I1    	=   L^2/4-j*(1/pi)*Quad(func1,0,L)/L;
func2  	=   @(l) (L-l).*log(k*0.5*gamma*(L-l))+l.*log(k*0.5*gamma*l);
I2    	=   1-j*(2/pi)*((1/L^2)*Quad(func2,0,L)-1);
I       =   j*k*eta*I1/(j*4)+j*(eta/k)*I2/(j*4);
end
%%
end
%% Solve MFIE for TE
function[Zmn,Vm]=MM_MFIE_TE(Data,Nss,phi_i)
%% Definitions
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
phi_i   =   deg2rad(phi_i);
%% Fill Zmn and Vm matrices
Zmn   	=   zeros(Nss,Nss);
Vm   	=   zeros(Nss,1);
for m=1:Nss
    for n=1:Nss
        %% Data m
        xm      =   Data(m,1);
        ym      =   Data(m,2);
        xmm    	=   Data(m,3);
        ymm    	=   Data(m,4);
        xmp    	=   Data(m,5);
        ymp    	=   Data(m,6);
        [lmm,lmp,phi_mm,phi_mp]=rn(xm,ym,xmm,ymm,xmp,ymp);
        %% Data n
       	xn      =   Data(n,1);
        yn      =   Data(n,2);
        xnm    	=   Data(n,3);
        ynm    	=   Data(n,4);
        xnp    	=   Data(n,5);
        ynp    	=   Data(n,6);
        [lnm,lnp,phi_nm,phi_np]=rn(xn,yn,xnm,ynm,xnp,ynp);
        %%
        if m==n
            Zmn(m,n)    =   Term_pm+Term_mp+0.5*(lmp/3+lmm/3);
        elseif (xnp==xm)&&(ynp==ym)
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+0.5*lmm/6; 
        elseif (xnm==xm)&&(ynm==ym)
            Zmn(m,n)    =   Term_pp+Term_mm+Term_mp+0.5*lmp/6;
        else
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+Term_mp;  
        end  
    end
  	Xmp     =   @(l) (1-l/lmp).*exp(+j*k*l*cos(phi_i-phi_mp));
  	Xmm     =   @(l) (1-l/lmm).*exp(-j*k*l*cos(phi_i-phi_mm));
   	Vm(m,1)	=   exp(j*k*(xm*cos(phi_i)+ym*sin(phi_i)))*...
                (Quad(Xmp,0,lmp)+Quad(Xmm,0,lmm))/eta;
end
%% Term ++
function[Termpp]=Term_pp()
Dmnx    =   @(l,l_) xm-xn+l*cos(phi_mp)-l_*cos(phi_np);
Dmny    =   @(l,l_) ym-yn+l*sin(phi_mp)-l_*sin(phi_np);
RR     	=   @(l,l_) sqrt(abs(Dmnx(l,l_)).^2+abs(Dmny(l,l_)).^2);
dotR1  	=   @(l,l_) ((sin(phi_mp)*Dmnx(l,l_)-cos(phi_mp)*Dmny(l,l_))./RR(l,l_))*cos(phi_mp-phi_np);
dotR2  	=   @(l,l_) ((cos(phi_mp)*Dmnx(l,l_)+sin(phi_mp)*Dmny(l,l_))./RR(l,l_))*sin(phi_mp-phi_np);
g_pp  	=   @(l,l_) besselh(1,2,k*RR(l,l_));
Z1_pp   =   @(l,l_) (j*k/4)*(dotR1(l,l_)-dotR2(l,l_)).*g_pp(l,l_).*(1-l/lmp).*(1-l_/lnp);
Termpp 	=   Quad2(Z1_pp,0,lmp,0,lnp);
end
%% Term --
function[Termmm]=Term_mm()
Dmnx    =   @(l,l_) xm-xn-l*cos(phi_mm)+l_*cos(phi_nm);
Dmny    =   @(l,l_) ym-yn-l*sin(phi_mm)+l_*sin(phi_nm);
RR     	=   @(l,l_) sqrt(abs(Dmnx(l,l_)).^2+abs(Dmny(l,l_)).^2);
dotR1  	=   @(l,l_) ((sin(phi_mm)*Dmnx(l,l_)-cos(phi_mm)*Dmny(l,l_))./RR(l,l_))*cos(phi_mm-phi_nm);
dotR2  	=   @(l,l_) ((cos(phi_mm)*Dmnx(l,l_)+sin(phi_mm)*Dmny(l,l_))./RR(l,l_))*sin(phi_mm-phi_nm);
g_mm  	=   @(l,l_) besselh(1,2,k*RR(l,l_));
Z1_mm   =   @(l,l_) (j*k/4)*(dotR1(l,l_)-dotR2(l,l_)).*g_mm(l,l_).*(1-l/lmm).*(1-l_/lnm);
Termmm 	=   Quad2(Z1_mm,0,lmm,0,lnm);
end
%% Term +-
function[Termpm]=Term_pm()
Dmnx    =   @(l,l_) xm-xn+l*cos(phi_mp)+l_*cos(phi_nm);
Dmny    =   @(l,l_) ym-yn+l*sin(phi_mp)+l_*sin(phi_nm);
RR     	=   @(l,l_) sqrt(abs(Dmnx(l,l_)).^2+abs(Dmny(l,l_)).^2);
dotR1  	=   @(l,l_) ((sin(phi_mp)*Dmnx(l,l_)-cos(phi_mp)*Dmny(l,l_))./RR(l,l_))*cos(phi_mp-phi_nm);
dotR2  	=   @(l,l_) ((cos(phi_mp)*Dmnx(l,l_)+sin(phi_mp)*Dmny(l,l_))./RR(l,l_))*sin(phi_mp-phi_nm);
g_pm  	=   @(l,l_) besselh(1,2,k*RR(l,l_));
Z1_pm   =   @(l,l_) (j*k/4)*(dotR1(l,l_)-dotR2(l,l_)).*g_pm(l,l_).*(1-l/lmp).*(1-l_/lnm);
Termpm 	=   Quad2(Z1_pm,0,lmp,0,lnm);
end
%% Term -+
function[Termmp]=Term_mp()
Dmnx    =   @(l,l_) xm-xn-l*cos(phi_mm)-l_*cos(phi_np);
Dmny    =   @(l,l_) ym-yn-l*sin(phi_mm)-l_*sin(phi_np);
RR     	=   @(l,l_) sqrt(abs(Dmnx(l,l_)).^2+abs(Dmny(l,l_)).^2);
dotR1  	=   @(l,l_) ((sin(phi_mm)*Dmnx(l,l_)-cos(phi_mm)*Dmny(l,l_))./RR(l,l_))*cos(phi_mm-phi_np);
dotR2  	=   @(l,l_) ((cos(phi_mm)*Dmnx(l,l_)+sin(phi_mm)*Dmny(l,l_))./RR(l,l_))*sin(phi_mm-phi_np);
g_mp  	=   @(l,l_) besselh(1,2,k*RR(l,l_));
Z1_mp   =   @(l,l_) (j*k/4)*(dotR1(l,l_)-dotR2(l,l_)).*g_mp(l,l_).*(1-l/lmm).*(1-l_/lnp);
Termmp 	=   Quad2(Z1_mp,0,lmm,0,lnp);
end
end
%% TM solution
function[In,phi,sigma]=TM(phi_i,Data)
%%
alpha   =   0.5;
%%
[Nss,~] =   size(Data);
%% Call EFIE and MFIE
eta     =   120*pi;
tic;
fprintf('Solving the EFIE...\n');
[Zmn_EFIE,Vm_EFIE]=MM_EFIE_TM(Data,Nss,phi_i);
toc;
tic;
fprintf('Solving the MFIE...\n');
[Zmn_MFIE,Vm_MFIE]=MM_MFIE_TM(Data,Nss,phi_i);
toc;
%% Solve for current
fprintf('Solving the CFIE...\n');
Zmn     =   alpha*Zmn_EFIE+(1-alpha)*eta*Zmn_MFIE;
Vm      =   alpha*Vm_EFIE+(1-alpha)*eta*Vm_MFIE;
In      =   Zmn\Vm;
%% Find RCS
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
%%
[Nss,~] =   size(Data);
phi     =   linspace(-pi,pi,1e3);
sigma   =   0;
for n=1:Nss
    xn      =   Data(n,1);
    yn      =   Data(n,2);
    xnm    	=   Data(n,3);
    ynm    	=   Data(n,4);
    xnp    	=   Data(n,5);
    ynp    	=   Data(n,6);
    lnm     =   sqrt((xn-xnm)^2+(yn-ynm)^2);
 	lnp     =   sqrt((xn-xnp)^2+(yn-ynp)^2);
    Lxnm   	=   xn-xnm;
    Lynm   	=   yn-ynm;
    Lxnp   	=   xnp-xn;
    Lynp   	=   ynp-yn;
    kap_np  =   exp(j*k*(xnp*cos(phi)+ynp*sin(phi)));
    kap_nm  =   exp(j*k*(xnm*cos(phi)+ynm*sin(phi)));
    ap      =   Lxnp*cos(phi)+Lynp*sin(phi);
    am      =   Lxnm*cos(phi)+Lynm*sin(phi);
    Tp      =   (exp(-j*k*ap).*(1+j*k*ap)-1)./(k^2*ap.^2);
    Tm      =   (exp(+j*k*am).*(1-j*k*am)-1)./(k^2*am.^2);
    Tp(abs(ap)<1e-5)=0.5;
    Tm(abs(am)<1e-5)=0.5;
    %% Standard
    Term    =   lnp*kap_np.*Tp+lnm*kap_nm.*Tm;
    sigma   =   sigma+(eta/2)*In(n,1)*Term;
end
sigma   =   2*pi*abs(sigma).^2;
%% Plot RCS
figure()
plot(phi*180/pi,10*log10(sigma),'-k','LineWidth',1)
xlabel('$\varphi$ [deg]','Interpret','Latex','FontSize',14)
ylabel('$\sigma/\lambda$ [dB]','Interpret','Latex','FontSize',14)
title('TM','Interpret','Latex','FontSize',14)
set(gca,'TickLabel','Latex','FontSize',14)
xlim([-180 180])
%%
end
%% Solve EFIE for TM
function[Zmn,Vm]=MM_EFIE_TM(Data,Nss,phi_i)
%% Definitions
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
phi_i   =   deg2rad(phi_i);
%% Fill Zmn and Vm matrices
Zmn   	=   zeros(Nss,Nss);
Vm   	=   zeros(Nss,1);
for m=1:Nss
    for n=1:Nss
        %% Data m
        xm      =   Data(m,1);
        ym      =   Data(m,2);
        xmm    	=   Data(m,3);
        ymm    	=   Data(m,4);
        xmp    	=   Data(m,5);
        ymp    	=   Data(m,6);
        [lmm,lmp,phi_mm,phi_mp]=rn(xm,ym,xmm,ymm,xmp,ymp);
        %% Data n
       	xn      =   Data(n,1);
        yn      =   Data(n,2);
        xnm    	=   Data(n,3);
        ynm    	=   Data(n,4);
        xnp    	=   Data(n,5);
        ynp    	=   Data(n,6);
        [lnm,lnp,phi_nm,phi_np]=rn(xn,yn,xnm,ynm,xnp,ynp);
        %%
        if m==n
            Zmn(m,n)    =   Term_singular(lmp)+Term_singular(lmm)+Term_pm+Term_mp;
        elseif (xnp==xm)&&(ynp==ym)
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+Term_singular_s(lmm); 
        elseif (xnm==xm)&&(ynm==ym)
            Zmn(m,n)    =   Term_pp+Term_mm+Term_mp+Term_singular_s(lmp);
        else
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+Term_mp;  
        end
    end
  	Xmp     =   @(l) (1-l/lmp).*exp(+j*k*l*cos(phi_i-phi_mp));
  	Xmm     =   @(l) (1-l/lmm).*exp(-j*k*l*cos(phi_i-phi_mm));
   	Vm(m,1)	=   exp(j*k*(xm*cos(phi_i)+ym*sin(phi_i)))*...
                (Quad(Xmp,0,lmp)+Quad(Xmm,0,lmm));
end
%% Term ++
function[Termpp]=Term_pp()
R_pp    =   @(l,l_) Rmn_pp(l,l_,xm,ym,xn,yn,phi_mp,phi_np);
g_pp  	=   @(l,l_) besselh(0,2,k*R_pp(l,l_))/(4*j);
Z1_pp   =   @(l,l_) j*k*eta*g_pp(l,l_).*(1-l/lmp).*(1-l_/lnp);
Termpp 	=   Quad2(Z1_pp,0,lmp,0,lnp);
end
%% Term --
function[Termmm]=Term_mm()
R_mm    =   @(l,l_) Rmn_mm(l,l_,xm,ym,xn,yn,phi_mm,phi_nm);
g_mm  	=   @(l,l_) besselh(0,2,k*R_mm(l,l_))/(4*j);
Z1_mm   =   @(l,l_) j*k*eta*g_mm(l,l_).*(1-l/lmm).*(1-l_/lnm);
Termmm 	=   Quad2(Z1_mm,0,lmm,0,lnm);
end
%% Term +-
function[Termpm]=Term_pm()
R_pm    =   @(l,l_) Rmn_pm(l,l_,xm,ym,xn,yn,phi_mp,phi_nm);
g_pm   	=   @(l,l_) besselh(0,2,k*R_pm(l,l_))/(4*j);
Z1_pm   =   @(l,l_) j*k*eta*g_pm(l,l_).*(1-l/lmp).*(1-l_/lnm);
Termpm 	=   Quad2(Z1_pm,0,lmp,0,lnm);
end
%% Term -+
function[Termmp]=Term_mp()
R_mp    =   @(l,l_) Rmn_mp(l,l_,xm,ym,xn,yn,phi_mm,phi_np);
g_mp  	=   @(l,l_) besselh(0,2,k*R_mp(l,l_))/(4*j);
Z1_mp   =   @(l,l_) j*k*eta*g_mp(l,l_).*(1-l/lmm).*(1-l_/lnp);
Termmp 	=   Quad2(Z1_mp,0,lmm,0,lnp);
end
%% R ++
function[R]=Rmn_pp(l,l_,xm,ym,xn,yn,phi_m,phi_n)
Dmnx    =   xm-xn+l*cos(phi_m)-l_*cos(phi_n);
Dmny    =   ym-yn+l*sin(phi_m)-l_*sin(phi_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2);
end
%% R --
function[R]=Rmn_mm(l,l_,xm,ym,xn,yn,phi_m,phi_n)
Dmnx    =   xm-xn-l*cos(phi_m)+l_*cos(phi_n);
Dmny    =   ym-yn-l*sin(phi_m)+l_*sin(phi_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2);
end
%% R +-
function[R]=Rmn_pm(l,l_,xm,ym,xn,yn,phi_m,phi_n)
Dmnx    =   xm-xn+l*cos(phi_m)+l_*cos(phi_n);
Dmny    =   ym-yn+l*sin(phi_m)+l_*sin(phi_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2);
end
%% R -+
function[R]=Rmn_mp(l,l_,xm,ym,xn,yn,phi_m,phi_n)
Dmnx    =   xm-xn-l*cos(phi_m)-l_*cos(phi_n);
Dmny    =   ym-yn-l*sin(phi_m)-l_*sin(phi_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2);
end
%% Singular term ++ or --
function[I]=Term_singular(L)
gamma   =   1.781072417990198;
func1  	=   @(l) ((1-l/L).^2).*l.*(log(k*0.5*gamma*l)-1)+...
            (1-l/L).*(0.25*(l.^2)/L).*(2*log(k*0.5*gamma*l)-1);
I1    	=   L^2/4-j*(2/pi)*2*Quad(func1,0,L);
I       =   j*k*eta*I1/(j*4);
end
%% Singular term +- or -+ 
function[I]=Term_singular_s(L)
gamma   =   1.781072417990198;
func1  	=   @(l) (1-l/L).*((L^2-l.^2).*log(0.5*k*gamma*(L-l))+...
                (l.^2).*log(0.5*k*gamma*l)-L*(L+2*l)/2);        
I1    	=   L^2/4-j*(1/pi)*Quad(func1,0,L)/L;
I       =   j*k*eta*I1/(j*4);
end
%%
end
%% Solve MFIE for TM
function[Zmn,Vm]=MM_MFIE_TM(Data,Nss,phi_i)
%% Definitions
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
phi_i   =   deg2rad(phi_i);
%% Fill Zmn and Vm matrices
Zmn   	=   zeros(Nss,Nss);
Vm   	=   zeros(Nss,1);
for m=1:Nss
    for n=1:Nss
        %% Data m
        xm      =   Data(m,1);
        ym      =   Data(m,2);
        xmm    	=   Data(m,3);
        ymm    	=   Data(m,4);
        xmp    	=   Data(m,5);
        ymp    	=   Data(m,6);
        [lmm,lmp,phi_mm,phi_mp]=rn(xm,ym,xmm,ymm,xmp,ymp);
        %% Data n
       	xn      =   Data(n,1);
        yn      =   Data(n,2);
        xnm    	=   Data(n,3);
        ynm    	=   Data(n,4);
        xnp    	=   Data(n,5);
        ynp    	=   Data(n,6);
        [lnm,lnp,phi_nm,phi_np]=rn(xn,yn,xnm,ynm,xnp,ynp);
        %%
        if m==n
            Zmn(m,n)    =   Term_pm+Term_mp+0.5*(lmp/3+lmm/3);
        elseif (xnp==xm)&&(ynp==ym)
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+0.5*lmm/6; 
        elseif (xnm==xm)&&(ynm==ym)
            Zmn(m,n)    =   Term_pp+Term_mm+Term_mp+0.5*lmp/6;
        else
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+Term_mp;  
        end  
    end
  	Xmp     =   @(l) (1-l/lmp).*exp(+j*k*l*cos(phi_i-phi_mp));
  	Xmm     =   @(l) (1-l/lmm).*exp(-j*k*l*cos(phi_i-phi_mm));
   	Vm(m,1)	=   exp(j*k*(xm*cos(phi_i)+ym*sin(phi_i)))*...
        (sin(phi_mp-phi_i)*Quad(Xmp,0,lmp)+sin(phi_mm-phi_i)*Quad(Xmm,0,lmm))/eta;
end
%% Term ++
function[Termpp]=Term_pp()
Dmnx    =   @(l,l_) xm-xn+l*cos(phi_mp)-l_*cos(phi_np);
Dmny    =   @(l,l_) ym-yn+l*sin(phi_mp)-l_*sin(phi_np);
RR     	=   @(l,l_) sqrt(abs(Dmnx(l,l_)).^2+abs(Dmny(l,l_)).^2);
dotR1  	=   @(l,l_) ((sin(phi_mp)*Dmnx(l,l_)-cos(phi_mp)*Dmny(l,l_))./RR(l,l_));
dotR2  	=   @(l,l_) 0;
g_pp  	=   @(l,l_) besselh(1,2,k*RR(l,l_));
Z1_pp   =   @(l,l_) (j*k/4)*(dotR1(l,l_)-dotR2(l,l_)).*g_pp(l,l_).*(1-l/lmp).*(1-l_/lnp);
Termpp 	=   Quad2(Z1_pp,0,lmp,0,lnp);
end
%% Term --
function[Termmm]=Term_mm()
Dmnx    =   @(l,l_) xm-xn-l*cos(phi_mm)+l_*cos(phi_nm);
Dmny    =   @(l,l_) ym-yn-l*sin(phi_mm)+l_*sin(phi_nm);
RR     	=   @(l,l_) sqrt(abs(Dmnx(l,l_)).^2+abs(Dmny(l,l_)).^2);
dotR1  	=   @(l,l_) ((sin(phi_mm)*Dmnx(l,l_)-cos(phi_mm)*Dmny(l,l_))./RR(l,l_));
dotR2  	=   @(l,l_) 0;
g_mm  	=   @(l,l_) besselh(1,2,k*RR(l,l_));
Z1_mm   =   @(l,l_) (j*k/4)*(dotR1(l,l_)-dotR2(l,l_)).*g_mm(l,l_).*(1-l/lmm).*(1-l_/lnm);
Termmm 	=   Quad2(Z1_mm,0,lmm,0,lnm);
end
%% Term +-
function[Termpm]=Term_pm()
Dmnx    =   @(l,l_) xm-xn+l*cos(phi_mp)+l_*cos(phi_nm);
Dmny    =   @(l,l_) ym-yn+l*sin(phi_mp)+l_*sin(phi_nm);
RR     	=   @(l,l_) sqrt(abs(Dmnx(l,l_)).^2+abs(Dmny(l,l_)).^2);
dotR1  	=   @(l,l_) ((sin(phi_mp)*Dmnx(l,l_)-cos(phi_mp)*Dmny(l,l_))./RR(l,l_));
dotR2  	=   @(l,l_) 0;
g_pm  	=   @(l,l_) besselh(1,2,k*RR(l,l_));
Z1_pm   =   @(l,l_) (j*k/4)*(dotR1(l,l_)-dotR2(l,l_)).*g_pm(l,l_).*(1-l/lmp).*(1-l_/lnm);
Termpm 	=   Quad2(Z1_pm,0,lmp,0,lnm);
end
%% Term -+
function[Termmp]=Term_mp()
Dmnx    =   @(l,l_) xm-xn-l*cos(phi_mm)-l_*cos(phi_np);
Dmny    =   @(l,l_) ym-yn-l*sin(phi_mm)-l_*sin(phi_np);
RR     	=   @(l,l_) sqrt(abs(Dmnx(l,l_)).^2+abs(Dmny(l,l_)).^2);
dotR1  	=   @(l,l_) ((sin(phi_mm)*Dmnx(l,l_)-cos(phi_mm)*Dmny(l,l_))./RR(l,l_));
dotR2  	=   @(l,l_) 0;
g_mp  	=   @(l,l_) besselh(1,2,k*RR(l,l_));
Z1_mp   =   @(l,l_) (j*k/4)*(dotR1(l,l_)-dotR2(l,l_)).*g_mp(l,l_).*(1-l/lmm).*(1-l_/lnp);
Termmp 	=   Quad2(Z1_mp,0,lmm,0,lnp);
end
end
%% Find segments orientation definitions
function[lnm,lnp,phi_nm,phi_np]=rn(xn,yn,xnm,ynm,xnp,ynp)
%% Find lengths 
lnm     =   sqrt((xn-xnm)^2+(yn-ynm)^2);
lnp     =   sqrt((xnp-xn)^2+(ynp-yn)^2);
%% Correct angles
phi_nm	=   atan((yn-ynm)/(xn-xnm));
phi_nm(isnan(phi_nm))=0;
if (xn-xnm)<0
    phi_nm	=   phi_nm+pi;
end
phi_np	=   atan((ynp-yn)/(xnp-xn));
phi_np(isnan(phi_np))=0;
if (xnp-xn)<0
    phi_np	=   phi_np+pi;
end
end
%% Quad 1D
function[I]=Quad(func,a,b)
%% 16 Levels
x       =   [-0.0950125098376374 0.0950125098376374 -0.2816035507792589 ...
            0.2816035507792589 -0.4580167776572274 0.4580167776572274 ...
            -0.6178762444026438 0.6178762444026438 -0.7554044083550030 ...
            0.7554044083550030 -0.8656312023878318 0.8656312023878318 ...
        	-0.9445750230732326 0.9445750230732326 -0.9894009349916499 ...
             0.9894009349916499];
w       =   [0.1894506104550685 0.1894506104550685 0.1826034150449236 ...
            0.1826034150449236 0.1691565193950025 0.1691565193950025 ...
            0.1495959888165767 0.1495959888165767 0.1246289712555339...
         	0.1246289712555339 0.0951585116824928 0.0951585116824928 ...
            0.0622535239386479 0.0622535239386479 0.0271524594117541 ...
            0.0271524594117541];
%% Find I
hm    	=   (b-a)/2;
hp     	=   (b+a)/2;
f       =   func(hm*x+hp);
I       =   hm*(f*w');
end
%% Quad 2D
function[I]=Quad2(func,a1,b1,a2,b2)
%% 16 Levels
x       =   [-0.0950125098376374 0.0950125098376374 -0.2816035507792589 ...
            0.2816035507792589 -0.4580167776572274 0.4580167776572274 ...
            -0.6178762444026438 0.6178762444026438 -0.7554044083550030 ...
            0.7554044083550030 -0.8656312023878318 0.8656312023878318 ...
        	-0.9445750230732326 0.9445750230732326 -0.9894009349916499 ...
             0.9894009349916499];
w       =   [0.1894506104550685 0.1894506104550685 0.1826034150449236 ...
            0.1826034150449236 0.1691565193950025 0.1691565193950025 ...
            0.1495959888165767 0.1495959888165767 0.1246289712555339...
         	0.1246289712555339 0.0951585116824928 0.0951585116824928 ...
            0.0622535239386479 0.0622535239386479 0.0271524594117541 ...
            0.0271524594117541];
%% Find I
[x,y] 	=   meshgrid(x,x);
hm1   	=   (b1-a1)/2;
hp1   	=   (b1+a1)/2;
hm2   	=   (b2-a2)/2;
hp2   	=   (b2+a2)/2;
f       =   func(hm1*x+hp1,hm2*y+hp2);
I       =   hm1*hm2*(w*f*w');
end
%%
















