function[I,Zin,theta_E,E_theta_E,phi_A,E_theta_A]=DipoleMM(L,a,NN)
%%
[Data,M]=VerticalDipole(NN,L);
[Nss,~] =   size(Data);
%%
tic;
I       =   MM(Data,M,Nss,a);
Zin     =   1/I(M);
toc;
if imag(1/I(M))>0
    fprintf('Zin\t=\t%0.2f\t+j\t%0.2f\n',real(1/I(M)),imag(1/I(M)));
else
    fprintf('Zin\t=\t%0.2f\t-j\t%0.2f\n',real(1/I(M)),-imag(1/I(M)));
end
I       =   [0;I;0];
%%
L_    	=   zeros(1,Nss);
xn      =   Data(1,1);
yn      =   Data(1,2);
zn      =   Data(1,3);
xnm    	=   Data(1,4);
ynm    	=   Data(1,5);
znm    	=   Data(1,6);
xnp    	=   Data(1,7);
ynp    	=   Data(1,8);
znp    	=   Data(1,9);
[lnm,~,~,~,~,~]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp);
sum     =   lnm;
for i=1:Nss
    xn      =   Data(i,1);
    yn      =   Data(i,2);
    zn      =   Data(i,3);
    xnm    	=   Data(i,4);
    ynm    	=   Data(i,5);
    znm    	=   Data(i,6);
    xnp    	=   Data(i,7);
    ynp    	=   Data(i,8);
    znp    	=   Data(i,9);
    [~,lnp,~,~,~,~]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp);
    L_(1,i) =   sum;
    sum     =   sum+lnp;
end
L_      =   [0 L_ sum];
%%
mA      =   1e-3;
figure()
hold on
plot(L_-L/2,real(I)/mA,'-k','LineWidth',1)
plot(L_-L/2,imag(I)/mA,'--k','LineWidth',1)
hold off
xlabel('$l/\lambda$','Interpret','Latex')
ylabel('$I_{n}(l)$ [mA]','Interpret','Latex')
legend('Re','Im','Interpreter','Latex','Location','NorthEast')
set(gca,'TickLabel','Latex')
%%
mA      =   1e-3;
[N,~]   =   size(I);
figure()
hold on
for i=2:N-1
    plot([L_(i-1) L_(i)]-L/2,[0 abs(I(i))]/mA,'-k','LineWidth',0.1)
    plot([L_(i) L_(i+1)]-L/2,[abs(I(i)) 0]/mA,'-k','LineWidth',0.1)
end
plot(L_-L/2,abs(I)/mA,'-ob','MarkerFace','k','MarkerSize',2,'LineWidth',1)
hold off
xlabel('$l/\lambda$','Interpret','Latex')
ylabel('$|I_{n}(l)|$ [mA]','Interpret','Latex')
set(gca,'TickLabel','Latex')
%%
I           =   I(2:end-1);
%%
[E_theta_A,~,~,phi_A]=FarField_PhiCut(90,Data,I);
E_theta_A 	=   E_theta_A./max(abs(E_theta_A));
figure()
PolardB(phi_A,20*log10(abs(E_theta_A)),[-25 0],6,'-k')
title('$20\log_{10}|E_{\theta}|$ [dB]','Interpret','Latex')
pax = gca;
pax.RAxisLocation   =   90;
pax.ThetaDir = 'CounterClockwise';
pax.ThetaZeroLocation = 'right';
%%
[E_theta_E,~,~,theta_E]=FarField_ThetaCut(0,Data,I);
E_theta_E 	=   E_theta_E./max(abs(E_theta_E));
figure()
PolardB(theta_E,20*log10(abs(E_theta_E)),[-25 0],6,'-k')
title('$20\log_{10}|E_{\theta}|$ [dB]','Interpret','Latex')
%%
end
%%
function[In]=MM(Data,M,Nss,a)
%%
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
%%
Zmn   	=   zeros(Nss,Nss);
Vm   	=   zeros(Nss,1);
for m=1:Nss
    for n=1:Nss
        %% Data m
        xm      =   Data(m,1);
        ym      =   Data(m,2);
        zm      =   Data(m,3);
        xmm    	=   Data(m,4);
        ymm    	=   Data(m,5);
        zmm    	=   Data(m,6);
        xmp    	=   Data(m,7);
        ymp    	=   Data(m,8);
        zmp    	=   Data(m,9);
        [lmm,lmp,theta_mm,phi_mm,theta_mp,phi_mp]=rn(xm,ym,zm,xmm,ymm,zmm,xmp,ymp,zmp);
        %% Data n
       	xn      =   Data(n,1);
        yn      =   Data(n,2);
        zn      =   Data(n,3);
        xnm    	=   Data(n,4);
        ynm    	=   Data(n,5);
        znm    	=   Data(n,6);
        xnp    	=   Data(n,7);
        ynp    	=   Data(n,8);
        znp    	=   Data(n,9);
        [lnm,lnp,theta_nm,phi_nm,theta_np,phi_np]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp);
        %%
        if m==n
            Zmn(m,n)    =   SingularTerm(lmp)+SingularTerm(lmm)+Term_pm+Term_mp;
        elseif m==n+1 
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+SingularTerm2(lmm); 
        elseif n==m+1 
            Zmn(m,n)    =   Term_pp+Term_mm+Term_mp+SingularTerm2(lmp);
        else
            Zmn(m,n)    =   Term_pp+Term_mm+Term_pm+Term_mp;  
        end  
    end
    if m==M
        Vm(m,1)	=   1;
    end
end
%%
In      =   Zmn\Vm;
%% Term ++
function[Termpp]=Term_pp()
dot_pp  =   sin(theta_mp)*cos(phi_mp)*sin(theta_np)*cos(phi_np)+...
            sin(theta_mp)*sin(phi_mp)*sin(theta_np)*sin(phi_np)+...
          	cos(theta_mp)*cos(theta_np);
R_pp    =   @(l,l_) Rmn_pp(l,l_,xm,ym,zm,xn,yn,zn,theta_mp,phi_mp,theta_np,phi_np);
g_pp  	=   @(l,l_) exp(-j*k*R_pp(l,l_))./(4*pi*R_pp(l,l_));
Z1_pp   =   @(l,l_) j*k*eta*dot_pp*g_pp(l,l_).*(1-l/lmp).*(1-l_/lnp);
Z2_pp   =   @(l,l_) -j*(eta/k)*(1/(lmp*lnp))*g_pp(l,l_);
Termpp 	=   Quad2(Z1_pp,0,lmp,0,lnp)+Quad2(Z2_pp,0,lmp,0,lnp);
end
%% Term --
function[Termmm]=Term_mm()
dot_mm  =   sin(theta_mm)*cos(phi_mm)*sin(theta_nm)*cos(phi_nm)+...
            sin(theta_mm)*sin(phi_mm)*sin(theta_nm)*sin(phi_nm)+...
            cos(theta_mm)*cos(theta_nm);
R_mm    =   @(l,l_) Rmn_mm(l,l_,xm,ym,zm,xn,yn,zn,theta_mm,phi_mm,theta_nm,phi_nm);
g_mm  	=   @(l,l_) exp(-j*k*R_mm(l,l_))./(4*pi*R_mm(l,l_));
Z1_mm   =   @(l,l_) j*k*eta*dot_mm*g_mm(l,l_).*(1-l/lmm).*(1-l_/lnm);
Z2_mm   =   @(l,l_) -j*(eta/k)*(1/(lmm*lnm))*g_mm(l,l_);
Termmm 	=   Quad2(Z1_mm,0,lmm,0,lnm)+Quad2(Z2_mm,0,lmm,0,lnm);
end
%% Term +-
function[Termpm]=Term_pm()
dot_pm  =   sin(theta_mp)*cos(phi_mp)*sin(theta_nm)*cos(phi_nm)+...
            sin(theta_mp)*sin(phi_mp)*sin(theta_nm)*sin(phi_nm)+...
            cos(theta_mp)*cos(theta_nm);
R_pm    =   @(l,l_) Rmn_pm(l,l_,xm,ym,zm,xn,yn,zn,theta_mp,phi_mp,theta_nm,phi_nm);
g_pm   	=   @(l,l_) exp(-j*k*R_pm(l,l_))./(4*pi*R_pm(l,l_));
Z1_pm   =   @(l,l_) j*k*eta*dot_pm*g_pm(l,l_).*(1-l/lmp).*(1-l_/lnm);
Z2_pm   =   @(l,l_) -j*(eta/k)*(1/(lmp*lnm))*g_pm(l,l_);
Termpm 	=   Quad2(Z1_pm,0,lmp,0,lnm)-Quad2(Z2_pm,0,lmp,0,lnm);
end
%% Term -+
function[Termmp]=Term_mp()
dot_mp  =   sin(theta_mm)*cos(phi_mm)*sin(theta_np)*cos(phi_np)+...
            sin(theta_mm)*sin(phi_mm)*sin(theta_np)*sin(phi_np)+...
            cos(theta_mm)*cos(theta_np);
R_mp    =   @(l,l_) Rmn_mp(l,l_,xm,ym,zm,xn,yn,zn,theta_mm,phi_mm,theta_np,phi_np);
g_mp  	=   @(l,l_) exp(-j*k*R_mp(l,l_))./(4*pi*R_mp(l,l_));
Z1_mp   =   @(l,l_) j*k*eta*dot_mp*g_mp(l,l_).*(1-l/lmm).*(1-l_/lnp);
Z2_mp   =   @(l,l_) -j*(eta/k)*(1/(lmm*lnp))*g_mp(l,l_);
Termmp 	=   Quad2(Z1_mp,0,lmm,0,lnp)-Quad2(Z2_mp,0,lmm,0,lnp);
end
%% R ++
function[R]=Rmn_pp(l,l_,xm,ym,zm,xn,yn,zn,theta_m,phi_m,theta_n,phi_n)
Dmnx    =   xm-xn+l*sin(theta_m)*cos(phi_m)-l_*sin(theta_n)*cos(phi_n);
Dmny    =   ym-yn+l*sin(theta_m)*sin(phi_m)-l_*sin(theta_n)*sin(phi_n);
Dmnz    =   zm-zn+l*cos(theta_m)-l_*cos(theta_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2+abs(Dmnz).^2+a^2);
end
%% R --
function[R]=Rmn_mm(l,l_,xm,ym,zm,xn,yn,zn,theta_m,phi_m,theta_n,phi_n)
Dmnx    =   xm-xn-l*sin(theta_m)*cos(phi_m)+l_*sin(theta_n)*cos(phi_n);
Dmny    =   ym-yn-l*sin(theta_m)*sin(phi_m)+l_*sin(theta_n)*sin(phi_n);
Dmnz    =   zm-zn-l*cos(theta_m)+l_*cos(theta_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2+abs(Dmnz).^2+a^2);
end
%% R +-
function[R]=Rmn_pm(l,l_,xm,ym,zm,xn,yn,zn,theta_m,phi_m,theta_n,phi_n)
Dmnx    =   xm-xn+l*sin(theta_m)*cos(phi_m)+l_*sin(theta_n)*cos(phi_n);
Dmny    =   ym-yn+l*sin(theta_m)*sin(phi_m)+l_*sin(theta_n)*sin(phi_n);
Dmnz    =   zm-zn+l*cos(theta_m)+l_*cos(theta_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2+abs(Dmnz).^2+a^2);
end
%% R -+
function[R]=Rmn_mp(l,l_,xm,ym,zm,xn,yn,zn,theta_m,phi_m,theta_n,phi_n)
Dmnx    =   xm-xn-l*sin(theta_m)*cos(phi_m)-l_*sin(theta_n)*cos(phi_n);
Dmny    =   ym-yn-l*sin(theta_m)*sin(phi_m)-l_*sin(theta_n)*sin(phi_n);
Dmnz    =   zm-zn-l*cos(theta_m)-l_*cos(theta_n);
R       =   sqrt(abs(Dmnx).^2+abs(Dmny).^2+abs(Dmnz).^2+a^2);
end
%% Singular term ++ or --
function[I]=SingularTerm(L)
func1  	=   @(l) (1-l/L).*(sqrt(l.^2+a^2)/L+(1-l/L).*0.5.*log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l)));
func2  	=   @(l) (1-l/L).*(-sqrt((L-l).^2+a^2)/L+(1-l/L).*0.5.*log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l))));
I1    	=   -j*k*L^2/(16*pi)+Quad(func1,0,L)/(4*pi)+Quad(func2,0,L)/(4*pi);
func3  	=   @(l) 0.5*log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l));
func4  	=   @(l) 0.5*log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l)));
I2    	=   -j*k*L^2/(4*pi)+Quad(func3,0,L)/(4*pi)+Quad(func4,0,L)/(4*pi);
I       =   j*k*eta*I1-j*(eta/k)*(1/(L^2))*I2;
end
%% Singular term +- or -+
function[I]=SingularTerm2(L)
func1  	=   @(l) (1/L)*(1-l/L).*(sqrt((L-l).^2+a^2)-sqrt(l.^2+a^2));
func2  	=   @(l) (0.5/L)*(1-l/L).*l.*(log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l))+log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l))));
I1    	=   -j*k*L^2/(16*pi)+Quad(func1,0,L)/(4*pi)+Quad(func2,0,L)/(4*pi);
func3  	=   @(l) 0.5*log((sqrt(l.^2+a^2)+l)./(sqrt(l.^2+a^2)-l));
func4  	=   @(l) 0.5*log((sqrt((L-l).^2+a^2)+(L-l))./(sqrt((L-l).^2+a^2)-(L-l)));
I2    	=   -j*k*L^2/(4*pi)+Quad(func3,0,L)/(4*pi)+Quad(func4,0,L)/(4*pi);
I       =   j*k*eta*I1+j*(eta/k)*(1/(L^2))*I2;
end
end
%%
function[lnm,lnp,theta_nm,phi_nm,theta_np,phi_np]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp)
%%
lnm     =   sqrt((xn-xnm)^2+(yn-ynm)^2+(zn-znm)^2);
lnp     =   sqrt((xnp-xn)^2+(ynp-yn)^2+(znp-zn)^2);
%%
theta_nm    =   acos((zn-znm)/lnm);
theta_nm(isnan(theta_nm))=0;
theta_np    =   acos((znp-zn)/lnp);
theta_np(isnan(theta_np))=0;
%%
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
%%
function[Data,M]=VerticalDipole(Ns,L)
%%
Ns      =   round(Ns,0);
if mod(Ns,2)~=0
    Ns      =   Ns+1;
end
%%
dz      =   L/Ns;
Data    =   zeros(Ns-1,9);
%%
for i=1:Ns-1
    xn      =   0;
    yn      =   0;
    zn      =   i*dz;
    xnm    	=   0;
    ynm    	=   0;
    znm    	=   (i-1)*dz;
    xnp    	=   0;
    ynp    	=   0;
    znp    	=   (i+1)*dz;
    Data(i,:)   =   [xn yn zn xnm ynm znm xnp ynp znp];
end
%%
Data(:,3)   =   Data(:,3)-L/2;  
Data(:,6)   =   Data(:,6)-L/2;
Data(:,9)   =   Data(:,9)-L/2;
%%
M       =   Ns/2;
end
%%
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
%%
hm    	=   (b-a)/2;
hp     	=   (b+a)/2;
f       =   func(hm*x+hp);
I       =   hm*(f*w');
end
%%
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
%%
[x,y] 	=   meshgrid(x,x);
%%
hm1   	=   (b1-a1)/2;
hp1   	=   (b1+a1)/2;
hm2   	=   (b2-a2)/2;
hp2   	=   (b2+a2)/2;
f       =   func(hm1*x+hp1,hm2*y+hp2);
I       =   hm1*hm2*(w*f*w');
end
%%
function[E_theta,E_phi,E,phi]=FarField_PhiCut(theta,Data,I)
%%
j       =   sqrt(-1);
Ns      =   1e3;
k       =   2*pi;
%%
theta 	=   deg2rad(theta);
phi     =   linspace(0,2*pi,Ns);
%%
[N,~]   =   size(I);
E_phi       =   0;
E_theta     =   0;
for i=1:N-2
    xn      =   Data(i,1);
    yn      =   Data(i,2);
    zn      =   Data(i,3);
    xnm    	=   Data(i,4);
    ynm    	=   Data(i,5);
    znm    	=   Data(i,6);
    xnp    	=   Data(i,7);
    ynp    	=   Data(i,8);
    znp    	=   Data(i,9);
    [lnm,lnp,theta_nm,phi_nm,theta_np,phi_np]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp);
    k_r     =   k*(xn*sin(theta)*cos(phi)+yn*sin(theta)*sin(phi)+zn*cos(theta));
    dot_p   =   lnp*sin(theta_np)*cos(phi_np)*cos(theta)*cos(phi)+...
                lnp*sin(theta_np)*sin(phi_np)*cos(theta)*sin(phi)-...
                lnp*cos(theta_np)*sin(theta);
    dot_m   =   lnm*sin(theta_nm)*cos(phi_nm)*cos(theta)*cos(phi)+...
                lnm*sin(theta_nm)*sin(phi_nm)*cos(theta)*sin(phi)-...
                lnm*cos(theta_nm)*sin(theta);
    E_theta =   E_theta+I(i+1,1)*(dot_p+dot_m).*exp(j*k_r);
    dot_p   =   -lnp*sin(theta_np)*cos(phi_np)*sin(phi)+...
                lnp*sin(theta_np)*sin(phi_np)*cos(phi);
    dot_m   =   -lnm*sin(theta_nm)*cos(phi_nm)*sin(phi)+...
                lnm*sin(theta_nm)*sin(phi_nm)*cos(phi);
    E_phi  	=   E_phi+I(i+1,1)*(dot_p+dot_m).*exp(j*k_r);
end
%%
E_theta(isnan(E_theta))=0;
E_theta(abs(E_theta)<1e-10)=0;
E_phi(isnan(E_phi))=0;
E_phi(abs(E_phi)<1e-10)=0;
E   	=   sqrt(abs(E_theta).^2+abs(E_phi).^2);
end
%%
function[E_theta,E_phi,E,theta]=FarField_ThetaCut(phi,Data,I)
%%
j       =   sqrt(-1);
Ns      =   1e3;
k       =   2*pi;
%%
phi     =   deg2rad(phi);
theta   =   linspace(-pi,pi,Ns);
%%
[N,~]   =   size(I);
E_phi       =   0;
E_theta     =   0;
for i=1:N-2
    xn      =   Data(i,1);
    yn      =   Data(i,2);
    zn      =   Data(i,3);
    xnm    	=   Data(i,4);
    ynm    	=   Data(i,5);
    znm    	=   Data(i,6);
    xnp    	=   Data(i,7);
    ynp    	=   Data(i,8);
    znp    	=   Data(i,9);
    [lnm,lnp,theta_nm,phi_nm,theta_np,phi_np]=rn(xn,yn,zn,xnm,ynm,znm,xnp,ynp,znp);
    k_r     =   k*(xn*sin(theta)*cos(phi)+yn*sin(theta)*sin(phi)+zn*cos(theta));
    dot_p   =   lnp*sin(theta_np)*cos(phi_np)*cos(theta)*cos(phi)+...
                lnp*sin(theta_np)*sin(phi_np)*cos(theta)*sin(phi)-...
                lnp*cos(theta_np)*sin(theta);
    dot_m   =   lnm*sin(theta_nm)*cos(phi_nm)*cos(theta)*cos(phi)+...
                lnm*sin(theta_nm)*sin(phi_nm)*cos(theta)*sin(phi)-...
                lnm*cos(theta_nm)*sin(theta);
    E_theta =   E_theta+I(i+1,1)*(dot_p+dot_m).*exp(j*k_r);
    dot_p   =   -lnp*sin(theta_np)*cos(phi_np)*sin(phi)+...
                lnp*sin(theta_np)*sin(phi_np)*cos(phi);
    dot_m   =   -lnm*sin(theta_nm)*cos(phi_nm)*sin(phi)+...
                lnm*sin(theta_nm)*sin(phi_nm)*cos(phi);
    E_phi  	=   E_phi+I(i+1,1)*(dot_p+dot_m).*exp(j*k_r);
end
%%
E_theta(isnan(E_theta))=0;
E_theta(abs(E_theta)<1e-10)=0;
E_phi(isnan(E_phi))=0;
E_phi(abs(E_phi)<1e-10)=0;
E       =   sqrt(abs(E_theta).^2+abs(E_phi).^2);
end
%%
function[]=PolardB(theta,data,range,Ntick,Line,option)
%%
switch nargin
    case 2
        range 			=	[-40 0]; % default
        Ntick 			= 	5;
		Line 			= 	'k';
		option 			= 	1;
	case 3
        Ntick 			= 	5;
		Line 			= 	'k';
		option 			= 	1;
	case 4
		Line 			= 	'k';
		option 			= 	1;
    case 5
		option 			= 	1;
end
data(isnan(data))       =   min(range);
data(data < min(range))	=   min(range);
polarplot(theta,data,Line,'LineWidth',1)
rlim(range)
%% 
if option == 1
thetaticks(0:15:345)
thetaticklabels({'$0^\circ$',' ','$30^\circ$',' ','$60^\circ$',' ','$90^\circ$',' ','$120^\circ$',' ','$150^\circ$',' ','$180^\circ$',' ','$-150^\circ$',' ','$-120^\circ$',' ','$-90^\circ$',' ','$-60^\circ$',' ','$-30^\circ$'})
end
%% 
if option == 2
thetalim([-90 90])
thetaticks(-90:15:90)
thetaticklabels({'$-90^\circ$',' ','$-60^\circ$',' ','$-30^\circ$',' ','$0^\circ$',' ','$30^\circ$',' ','$60^\circ$',' ','$90^\circ$'})
end
%%
range_                  =   linspace(min(range),max(range),Ntick); 
rticks(range_)
Labels                  =   cell(1,Ntick);
for i=1:Ntick
    %% Default
  	Labels{1,i}    	=   num2str(round(range_(i),2));
end
rticklabels(Labels)
set(gca,'TickLabel','Latex')
%%
pax = gca;
pax.RAxisLocation   =   0;
pax.ThetaDir = 'clockwise';
pax.ThetaZeroLocation = 'top';
end
%%

