function[E]=PlotTM(phi_i,In_TM,Data_TM,r_)
Ns      =   600;
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
N       =   length(In_TM);
[x,y]   =   meshgrid(linspace(-r_,r_,Ns),linspace(-r_,r_,Ns));
Es      =   0;
fprintf('\nPlotting TM fields...\n');
for n=1:N
    xn      =   Data_TM(n,1);
    yn      =   Data_TM(n,2);
    xnm    	=   Data_TM(n,3);
    ynm    	=   Data_TM(n,4);
    xnp    	=   Data_TM(n,5);
    ynp    	=   Data_TM(n,6);
    lnm     =   sqrt((xn-xnm)^2+(yn-ynm)^2);
 	lnp     =   sqrt((xn-xnp)^2+(yn-ynp)^2);
    R       =   sqrt((x-xn).^2+(y-yn).^2);
    val     =   0.1;
    R       =   sqrt((x-xn).^2+(y-yn).^2).*(R>val)+sqrt((x-xn).^2+(y-yn).^2+1E-4./R).*(R<=val);
    Es      =   Es-(k*eta/4)*((lnm+lnp)/2)*In_TM(n,1)*besselh(0,2,k*R);
end
phi_i   =   phi_i*pi/180;
Ei      =   exp(j*k*(x*cos(phi_i)+y*sin(phi_i)));
E       =   Ei+Es;
%%
figure()
pcolor(x,y,abs(E))
hold on
Plot(Data_TM)
hold off
xlabel('$x/\lambda$','Interpret','Latex','FontSize',14)
ylabel('$y/\lambda$','Interpret','Latex','FontSize',14)
title('$|E|$ for TM','Interpret','Latex','FontSize',14)
shading flat 
colormap jet 
colorbar 
axis equal
set(gca,'TickLabel','Latex','FontSize',15)
set(colorbar,'TickLabelInterpret','Latex','FontSize',14)
axis([-1 +1 -1 +1]*r_)
caxis([0 2])
end
%%
function[]=Plot(Data)
[N,~]   =   size(Data);
for i=1:N
    plot(Data(i,1),Data(i,2),'.k')
    plot([Data(i,1) Data(i,3)],[Data(i,2) Data(i,4)],'-k','LineWidth',1)
    plot([Data(i,1) Data(i,5)],[Data(i,2) Data(i,6)],'-k','LineWidth',1)
end
end   
%%