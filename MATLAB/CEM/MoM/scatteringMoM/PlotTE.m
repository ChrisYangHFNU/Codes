function[H]=PlotTE(phi_i,In_TE,Data_TE,r_)
Ns      =   600;
j       =   sqrt(-1);
k       =   2*pi;
eta     =   120*pi;
N       =   length(In_TE);
[x,y]   =   meshgrid(linspace(-r_,r_,Ns),linspace(-r_,r_,Ns));
Hs      =   0;
fprintf('\nPlotting TE fields...\n');
for n=1:N
    xn      =   Data_TE(n,1);
    yn      =   Data_TE(n,2);
    xnm    	=   Data_TE(n,3);
    ynm    	=   Data_TE(n,4);
    xnp    	=   Data_TE(n,5);
    ynp    	=   Data_TE(n,6);
    Lxnm   	=   xn-xnm;
    Lynm   	=   yn-ynm;
    Lxnp   	=   xnp-xn;
    Lynp   	=   ynp-yn;
    R       =   sqrt((x-xn).^2+(y-yn).^2);
    val     =   0.1;
    R       =   sqrt((x-xn).^2+(y-yn).^2).*(R>val)+sqrt((x-xn).^2+(y-yn).^2+1E-4./R).*(R<=val);
    Rx      =   (x-xn)./R;
    Ry      =   (y-yn)./R;
    Cross   =   Rx*(Lynm+Lynp)-Ry*(Lxnm+Lxnp);
    Hs      =   Hs-(j*k/4)*In_TE(n,1)*(Cross/2).*besselh(1,2,k*R);
end
phi_i   =   phi_i*pi/180;
Hi      =   exp(j*k*(x*cos(phi_i)+y*sin(phi_i)))/eta;
H       =   Hi+Hs;
%%
figure()
pcolor(x,y,abs(H)*eta)
hold on
Plot(Data_TE)
hold off
xlabel('$x/\lambda$','Interpret','Latex','FontSize',14)
ylabel('$y/\lambda$','Interpret','Latex','FontSize',14)
title('$\eta|H|$ for TE','Interpret','Latex','FontSize',14)
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
    plot([Data(i,1) Data(i,3)],[Data(i,2) Data(i,4)],'-k','LineWidth',1)
    plot([Data(i,1) Data(i,5)],[Data(i,2) Data(i,6)],'-k','LineWidth',1)
end
end   
%%