function a=Smith_plotOpRefLine2PhaseCircle(ax,Z,Z0)
% plot reference line from origin oppsite to ZL to phase reference circle
% it plots ref line through YL

gamma_L=(Z-Z0)/(Z+Z0)

a=angle(gamma_L)
plot(ax,[0 -real(exp(1j*a))],[0 -imag(exp(1j*a))],'Color',[0 0.7 0],'LineWidth',1)           

end

