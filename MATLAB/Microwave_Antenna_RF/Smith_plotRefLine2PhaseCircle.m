function Smith_plotRefLine2PhaseCircle(ax,Z,Z0)
% plot reference line from origin through ZL to phase reference circle

gamma_L=z2gamma(Z,Z0);                                
a=atan(imag(gamma_L)/real(gamma_L));
plot(ax,[0 cos(a)],[0 sin(a)],'Color',[0 0.7 0],'LineWidth',1)           

end

