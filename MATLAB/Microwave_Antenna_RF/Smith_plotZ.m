function Smith_plotZ(ax,Z,Z0)
% add text showing ZL R+jX data, Smith Chart generated externally

gamma=z2gamma(Z,Z0);  
plot(ax,real(gamma),imag(gamma),'ro','LineWidth',1.5)

if imag(Z)<0
sign1='-';
else
    sign1='+';
end
hold all;plot(ax,real(gamma),imag(gamma),'ro','LineWidth',1.5);
str1=['ZL =' num2str(real(Z)) sign1 'j' num2str(abs(imag(Z))) ' \rightarrow'];
text(ax,real(gamma),imag(gamma)+.01,str1,'Color','blue','FontSize',20,'HorizontalAlignment','right','VerticalAlignment','middle');

end