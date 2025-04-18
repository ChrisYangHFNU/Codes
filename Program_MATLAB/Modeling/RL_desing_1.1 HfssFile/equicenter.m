function [XB, YB] = equicenter(XA,YA,XC,YC,P)

f =  @(x) sqrt((x-XA)^2+(polyval(P,x)-YA)^2)-sqrt((XC-x)^2+(YC-polyval(P,x))^2);XB=fzero(f,XC+XA/2);
YB=(polyval(P,XB));