function [Pinfo,Xite,Wtrans]=strucLens(Dspacing,CenterF,WidthTrans,SidCurvature,LB,LD,LR,Nbh,Ndh,Nrh,e,gama,ap,g,phim,r,Nbeam,Ndum,Nside,Nrcv,CtrN1,CtrN2)
Lamda = 0.3/CenterF;
Np = Nbeams*(Ndum+1);
Nr = 10*Nrcv;
w = WidthTrans;
W = ones(1,Np)*w;
L = ones(1,Np)*LB*Lamda;
for ii=1:Ndum+1:Np
    L(ii+1:ii+Ndum) = LD*Lamda;
end
Ne = Nbeam*(Ndum+1)*10;
[xb,yb,xr,yr,Wtran] = RotmBRW(eps,e,gama,r,Nr,ap,g,phim,Ne);
F1 = (Nrcv-1)*0.6*r/CenterF*Dspacing;
xb = F1*xb;
yb = F1*yb;
xr = F1*xr;
yr = F1*yr;
Wtran = F1*Wtran;
index = (1:Nrcv)*10-9;
Wtrans = Wtran(index);
Wtrans = [fliplr(Wtrans),Wtrans(2:end)];
yb(1) = 10e-10;
yr(1)=10e-9;

K = zeros(1,Np);
for jj=1:Np
    if CtrN1==0
        K(jj) = yb(10*jj-9)/xb(10*jj-9);
    else
        if CtrN1==1
            K(jj) = -(xb(10*jj-8)-xb(10*jj-9))/(yb(10*jj-8)-yb(10*jj-9));
        end
    end
end
Edge = zeros(2,2*Np);
for jj=1:Np
    Pout = PCarmcal([xb(10*jj-9),yb(10*jj-9),K(jj),L(jj),W(jj),0,xb(1)/2,0]);
    Edge(1,(2*jj-1):(2*jj)) = Pout(1:2);
    Edge(2,(2*jj-1):(2*jj)) = Pout(3:4);
end
