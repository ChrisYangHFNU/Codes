a = 0.35;
[X,Y] = meshgrid(0:0.1:1);
[theta,r] = cart2pol(Y,X);
rout = r; rin = r;
rin(find(rin>a))=NaN;
rout(find(rout<a))=NaN;
Boutr = 0;
Binr = 0;
BoutQ = 0 ;
BinQ = 0;

for k=1:2:43
    fun10 = legendre(k,0);
    fun0 = fun10(2);
    fun = legendre(k,cos(theta));
    rfun = squeeze(fun(1,:,:));
    Qfun = squeeze(fun(2,:,:));
    % -----------------------------------------
    Cr = 1/a^k*(rin).^(k-1); % 球内 r 分量
    Bkinr = -1/2*fun0*Cr.*rfun; % 第 k 项
    Binr = Binr+Bkinr; % 总和
    CQ = 1/k/a^k*rin.^(k-1); % 球内 theta 分量
    BkinQ = -1/2*fun0*CQ.*Qfun;
    BinQ = BinQ+BkinQ;
    % -----------------------------------------
    Dr = a^(k+1)./rout.^(k+2); % 球外 r 分量
    Bkoutr = -1/2*fun0*Dr.*rfun; 
    Boutr = Boutr + Bkoutr;
    DQ = 1/(k+1)*a^(k+1)./rout.^(k+2); % 球外 theta 分量
    BkoutQ = 1/2*fun0*DQ.*Qfun;
    BoutQ = BoutQ+BkoutQ;
end

L1 = isnan(Binr); % 将非数元素找出并赋为零
Binr(find(L1==1))=0; BinQ(find(L1==1))=0;
L2 = isnan(Boutr); % 将非数元素找出并赋为零
Boutr(find(L2==1))=0; BoutQ(find(L2==1))=0;
Br = Binr+Boutr; BQ = BinQ+BoutQ; % 合并成一个场 B
By = Br.*cos(theta)-BQ.*sin(theta); % 换成直角坐标以作流线图
Bx = Br.*sin(theta)+BQ.*cos(theta);
vy = 0; vx = [0:0.01:0.35];
[Vx,Vy] = meshgrid(vx,vy);
cla
streamline(X,Y,Bx,By,Vx,Vy)
hold on
streamline(-X,Y,-Bx,By,-Vx,Vy) % 利用对称性画其余象限
streamline(-X,-Y,-Bx,-By,-Vx,-Vy)
streamline(X,-Y,Bx,-By,Vx,-Vy)
title('B 的场线')
figure
subplot(221)
surf(X,Y,Br)
title('Br 分量')
view(-130,20)
subplot(222)
contour(X,Y,Br,10)
title('Br 的等值线')
subplot(223)
surf(X,Y,BQ)
title('B_{\theta}分量')
subplot(224)
contour(X,Y,BQ,10)
title('B_{\theta}的等值线')
