function randgui(randfun)
%RANDGUI   Monte Carlo computation of pi.
%  Generate random points in a cube and count the portion that are
%  also in the inscribed sphere.  The ratio of the volume of the
%  sphere to the volume of the cube is pi/6.  
%
%  RANDGUI with no arguments or RANDGUI(＠rand) uses MATLAB's
%  built-in random number generator.  RANDGUI(＠randtx) uses our
%  textbook version of the built-in generator.  RANDGUI(＠randmcg)
%  and RANDGUI(＠randssp) use different Lehmer congruential generators,
%  one with good parameters and one with the parameters used years
%  ago by IBM's "RANDU" function.
%
% See also RAND, RANDTX, RANDMCG, RANDSSP.
% Copyright 6/Nov/2013 Cleve Moler and The MathWorks, Inc.
% Revised by 张志涌（Zhiyong Zhang）,2022
%   参阅《MATLAB数值计算（2021修订版）》（张志涌等译）的第9.2节，有助于理解本文件

nmax = 10000;  % Number of samples
m = 25;        % Samples per plot

if nargin < 1, randfun = @rand; end

shg
clf reset
set(gcf,'doublebuffer','on','pos',[232 100 560 580],'toolbar','figure', ...
   'MenuBar','none','numbertitle','off','name','randgui','Color','w')
[xx,yy,zz]=sphere(50);
hs=surf(xx,yy,zz);
set(hs,'FaceColor','y','FaceAlpha',0.3,...
    'EdgeColor',[0.3,0.3,0.3],'EdgeAlpha',0.5)

ax1=gca;
set(ax1,'pos',[.130 .360 .775 .620],'view',[-39.5 29.4],  ...
   'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1],'box','on','BoxStyle','full', ...
   'GridColor','w','plotboxaspectratiomode','manual','LineWidth',2);

h1 = animatedline(NaN,NaN,NaN,'color','red','linestyle','none', ...
   'marker','.');
h2 = animatedline(NaN,NaN,NaN,'color','c','linestyle','none', ...
   'marker','.');
ax2 = axes('pos',[.130 .110 .775 .200],'xlim',[0 nmax],'ylim',[3 3.3],'Box','on');
h3 = animatedline(NaN,NaN,'color',[0 2/3 0],'Marker','.','MarkerSize',1);
h4 = text(.8*nmax,3.25,'','fontsize',14);
line([0 nmax],[pi pi],'color','black','linestyle',':');
line([0 nmax],[3.3 3.3],'color','black','linestyle','-');
rpt = uicontrol('units','norm','pos',[.02 .01 .10 .05], ...
   'style','push','string','repeat','value',0,'userdata',randfun, ...
   'callback','randgui(get(gcbo,''userdata''))');
stop = uicontrol('units','norm','pos',[.14 .01 .10 .05], ...
   'style','toggle','string','stop','value',0);

n = 0;
s = 0;
while n < nmax && get(stop,'value') == 0
   X = 2*randfun(3,m)-1;
   r = sum(X.^2);
   k = (r <= 1);
   pie = 6*(s+cumsum(k))./(n+1:n+m);
   addpoints(h1,X(1,k),X(2,k),X(3,k));
   addpoints(h2,X(1,~k),X(2,~k),X(3,~k));
   addpoints(h3,n+1:n+m,pie);
   set(h4,'string',sprintf('%7.4f',pie(m)))
   drawnow
   n = n + m;
   s = s + sum(k);
end

set(stop,'style','push','userdata',randfun,'string','close', ...
   'callback','close(gcf)')
