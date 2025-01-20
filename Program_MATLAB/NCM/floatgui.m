function floatgui(~)
%FLOATGUI  Show structure of floating point numbers.
%  The set of positive model floating point numbers is determined
%  by three parameters: t, emin, and emax.  It is the set of rational
%  numbers of the form x = (1+f)*2^e where f = (integer)/2^t,
%  0 <= f < 1, e = integer, and emin <= e <= emax.
%  波浪符 ~ 函数输入量,用来接收本程序形成界面上控件的回调输入量
%
%  IEEE 754 double precision has t = 52, emin = -1022, emax = 1023.

%   Copyright 2013 Cleve Moler and The MathWorks, Inc.
%   Revised by 张志涌（Zhiyong Zhang）,2022

% Initialize parameters

if nargin == 0
   t = 3;
   emin = -4;
   emax = 2;
   logscale = 0;
   Fpos=[50 300 900 250];
else
   t = round(get(findobj('tag','t'),'value'));
   emin = round(get(findobj('tag','emin'),'value'));
   emax = round(get(findobj('tag','emax'),'value'));
   logscale = get(findobj('style','check'),'value');
   Fpos=get(gcf,'pos');
end

% Position figure window

shg
clf reset
set(gcf,'pos',Fpos,'name','浮点体系演示模型', ...
   'resize','on','defaultuicontrolunits','normalized',...
   'numbertitle','off','menubar','none')%

% Generate and plot floating point numbers

f = (0:2^t-1)/2^t;
F = [];
for e = emin:emax
   F = [F (1+f)*2^e];
end
for x = F
   text(x,0,'|','fontunits','normalized','fontsize',0.3,...
       'hori','center')
end

% Set axes

set(gca,'pos',[.05 .35 .9 .2],'fontunits','normalized',...
    'fontsize',0.22,'linewidth',2,'box','on')
if logscale
   set(gca,'xscale','log')
   xmin = 1/2^(-emin+.5);
   xmax = 2^(emax+1.5);
   XM=log(xmax/4);
else
   set(gca,'xscale','linear')
   xmin = 0;
   xmax = 2^(emax+1);
   XM=xmax/2;
end
axis([xmin xmax -1 1])

% Set tick marks

fmin = min(F);
fmax = max(F);
xtick = 1;
xticklab = {'1'};
if fmin < 1
   xtick = [1/2 xtick];
   xticklab = [{'1/2'}, xticklab];
end
if logscale && (fmin < 1/4)
   xtick = [1/4 xtick];
   xticklab = [{'1/4'}, xticklab];
end
if fmin < 1/2
   xtick = [fmin xtick];
   xticklab = [{['1/' int2str(1/fmin)]}, xticklab];
end
for kk=1:emax
    if 2^kk<fmax
        xtick=[xtick,2^kk];
        xticklab = [xticklab,{int2str(2^kk)}];
    end
end
    
    
    
% if 2 < fmax
%    xtick = [xtick 2];
%    xticklab = [xticklab,{'2'}];
% end
% if 4 < fmax
%    xtick = [xtick 4];
%    xticklab = [xticklab,{'4'}];
% end
if max(xtick) < fmax
   xtick = [xtick fmax];
   if fmax == round(fmax)
      fmaxlab = int2str(fmax);
   else
      over = 2^(emax+1);
      fmaxlab = [int2str(over) '-1/' int2str(1/(over-fmax))];
   end
   xticklab = [xticklab,{fmaxlab}];
end
set(gca,'xtick',xtick,'xticklabel',xticklab,'xminortick','off','ytick',[])

text(XM,4.9,'浮点正数结构模型：x=(1+f)2^e ；f=0:2^{-t}:(1-2^{-t})',...
    'fontsize',20,'hori','center','fontweight','bold',...
    'fontunits','normalized','fontsize',0.35)
% Create uicontrols

uicontrol('style','slider','tag','emin','value',emin, ...
   'min',-8,'max',0,...
   'pos',[0.15 0.12 0.13 0.07],'sliderstep',[1/8 1/8], ...
   'callback','floatgui(1)');
uicontrol('style','slider','tag','t','value',t, ...
   'min',0,'max',8,...
   'pos',[0.435 0.12 0.13 0.07],'sliderstep',[1/8 1/8], ...
   'callback','floatgui(1)');
uicontrol('style','slider','tag','emax','value',emax, ...
   'min',0,'max',8,...
   'pos',[0.72 0.12 0.13 0.07],'sliderstep',[1/8 1/8], ...
   'callback','floatgui(1)');%[0.72 0.26 0.13 0.07]
uicontrol('style','text','string',['emin = ' int2str(emin)], ...
   'pos',[0.15 0.2 0.13 0.07],'fontunits','normalized','fontsize',0.7)
uicontrol('style','text','string',['t = ' int2str(t)], ...
   'pos',[0.435 0.2 0.13 0.07],'fontunits','normalized','fontsize',0.7)
uicontrol('style','text','string',['emax = ' int2str(emax)], ...
   'pos',[0.72 0.2 0.13 0.07],'fontunits','normalized','fontsize',0.7)%[0.72 0.35 0.13 0.07]
uicontrol('style','check','string','log scale','value',logscale, ...
   'pos',[0.435 0.03 0.13 0.07],'fontunits','normalized','fontsize',0.7, ...
   'callback','floatgui(1)');%[0.435 0.15 0.13 0.07]
uicontrol('style','push','pos',[0.88 0.03 0.07 0.07], ...
   'fontunits','normalized','fontsize',0.7, ...
   'string','close','callback','close(gcf)')%[0.88 0.1 0.07 0.07]

% eps

if fmax > 1
   eps = 2^(-t);
%    text(1,0,'|','color','r','fontunits','normalized','fontsize',0.6,'hor','center')
   text(1,0.7,'|','color','r','fontunits','normalized','fontsize',0.8,'hor','center')
%    text(1+eps,0,'|','color','r','fontunits','normalized','fontsize',0.4,'hor','center')
   text(1+eps,0.7,'|','color','r','fontunits','normalized','fontsize',0.8,'hor','center')
   xx=2^emax;Dx=eps*2^emax;
   text(xx,0.7,'|','color','m','fontunits','normalized','fontsize',0.8,'hor','center')
   text(xx+Dx,0.7,'|','color','m','fontunits','normalized','fontsize',0.8,'hor','center')
   if eps < 1
      text(2^emin,2,{'realmin','= 2^{emin}',['= ',num2str(2^emin,'%.4f')]},'hori','left',...
      'fontunits','normalized','fontsize',0.3,'rotation',45)%,'fontweight','bold'
      text(1,2,{'eps','= 1/2^{- t}',[' = ',num2str(eps,'%.4f')]},'color','r',...
      'fontunits','normalized','fontsize',0.3,'hori','left',...
      'rotation',45)%,'fontweight','bold'
      text(xx+Dx/2,2.6,{'eps*2^{emax}',['= ',num2str(eps*2^emax,'%.4f')]},'color','m',...
      'fontunits','normalized','fontsize',0.3,'rotation',45,'hori','center')%,'fontweight','bold'
      text((2-1.5*eps)*2^emax,2,{'realmax','=(2-eps)*2^{emax}',['= ',num2str((2-eps)*2^emax,'%.4f')]},...
      'fontunits','normalized','fontsize',0.3,'rotation',45,'hori','left')%,'fontweight','bold'
   else
      text(1.0,1.5,'eps = 1', 'fontunits','normalized','fontsize',0.3,'fontweight','bold')
   end
end

% Number of numbers

% Exercise:
% How many "floating point" numbers are in the set?
% Complete this statement.
% text(.9*xmax,2,num2str(???)
