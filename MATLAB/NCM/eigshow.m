function eigshow(arg)
%EIGSHOW Graphical demonstration of eigenvalues and singular values.
%
%   EIGSHOW presents a graphical experiment showing the effect on the
%   the unit circle of the mapping induced by various 2-by-2 matrices.
%   A pushbutton allows the choice of "eig" mode or "svd" mode.
%
%   In eig mode, the mouse can be used to move the vector x around the
%   unit circle.  The resulting trajectory of A*x is plotted.  The object
%   is to find vectors x so that A*x is parallel to x.  Each such x is an
%   eigenvector of A.  The length of A*x is the corresponding eigenvalue.
%
%   In svd mode, the mouse moves two perpendicular unit vectors, x and y.
%   The resulting A*x and A*y are plotted.  When A*x is perpendicular to
%   A*y, then x and y are right singular vectors, A*x and A*y are
%   multiples of left singular vectors, and the lengths of A*x and A*y
%   are the corresponding singular values.
%
%   The figure title is a menu of selected matrices, including some
%   with fewer than two real eigenvectors.  EIGSHOW(A) inserts A,
%   which must be 2-by-2, in the menu.
%
%   Here are some questions to consider:
%      Which matrices are singular?
%      Which matrices have complex eigenvalues?
%      Which matrices have double eigenvalues?
%      Which matrices have eigenvalues equal to singular values?
%      Which matrices have nondiagonal Jordan canonical forms?

%   Copyright 1984-2007 The MathWorks, Inc.
%   Revised by 张志涌（Zhiyong Zhang）,2022


if nargin == 0
   initialize
elseif arg == 0
   action
elseif arg < 0
   setmode(arg)
else
   initialize(arg);
end

%------------------

function initialize(arg)

if nargin == 0
   arg = 6;
end

if isequal(get(gcf,'tag'),'eigshow')
   h = get(gcf,'userdata');
   mats = h.mats;
else
   set(gcf,'numbertitle','off','menubar','none','Color','w')
   h.svd = 0;
   mats = {
      '[5/4 0; 0 3/4]'
      '[5/4 0; 0 -3/4]'
      '[1 0; 0 1]'
      '[0 1; 1 0]'
      '[0 1; -1 0]'
      '[1 3; 4 2]/4'
      '[1 3; 2 4]/4'
      '[3 1; 4 2]/4'
      '[3 1; -2 4]/4'
      '[2 4; 2 4]/4'
      '[2 4; -1 -2]/4'
      '[6 4; -1 2]/4'
      'randn(2,2)'};
end

if all(size(arg)==1)
   if (arg < length(mats))
      mindex = arg;
      A = eval(mats{mindex});
   else
      A = randn(2,2);
      S = ['[' sprintf('%4.2f %4.2f; %4.2f %4.2f',A.') ']'];
      mindex = length(mats);
      mats = [mats(1:mindex-1); {S}; mats(mindex)];
   end
else
   A = arg;
   if ischar(A)
      S = A;
      A = eval(A);
   else
      S = ['[' sprintf('%4.2f %4.2f; %4.2f %4.2f',A.') ']'];
   end
   if any(size(A) ~= 2)
       error('MATLAB:eigshow:InputSizeIncorrect','%s',...
           getString(message('MATLAB:demos:eigshow:InputSizeIncorrect')))
   end
   mats = [{S};  mats];
   mindex = 1;
end

clf
if h.svd
    t = 'svd / (eig)';
else
    t = 'eig / (svd)';
end
uicontrol(...
   'style','pushbutton', ...
   'units','normalized', ...
   'position',[.86 .60 .12 .06], ...
   'string',t, ...
   'value',h.svd, ...
   'callback','eigshow(-1)');
uicontrol(...
   'style','pushbutton', ...
   'units','normalized', ...
   'position',[.86 .50 .12 .06], ...
   'string','帮助',...
   'callback','helpwin eigshow')
uicontrol(...
   'style','pushbutton', ...
   'units','normalized', ...
   'position',[.86 .40 .12 .06], ...
   'string','关闭',...
   'callback','close(gcf)')
uicontrol(...
   'style','popup', ...
   'units','normalized', ...
   'position',[.38 .92 .35 .075], ...
   'string',mats, ...
   'tag','mats', ...
   'fontname','courier', ...
   'fontweight','bold', ...
   'fontsize',12, ...
   'value',mindex, ...
   'callback','eigshow(get(gco,''value''))');
uicontrol(...
   'style','text', ...
   'units','normalized', ...
   'position',[.28 .925 .10 .06], ...
   'BackgroundColor',[1,1,1],...
   'string','A 阵',...
   'fontsize',11);


s = 1.1*max(1,norm(A));
axis([-s s -s s])
axis square;box on
line(0,0,'Color','k','Marker','o')

xcolor = [0 .6 0];
Axcolor = [0 0 .8];
h.A = A;
h.mats = mats;
h.x = initv([1 0]','x',xcolor);
ax=A(:,1);
% tt=h.text.string;
sax=['Ax ',num2str(norm(ax),'%.3f')];
h.Ax = initv(ax,sax,Axcolor);
if h.svd
   h.y = initv([0 1]','y',xcolor);
   ax=A(:,2);
   sax=['Ay ',num2str(norm(ax),'%.3f')];
   h.Ay = initv(A(:,2),sax,Axcolor);
   D=svd(A);Ds=num2str(D,'%.4f'); 
   text(-2,0,{'A 阵奇异值',Ds},'HorizontalAlignment','center','fontweight','bold')
   xlabel('用鼠标转动"x-y 对"，使 Ax 与 Ay 相互垂直，此Ax-Ay向量对就是奇异向量','fontweight','bold')
%    xlabel(getString(message('MATLAB:demos:eigshow:LabelMakeAxPerpendicularToAy')),'fontweight','bold')
   set(gcf,'name','svdshow')
else
   D=eig(A);Ds=num2str(D,'%.4f'); 
   text(-1.4*s,0,{'A 阵特征值',Ds},'HorizontalAlignment','center','fontweight','bold')
   xlabel('用鼠标转动 x，使 Ax 同向或反向于 x，此 Ax 就是特征向量 ','fontweight','bold')
%   xlabel(getString(message('MATLAB:demos:eigshow:LabelMakeAxParallelToX')),'fontweight','bold')
   set(gcf,'name','eigshow')
end
set(gcf,'tag','eigshow', ...
   'userdata',h, ...
   'windowbuttondownfcn', ...
      'eigshow(0); set(gcf,''windowbuttonmotionfcn'',''eigshow(0)'')', ...
   'windowbuttonupfcn', ...
      'set(gcf,''windowbuttonmotionfcn'','''')')

%------------------

function h = initv(v,t,color)
h.mark =animatedline(v(1),v(2),'marker','.','color',color);%,'erase','none');
if color(3)==0
    h.line = line([0 v(1)],[0 v(2)],'color',color,'LineStyle','--');
else
    h.line = line([0 v(1)],[0 v(2)],'color',color);
end
h.text = text(v(1)/2,v(2)/2,t,'fontsize',12,'color',color);%'erase','xor',

%------------------

function action
h = get(gcf,'userdata');
pt = get(gca,'currentpoint');
x = pt(1,1:2)';
x = x/norm(x);
movev(h.x,x);
A = h.A;
movev(h.Ax,A*x);
if h.svd
   y = [-x(2); x(1)];
   movev(h.y,y);
   movev(h.Ay,A*y);
end

%------------------

function movev(h,v)
addpoints(h.mark,v(1),v(2))
%set(h.mark,'xdata',v(1),'ydata',v(2));
set(h.line,'xdata',[0 v(1)],'ydata',[0 v(2)]);
if h.text.Color(3)==0
    set(h.text,'pos',v/2);
else
    tt=h.text.String(1:3);
    set(h.text,'pos',v/2,'string',[tt,num2str(norm(v),'%.3f')])
%     set(h.text,'pos',v/2,'string',[tt,' 长度 ',nu])
end
%------------------
function setmode(~)
h = get(gcf,'userdata');
h.svd = ~h.svd;
set(gcf,'userdata',h)
initialize(get(findobj(gcf,'tag','mats'),'value'))

