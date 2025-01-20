function varargout = exm070701(N,varargin)
% exm070701  Plot N-level Sigmoid function and its 1st derivative.
% 最完整调用格式 [x,y,dy1] = exm070701(N,a,nx,der1,flagplot)
% 输入量
%       N                      电平数N不得小于2
%       a                      指数函数的衰减系数，取值不小于2
%       nx                     Sigmoid函数自变量的采样点数
%       der1                   der1非零时，计算一阶导函数；否则，不算
%       flagplot               flagplot取0时，不画曲线；非零或缺省时，画曲线。
% 输出量
%       x                      Sigmoid函数自变量
%       y                      Sigmoid函数值
%       dy1                    Sigmoid的一阶导数
% 调用格式
% exm070701                    在a=2，nx=100绘出默认的2值Sigmoid函数
% exm070701(N)                 按a=2，nx=100默认值绘出N值Sigmoid函数
% exm070701(N,a)               按输入的a，采用nx=100采样点绘制N值Sigmoid函数
% exm070701(N,a,nx)            按输入的a，nx绘制N值Sigmoid函数
% exm070701(N,a,nx,der1)       按输入的a，nx绘制N值Sigmoid函数及其导函数
% [x,y] = exm070701(N,a,nx)    按输入的a，nx绘制N值Sigmoid函数,并输出其数据
% [x,y,dy1] = exm070701(N,a,nx,der1)
%             按输入参数N，a，nx，绘制Sigmoid函数、其导函数，并输出这两条曲线的数据
% [x,y,dy1] = exm070701(N,a,nx,der1,flagplot)
%             按输入参数N，a，nx计算并输出Sigmoid函数及其导函数；flagplot控制绘图

% Written by Chris Yang, 2022-07-19

if ~any(nargout==[0,2,3])     % 该函数被调用时，若输出量数目不是0，2，3，则给出“出错警告”，并终止程序的执行
    error('输出量数目必须在集合{0 2 3}中')
end
if nargin<4 && nargout==3      % && 为“与”运算快捷方式
    error('想获得导函数数据，输入量数目必须大于等于4')
end
flagplot = 1;

switch nargin
    case 0
        N = 2;a = 2;nx = 100;der1 = 0;
    case 1
        a = 2;nx = 100; der1 = 0;
    case 2
        a = varargin{1}; % 从第一胞元获取a参考值，注意：花括号
        nx = 100;der1 = 0;
    case 3
        a = varargin{1};
        nx = varargin{2};
        der1 = 0;
    case {4,5}
        a = varargin{1};
        nx = varargin{2};
        der1 = varargin{3};
        if nargin == 5
            flagplot = varargin{4};
        end
end

N1 = N-1;
x = linspace(-N,N,nx);
dx = 2*N/nx;
b = -(N-2):2:(N-2);          % N级电平Sigmoid函数的拐点位置
fsum = 0;
for ii=1:N1
    fsum = fsum+2./(1+exp(-a*(x+b(ii))));
end
y = fsum-N1;

if der1 == 0
    h0 = plot(x,y,'-r',x,x,'-g');    %画Sigmoid函数和y=x参照线
    set(h0(1),'LineWidth',2)
    axis([-N1,N1,-N1,N1]),axis equal,grid on
    xlabel('x'),ylabel('y')
    STR = '-level Sigmoid Function with a = ';
    title([num2str(N),STR,num2str(a)])
    legend('Sigmoid','y=x','Location','east')
elseif flagplot ~=0
    dy1 = gradient(y)/dx;           %计算近似导函数
    [ax,h1,h2] = plotyy(x,y,x,dy1); %采用“双纵轴坐标”表现函数及其导数
    grid on
    set(ax(1),'YColor','r')         %设置左纵轴及刻度的颜色
    set(get(ax(1),'Ylabel'),'String','y')
    line(x,x,'Color','g')
    set(h1,'Color','r','LineWidth',2) %设置Sigmoid曲线的颜色及线宽
    set(h2,'LineWidth',2)             %设置导函数曲线线宽
    set(get(ax(2),'Ylabel'),'String','dydx')
    xlabel('x')
    STR = '-level Sigmoid Function and its 1st derivative with a = ';
    title([num2str(N),STR,num2str(a)])
    legend('Sigmoid','y=x','dy1dx','Location','east')
else
    dy1 = gradient(y)/dx;
end

if nargout == 2
    varargout{1} = x;
    varargout{2} = y;
elseif nargout == 3
    varargout{1} = x;
    varargout{2} = y;
    varargout{3} = dy1;
end