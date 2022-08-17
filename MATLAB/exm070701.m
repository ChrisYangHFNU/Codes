function varargout = exm070701(N,varargin)
% exm070701  Plot N-level Sigmoid function and its 1st derivative.
% ���������ø�ʽ [x,y,dy1] = exm070701(N,a,nx,der1,flagplot)
% ������
%       N                      ��ƽ��N����С��2
%       a                      ָ��������˥��ϵ����ȡֵ��С��2
%       nx                     Sigmoid�����Ա����Ĳ�������
%       der1                   der1����ʱ������һ�׵����������򣬲���
%       flagplot               flagplotȡ0ʱ���������ߣ������ȱʡʱ�������ߡ�
% �����
%       x                      Sigmoid�����Ա���
%       y                      Sigmoid����ֵ
%       dy1                    Sigmoid��һ�׵���
% ���ø�ʽ
% exm070701                    ��a=2��nx=100���Ĭ�ϵ�2ֵSigmoid����
% exm070701(N)                 ��a=2��nx=100Ĭ��ֵ���NֵSigmoid����
% exm070701(N,a)               �������a������nx=100���������NֵSigmoid����
% exm070701(N,a,nx)            �������a��nx����NֵSigmoid����
% exm070701(N,a,nx,der1)       �������a��nx����NֵSigmoid�������䵼����
% [x,y] = exm070701(N,a,nx)    �������a��nx����NֵSigmoid����,�����������
% [x,y,dy1] = exm070701(N,a,nx,der1)
%             ���������N��a��nx������Sigmoid�������䵼��������������������ߵ�����
% [x,y,dy1] = exm070701(N,a,nx,der1,flagplot)
%             ���������N��a��nx���㲢���Sigmoid�������䵼������flagplot���ƻ�ͼ

% Written by Chris Yang, 2022-07-19

if ~any(nargout==[0,2,3])     % �ú���������ʱ�����������Ŀ����0��2��3��������������桱������ֹ�����ִ��
    error('�������Ŀ�����ڼ���{0 2 3}��')
end
if nargin<4 && nargout==3      % && Ϊ���롱�����ݷ�ʽ
    error('���õ��������ݣ���������Ŀ������ڵ���4')
end
flagplot = 1;

switch nargin
    case 0
        N = 2;a = 2;nx = 100;der1 = 0;
    case 1
        a = 2;nx = 100; der1 = 0;
    case 2
        a = varargin{1}; % �ӵ�һ��Ԫ��ȡa�ο�ֵ��ע�⣺������
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
b = -(N-2):2:(N-2);          % N����ƽSigmoid�����Ĺյ�λ��
fsum = 0;
for ii=1:N1
    fsum = fsum+2./(1+exp(-a*(x+b(ii))));
end
y = fsum-N1;

if der1 == 0
    h0 = plot(x,y,'-r',x,x,'-g');    %��Sigmoid������y=x������
    set(h0(1),'LineWidth',2)
    axis([-N1,N1,-N1,N1]),axis equal,grid on
    xlabel('x'),ylabel('y')
    STR = '-level Sigmoid Function with a = ';
    title([num2str(N),STR,num2str(a)])
    legend('Sigmoid','y=x','Location','east')
elseif flagplot ~=0
    dy1 = gradient(y)/dx;           %������Ƶ�����
    [ax,h1,h2] = plotyy(x,y,x,dy1); %���á�˫�������ꡱ���ֺ������䵼��
    grid on
    set(ax(1),'YColor','r')         %���������ἰ�̶ȵ���ɫ
    set(get(ax(1),'Ylabel'),'String','y')
    line(x,x,'Color','g')
    set(h1,'Color','r','LineWidth',2) %����Sigmoid���ߵ���ɫ���߿�
    set(h2,'LineWidth',2)             %���õ����������߿�
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