%ccly_Patch_Antenna
%%MATLAB控制CST的思想基础：CST的建模仿真的命令都是VB命令，而要想使用MATLAB来控制CST就需要将应用之间的通信接口作为一个桥梁架接在MATLAB与CST之间，COM组件就是做这样的事情。
%%由于CST程序中的帮助文件关于VB语言的使用非常详细，所以我们需要做的就是将这座桥梁搭好。注：理论上所有手动操作CST的VB命令能在帮助文件里面找到对应代码
%%将VB命令以字符串的方式从MATLAB传递到CST中运行即可，以下采用的方式可以将关键的建模和设置命令都保存进CST的History List里面，这样就避免了之前采用的MATLABA-CST联合控制方式导致CST文件再次打开或者更新之后变成无用的尴尬。
%%典型的VB命令转化为CST代码如下所示[1]:
%                                               sCommand = '';
% with Brick                                    sCommand = [sCommand 'With Brick'];
%     .Reset                                    sCommand = [sCommand 10 '.Reset'];
%     .Name ("brick1")                          sCommand = [sCommand 10 '.Name ', '"brick1"'];
%     .Component ("component1")                 sCommand = [sCommand 10 '.Component ', '"component1"'];
%     .Material ("PEC")               ==》      sCommand = [sCommand 10 '.Material ', '"PEC"'];
%     .Xrange (0, 2)                            sCommand = [sCommand 10 '.Xrange ', '0',',2'];
%     .Yrange (0, 3)                            sCommand = [sCommand 10 '.Yrange ', '0',',3'];
%     .Zrange (0, "a+3")                        sCommand = [sCommand 10 '.Zrange ', '0', ',"',num2str(a+3),'"');
%     .Create                                   sCommand = [sCommand 10 '.Create'];
% End With                                      sCommand = [sCommand 10 'End With'] ;
%                                               invoke(mws, 'AddToHistory','define brick', sCommand);
%基本思路就是找出mws中需要操作的对象，然后对其所有参数进行赋值或者执行命令，然后释放对象

%%参考程序为CST自带的Macro，点击Home->Macros->Edit Macro->Edit/Move/Delete VBA Macro...->Construc\Demo Examples\Hairpin Filter

%%使用中一定注意单引号双引号的位置，调试中可以先看看MATLAB生成的VBA字符串指令和CST中History
%%List中的VBA命令是否一致（包括所有有用的符号）。

clc;
clear;
close all;

%%CST文件初始化
cst = actxserver('CSTStudio.application');%首先载入CST应用控件
mws = invoke(cst, 'NewMWS');%新建一个MWS项目
app = invoke(mws, 'GetApplicationName');%获取当前应用名称
ver = invoke(mws, 'GetApplicationVersion');%获取当前应用版本号
invoke(mws, 'FileNew');%新建一个CST文件
path=pwd;%获取当前m文件夹路径
filename='\Patch_Antenna.cst';%新建的CST文件名字
fullname=[path filename];
invoke(mws, 'SaveAs', fullname, 'True');%True表示保存到目前为止的结果
invoke(mws, 'DeleteResults');%删除之前的结果。注：在有结果的情况下修改模型会出现弹窗提示是否删除结果，这样运行的程序会停止，需等待手动点击弹窗使之消失
%%CST文件初始化结束

%%贴片天线建模基本参数
a=38.6;%贴片长
b=38;%贴片宽
w=1.46;%馈线宽，100欧姆传输线
l=40;%馈线长
lx=100;%基板长
ly=100;%基板宽
ts=2;%基板厚
tm=0.035;%金属层厚
Frq=[2,2.7];%工作频率，单位：GHz
%在CST中加入结构参数，方便后续手动在CST文件中进行操作
invoke(mws, 'StoreParameter','a',a);
invoke(mws, 'StoreParameter','b',b);
invoke(mws, 'StoreParameter','w',w);
invoke(mws, 'StoreParameter','l',l);
invoke(mws, 'StoreParameter','lx',lx);
invoke(mws, 'StoreParameter','ly',ly);
invoke(mws, 'StoreParameter','ts',ts);
invoke(mws, 'StoreParameter','tm',tm);
%%建模基本参数设置结束

%%全局单位初始化
sCommand = '';
sCommand = [sCommand 'With Units' ];
sCommand = [sCommand 10 '.Geometry "mm"'];%10在这里是换行的作用
sCommand = [sCommand 10 '.Frequency "ghz"' ];
sCommand = [sCommand 10 '.Time "ns"'];
sCommand = [sCommand 10 'End With'] ;
invoke(mws, 'AddToHistory','define units', sCommand);
%%全局单位初始化结束

%%工作频率设置
sCommand = '';
sCommand = [sCommand 'Solver.FrequencyRange '  num2str(Frq(1)) ',' num2str(Frq(2)) ];
invoke(mws, 'AddToHistory','define frequency range', sCommand);
%%工作频率设置结束

%%背景材料设置
sCommand = '';
sCommand = [sCommand 'With Background' ];
sCommand = [sCommand 10 '.ResetBackground'];
sCommand = [sCommand 10 '.Type "Normal"' ];
sCommand = [sCommand 10 'End With'] ;
invoke(mws, 'AddToHistory','define background', sCommand);
%%背景材料设置结束

%%边界条件设置。
sCommand = '';
sCommand = [sCommand 'With Boundary' ];
sCommand = [sCommand 10 '.Xmin "expanded open"'];
sCommand = [sCommand 10 '.Xmax "expanded open"'];
sCommand = [sCommand 10 '.Ymin "expanded open"'];
sCommand = [sCommand 10 '.Ymax "expanded open"'];
sCommand = [sCommand 10 '.Zmin "expanded open"'];
sCommand = [sCommand 10 '.Zmax "expanded open"'];
sCommand = [sCommand 10 '.Xsymmetry "none"' ];
sCommand = [sCommand 10 '.Ysymmetry "none"' ];
sCommand = [sCommand 10 '.Zsymmetry "none"' ];
sCommand = [sCommand 10 'End With'] ;
invoke(mws, 'AddToHistory','define boundary', sCommand);
%%边界条件设置结束

%%新建所需介质材料
er1 = 2.65;
sCommand = '';
sCommand = [sCommand 'With Material' ];
sCommand = [sCommand 10 '.Reset'];
sCommand = [sCommand 10 '.Name "material1"' ];
sCommand = [sCommand 10 '.FrqType "all"' ];
sCommand = [sCommand 10 '.Type "Normal"' ];
sCommand = [sCommand 10 '.Epsilon ' num2str(er1) ];
sCommand = [sCommand 10 '.Create'];
sCommand = [sCommand 10 'End With'] ;
invoke(mws, 'AddToHistory','define material: material265', sCommand);
%%新建所需介质材料结束

%%使Bounding Box显示,这段代码不用保存进历史树
plot = invoke(mws, 'Plot');
invoke(plot, 'DrawBox', 'True');
%%使Bounding Box显示结束

%%建模开始
%%调用Brick对象开始
Str_Name='patch';
Str_Component='Patch';
Str_Material='PEC';
%以下这一串可以写成函数
sCommand = '';
sCommand = [sCommand 'With Brick'];
sCommand = [sCommand 10 '.Reset'];
sCommand = [sCommand 10 '.Name "',Str_Name, '"'];
sCommand = [sCommand 10 '.Component "', Str_Component, '"'];
sCommand = [sCommand 10 '.Material "', Str_Material, '"'];
sCommand = [sCommand 10 '.Xrange ', '"-a/2", "a/2"'];%这里的变两名作为字符串对应CST里面参数列表的变量了,这里的变量名假如不加双引号，那么在CST里面手动修改参数列表里的数据后不会出现提示让你按F7更新，而是需要自己打开历史列表点击更新
sCommand = [sCommand 10 '.Yrange ', '"-b/2","b/2"'];
sCommand = [sCommand 10 '.Zrange ', '"0","tm"'];
sCommand = [sCommand 10 '.Create'];
sCommand = [sCommand 10 'End With'] ;
invoke(mws, 'AddToHistory',['define brick:',Str_Component,':',Str_Name], sCommand);
%以上这一串可以写成函数

invoke(plot, 'ZoomToStructure');%缩放到适合大小，就和在CST里面按空格是一个效果，没事按按空格看下整体结构很重要

Str_Name='line1';
Str_Component='Feed';
Str_Material='PEC';
sCommand = '';
sCommand = [sCommand 'With Brick'];
sCommand = [sCommand 10 '.Reset'];
sCommand = [sCommand 10 '.Name "',Str_Name, '"'];
sCommand = [sCommand 10 '.Component "', Str_Component, '"'];
sCommand = [sCommand 10 '.Material "', Str_Material, '"'];
sCommand = [sCommand 10 '.Xrange ', '"-lx/2","-a/2"'];%这里的变两名作为字符串对应CST里面参数列表的变量了
sCommand = [sCommand 10 '.Yrange ', '"-w/2","w/2"'];
sCommand = [sCommand 10 '.Zrange ', '"0","tm"'];
sCommand = [sCommand 10 '.Create'];
sCommand = [sCommand 10 'End With'] ;
invoke(mws, 'AddToHistory',['define brick:',Str_Component,':',Str_Name], sCommand);

Str_Name='line2';
Str_Component='Feed';
Str_Material='PEC';
sCommand = '';
sCommand = [sCommand 'With Brick'];
sCommand = [sCommand 10 '.Reset'];
sCommand = [sCommand 10 '.Name "',Str_Name, '"'];
sCommand = [sCommand 10 '.Component "', Str_Component, '"'];
sCommand = [sCommand 10 '.Material "', Str_Material, '"'];
sCommand = [sCommand 10 '.Xrange ', '"a/2","lx/2"'];%这里的变两名作为字符串对应CST里面参数列表的变量了
sCommand = [sCommand 10 '.Yrange ', '"-w/2","w/2"'];
sCommand = [sCommand 10 '.Zrange ', '"0","tm"'];
sCommand = [sCommand 10 '.Create'];
sCommand = [sCommand 10 'End With'] ;
invoke(mws, 'AddToHistory',['define brick:',Str_Component,':',Str_Name], sCommand);

Str_Name='bottom';
Str_Component='Bottom';
Str_Material='PEC';
sCommand = '';
sCommand = [sCommand 'With Brick'];
sCommand = [sCommand 10 '.Reset'];
sCommand = [sCommand 10 '.Name "',Str_Name, '"'];
sCommand = [sCommand 10 '.Component "', Str_Component, '"'];
sCommand = [sCommand 10 '.Material "', Str_Material, '"'];
sCommand = [sCommand 10 '.Xrange ', '"-lx/2","lx/2"'];%这里的变两名作为字符串对应CST里面参数列表的变量了
sCommand = [sCommand 10 '.Yrange ', '"-ly/2","ly/2"'];
sCommand = [sCommand 10 '.Zrange ', '"-ts-tm","-ts"'];
sCommand = [sCommand 10 '.Create'];
sCommand = [sCommand 10 'End With'] ;
invoke(mws, 'AddToHistory',['define brick:',Str_Component,':',Str_Name], sCommand);

invoke(plot, 'ZoomToStructure');

Str_Name='sub';
Str_Component='Sub';
Str_Material='material1';
sCommand = '';
sCommand = [sCommand 'With Brick'];
sCommand = [sCommand 10 '.Reset'];
sCommand = [sCommand 10 '.Name "',Str_Name, '"'];
sCommand = [sCommand 10 '.Component "', Str_Component, '"'];
sCommand = [sCommand 10 '.Material "', Str_Material, '"'];
sCommand = [sCommand 10 '.Xrange ', '"-lx/2","lx/2"'];%这里的变两名作为字符串对应CST里面参数列表的变量了
sCommand = [sCommand 10 '.Yrange ', '"-ly/2","ly/2"'];
sCommand = [sCommand 10 '.Zrange ', '"-ts","0"'];
sCommand = [sCommand 10 '.Create'];
sCommand = [sCommand 10 'End With'] ;
invoke(mws, 'AddToHistory',['define brick:',Str_Component,':',Str_Name], sCommand);

%%建模结束

invoke(plot, 'ZoomToStructure');

%%端口设置，采用的方法和在CST里面选中一个面然后设置端口是一样的操作，这里完全复现
%端口1
sCommand = '';
sCommand = [sCommand 'Pick.PickFaceFromId "','Feed:line1"',',4'];
invoke(mws, 'AddToHistory','pick face', sCommand);
sCommand = '';
sCommand = [sCommand 10 'With Port'];
sCommand = [sCommand 10 '.Reset'];
sCommand = [sCommand 10 '.PortNumber ', '1'];
sCommand = [sCommand 10 '.Label ', '""'];
sCommand = [sCommand 10 '.NumberOfModes ', '1'];
sCommand = [sCommand 10 '.AdjustPolarization ', '"False"'];
sCommand = [sCommand 10 '.PolarizationAngle ', '0.0'];
sCommand = [sCommand 10 '.ReferencePlaneDistance ', '0'];
sCommand = [sCommand 10 '.TextSize ', '50'];
sCommand = [sCommand 10 '.TextMaxLimit ', '0'];
sCommand = [sCommand 10 '.Coordinates ', '"Picks"'];
sCommand = [sCommand 10 '.Orientation ', '"positive"'];
sCommand = [sCommand 10 '.PortOnBound ', '"False"'];
sCommand = [sCommand 10 '.ClipPickedPortToBound ', '"False"'];
sCommand = [sCommand 10 '.Xrange ', '"-lx/2","-lx/2"'];
sCommand = [sCommand 10 '.Yrange ', '"-w/2","w/2"'];
sCommand = [sCommand 10 '.Zrange ', '"0","tm"'];
sCommand = [sCommand 10 '.XrangeAdd ', '"0.0","0.0"'];
sCommand = [sCommand 10 '.YrangeAdd ', '"3*ts","3*ts"'];
sCommand = [sCommand 10 '.ZrangeAdd ', '"ts","3*ts"'];
sCommand = [sCommand 10 '.SingleEnded ', '"False"'];
sCommand = [sCommand 10 '.Create'];
sCommand = [sCommand 10 'End With'];
invoke(mws, 'AddToHistory','define port1', sCommand);
%端口2
sCommand = '';
sCommand = [sCommand 'Pick.PickFaceFromId "','Feed:line2"',',6'];
invoke(mws, 'AddToHistory','pick face', sCommand);
sCommand = '';
sCommand = [sCommand 10 'With Port'];
sCommand = [sCommand 10 '.Reset'];
sCommand = [sCommand 10 '.PortNumber ', '2'];
sCommand = [sCommand 10 '.Label ', '""'];
sCommand = [sCommand 10 '.NumberOfModes ', '1'];
sCommand = [sCommand 10 '.AdjustPolarization ', '"False"'];
sCommand = [sCommand 10 '.PolarizationAngle ', '0.0'];
sCommand = [sCommand 10 '.ReferencePlaneDistance ', '0'];
sCommand = [sCommand 10 '.TextSize ', '50'];
sCommand = [sCommand 10 '.TextMaxLimit ', '0'];
sCommand = [sCommand 10 '.Coordinates ', '"Picks"'];
sCommand = [sCommand 10 '.Orientation ', '"positive"'];
sCommand = [sCommand 10 '.PortOnBound ', '"False"'];
sCommand = [sCommand 10 '.ClipPickedPortToBound ', '"False"'];
sCommand = [sCommand 10 '.Xrange ', '"lx/2","lx/2"'];
sCommand = [sCommand 10 '.Yrange ', '"-w/2","w/2"'];
sCommand = [sCommand 10 '.Zrange ', '"0","tm"'];
sCommand = [sCommand 10 '.XrangeAdd ', '"0.0","0.0"'];
sCommand = [sCommand 10 '.YrangeAdd ', '"3*ts","3*ts"'];
sCommand = [sCommand 10 '.ZrangeAdd ', '"ts","3*ts"'];
sCommand = [sCommand 10 '.SingleEnded ', '"False"'];
sCommand = [sCommand 10 '.Create'];
sCommand = [sCommand 10 'End With'];
invoke(mws, 'AddToHistory','define port1', sCommand);
%%端口设置结束

%%设置远场方向图的Monitor
monitor = invoke(mws, 'Monitor');
farfield_monitor = 2.1:0.05:2.6;
for i = 1:length(farfield_monitor)
    Str_name = ['Farfield (f=',num2str(farfield_monitor(i)),')'];
    sCommand = '';
    sCommand = [sCommand 'With Monitor'];
    sCommand = [sCommand 10 '.Reset'];
    sCommand = [sCommand 10 '.Name "', Str_name,'"'];
    sCommand = [sCommand 10 '.Dimension "Volume"'];
    sCommand = [sCommand 10 '.Domain "Frequency"'];
    sCommand = [sCommand 10 '.FieldType "Farfield"'];
    sCommand = [sCommand 10 '.Frequency ', num2str(farfield_monitor(i))];%这里不是引用列入CST参数列表的变量，因此不用在变量名两侧加入双引号
    sCommand = [sCommand 10 '.Create'];
    sCommand = [sCommand 10 'End With'];
    invoke(mws, 'AddToHistory',['define farfield monitor:farfiled (f=', num2str(farfield_monitor(i)),')'],sCommand);
end
%%设置远场方向图的Monitor结束

%%仿真开始
solver = invoke(mws, 'Solver');
invoke(solver, 'Start');
%%仿真结束

invoke(mws, 'Save');%保存
%invoke(mws, 'Quit');%退出
