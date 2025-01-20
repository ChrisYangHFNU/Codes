#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@ccly
"""
import os
#import numpy as np
import cst.interface
path = os.getcwd()#获取当前py文件所在文件夹
filename = 'Patch_Antenna.cst'
fullname = os.path.join(path,filename)
print(fullname)
cst = cst.interface.DesignEnvironment()
mws = cst.new_mws()
mws.save(fullname)
modeler = mws.modeler
#贴片天线建模基本参数
a = 38.6 #贴片长
b = 38 #贴片宽
w = 1.46 #馈线宽，100欧姆传输线
l = 40 #馈线长
lx = 100 #基板长
ly = 100 #基板宽
ts = 2 #基板厚
tm = 0.035 #金属层厚
Frq = [2,2.7] #工作频率，单位：GHz

#在CST中加入结构参数，方便后续手动在CST文件中进行操作
modeler.add_to_history('StoreParameter','MakeSureParameterExists("a", "%f")' % a)
modeler.add_to_history('StoreParameter','MakeSureParameterExists("b", "%f")' % b)
modeler.add_to_history('StoreParameter','MakeSureParameterExists("w", "%f")' % w)
modeler.add_to_history('StoreParameter','MakeSureParameterExists("l", "%f")' % l)
modeler.add_to_history('StoreParameter','MakeSureParameterExists("lx", "%f")' % lx)
modeler.add_to_history('StoreParameter','MakeSureParameterExists("ly", "%f")' % ly)
modeler.add_to_history('StoreParameter','MakeSureParameterExists("ts", "%f")' % ts)
modeler.add_to_history('StoreParameter','MakeSureParameterExists("tm", "%f")' % tm)
#建模基本参数设置结束
line_break = '\n'#换行符，后面用于VBA代买的拼接用
#全局单位初始化
sCommand = ['With Units',
            '.Geometry "mm"',
            '.Frequency "ghz"',
            '.Time "ns"',
            'End With']
sCommand = line_break.join(sCommand)
modeler.add_to_history('define units', sCommand)
#全局单位初始化结束

#工作频率设置
sCommand = 'Solver.FrequencyRange "%f","%f"'  % (Frq[0],Frq[1])
modeler.add_to_history('define frequency range', sCommand)
#工作频率设置结束

#背景材料设置
sCommand = ['With Background',
            '.ResetBackground',
            '.Type "Normal"',
            'End With']
sCommand = line_break.join(sCommand)
modeler.add_to_history('define background', sCommand)
#背景材料设置结束

#边界条件设置。
sCommand = ['With Boundary',
            '.Xmin "expanded open"',
            '.Xmax "expanded open"',
            '.Ymin "expanded open"',
            '.Ymax "expanded open"',
            '.Zmin "expanded open"',
            '.Zmax "expanded open"',
            '.Xsymmetry "none"',
            '.Ysymmetry "none"',
            '.Zsymmetry "none"',
            'End With']
sCommand = line_break.join(sCommand)
modeler.add_to_history('define boundary', sCommand)
#边界条件设置结束

#新建所需介质材料
er1 = 2.65
sCommand = ['With Material',
            '.Reset',
            '.Name "material1"',
            '.FrqType "all"',
            '.Type "Normal"',
            '.Epsilon %f' %er1,
            '.Create',
            'End With']
sCommand = line_break.join(sCommand)
modeler.add_to_history('define material: material265', sCommand)
#新建所需介质材料结束

#使Bounding Box显示
sCommand = 'Plot.DrawBox "True"'
modeler.add_to_history('switch bounding box', sCommand)
#使Bounding Box显示结束

#建模开始
#调用Brick对象开始
Str_Name='patch'
Str_Component='Patch'
Str_Material='PEC'
#以下这一串可以写成函数
sCommand = ['With Brick',
            '.Reset',
            '.Name "%s"' % Str_Name,
            '.Component "%s"' % Str_Component,
            '.Material "%s"' % Str_Material,
            '.Xrange "-a/2","a/2"',
            '.Yrange "-b/2","b/2"',
            '.Zrange "0","tm"',
            '.Create',
            'End With']
sCommand = line_break.join(sCommand)
modeler.add_to_history('define brick:%s:%s' % (Str_Component,Str_Name,), sCommand)
#以上这一串可以写成函数

sCommand = 'Plot.ZoomToStructure'
modeler.add_to_history('ZoomToStructure', sCommand)#缩放到适合大小，就和在CST里面按空格是一个效果

Str_Name='line1'
Str_Component='Feed'
Str_Material='PEC'
sCommand = ['With Brick',
            '.Reset',
            '.Name "%s"' % Str_Name,
            '.Component "%s"' % Str_Component,
            '.Material "%s"' % Str_Material,
            '.Xrange "-lx/2","-a/2"',
            '.Yrange "-w/2","w/2"',
            '.Zrange "0","tm"',
            '.Create',
            'End With']
sCommand = line_break.join(sCommand)
modeler.add_to_history('define brick:%s:%s' % (Str_Component,Str_Name,), sCommand)

Str_Name='line2'
Str_Component='Feed'
Str_Material='PEC'
sCommand = ['With Brick',
            '.Reset',
            '.Name "%s"' % Str_Name,
            '.Component "%s"' % Str_Component,
            '.Material "%s"' % Str_Material,
            '.Xrange "a/2","lx/2"',
            '.Yrange "-w/2","w/2"',
            '.Zrange "0","tm"',
            '.Create',
            'End With']
sCommand = line_break.join(sCommand)
modeler.add_to_history('define brick:%s:%s' % (Str_Component,Str_Name,), sCommand)

Str_Name='bottom'
Str_Component='Bottom'
Str_Material='PEC'
sCommand = ['With Brick',
            '.Reset',
            '.Name "%s"' % Str_Name,
            '.Component "%s"' % Str_Component,
            '.Material "%s"' % Str_Material,
            '.Xrange "-lx/2","lx/2"',
            '.Yrange "-ly/2","ly/2"',
            '.Zrange "-ts-tm","-ts"',
            '.Create',
            'End With']
sCommand = line_break.join(sCommand)
modeler.add_to_history('define brick:%s:%s' % (Str_Component,Str_Name,), sCommand)

Str_Name='sub'
Str_Component='Sub'
Str_Material='material1'
sCommand = ['With Brick',
            '.Reset',
            '.Name "%s"' % Str_Name,
            '.Component "%s"' % Str_Component,
            '.Material "%s"' % Str_Material,
            '.Xrange "-lx/2","lx/2"',
            '.Yrange "-ly/2","ly/2"',
            '.Zrange "-ts","0"',
            '.Create',
            'End With']
sCommand = line_break.join(sCommand)
modeler.add_to_history('define brick:%s:%s' % (Str_Component,Str_Name,), sCommand)
#建模结束

sCommand = 'Plot.ZoomToStructure'
modeler.add_to_history('ZoomToStructure', sCommand)#缩放到适合大小，就和在CST里面按空格是一个效果

#端口设置，采用的方法和在CST里面选中一个面然后设置端口是一样的操作，这里完全复现
#端口1
sCommand = 'Pick.PickFaceFromId "Feed:line1",4'
modeler.add_to_history('pick face', sCommand)
sCommand = ['With Port',
            '.Reset',
            '.PortNumber 1',
            '.Label  ""',
            '.NumberOfModes 1',
            '.AdjustPolarization "False"',
            '.PolarizationAngle 0.0',
            '.ReferencePlaneDistance 0',
            '.TextSize 50',
            '.TextMaxLimit 0',
            '.Coordinates "Picks"',
            '.Orientation "positive"',
            '.PortOnBound "False"',
            '.ClipPickedPortToBound "False"',
            '.Xrange "-lx/2","-lx/2"',
            '.Yrange "-w/2","w/2"',
            '.Zrange "0","tm"',
            '.XrangeAdd "0.0","0.0"',
            '.YrangeAdd "3*ts","3*ts"',
            '.ZrangeAdd "ts","3*ts"',
            '.SingleEnded "False"',
            '.Create',
            'End With']
sCommand = line_break.join(sCommand)     
modeler.add_to_history('define port1', sCommand)

#端口2
sCommand = 'Pick.PickFaceFromId "Feed:line2",6'
modeler.add_to_history('pick face', sCommand)
sCommand = ['With Port',
            '.Reset',
            '.PortNumber 2',
            '.Label  ""',
            '.NumberOfModes 1',
            '.AdjustPolarization "False"',
            '.PolarizationAngle 0.0',
            '.ReferencePlaneDistance 0',
            '.TextSize 50',
            '.TextMaxLimit 0',
            '.Coordinates "Picks"',
            '.Orientation "positive"',
            '.PortOnBound "False"',
            '.ClipPickedPortToBound "False"',
            '.Xrange "lx/2","lx/2"',
            '.Yrange "-w/2","w/2"',
            '.Zrange "0","tm"',
            '.XrangeAdd "0.0","0.0"',
            '.YrangeAdd "3*ts","3*ts"',
            '.ZrangeAdd "ts","3*ts"',
            '.SingleEnded "False"',
            '.Create',
            'End With']
sCommand = line_break.join(sCommand)     
modeler.add_to_history('define port2', sCommand)
#端口设置结束

#设置远场方向图的Monitor
sCommand = ['With Monitor',
            '.Reset',
            '.Domain "Frequency"',
            '.FieldType "Farfield"',
            '.ExportFarfieldSource "False"',
            '.UseSubvolume "False"',
            '.Coordinates "Picks"',
            '.SetSubvolume "50", "50", "-0.73", "0.73", "0", "0.035"',
            '.SetSubvolumeOffset "10", "10", "10", "10", "10", "10" ',
            '.SetSubvolumeInflateWithOffset "False" ',
            '.SetSubvolumeOffsetType "FractionOfWavelength" ',
            '.EnableNearfieldCalculation "True" ',
            '.CreateUsingLinearStep "%f", "%f", "%f"' % (Frq[0],Frq[1],0.05),
            'End With']
sCommand = line_break.join(sCommand) 
modeler.add_to_history('define farfield monitor (using linear step)',sCommand)
#设置远场方向图的Monitor结束

#仿真开始
modeler. run_solver()
#仿真结束

mws.save(fullname)#保存