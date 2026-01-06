# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 18:45:06 2025

@author: chris
"""

import cst.interface 
cst = cst.interface.DesignEnvironment()
mws = cst.new_mws()
filename = "Patch_Antenna_Python.cst"
mws.save(filename)
modeler = mws.modeler

# 贴片天线建模基本参数
L = 27.3
W = 27.3
Lg = 80
Wg = 80
h = 1.57
d = 2
Yf = 4.2
Frq = [3.2,5.2]
R_pin = 0.2
offset_pin = 0.7*0.5*W

# 在CST中加入结构参数，方便后续手动在CST文件中进行操作
modeler.add_to_history('StoreParameter','StoreParameter("L",%f)'%L)
modeler.add_to_history('StoreParameter','StoreParameter("W",%f)'%W)
modeler.add_to_history('StoreParameter','StoreParameter("Lg",%f)'%Lg)
modeler.add_to_history('StoreParameter','StoreParameter("Wg",%f)'%Wg)
modeler.add_to_history('StoreParameter','StoreParameter("h",%f)'%h)
modeler.add_to_history('StoreParameter','StoreParameter("d",%f)'%d)
modeler.add_to_history('StoreParameter','StoreParameter("Yf",%f)'%Yf)
modeler.add_to_history('StoreParameter','StoreParameter("R_pin",%f)'%R_pin)
modeler.add_to_history('StoreParameter','StoreParameter("offset_pin",%f)'%offset_pin)

line_break = '\n' #换行符，后续用于VBA代码拼接
tab = ' ' #空格符，后续用于VBA代码拼接

# 设置全局单位
VBA = ['With Units',
       '.Geometry "mm"',
       '.Frequency "GHz"',
       '.Time "ns"',
       'End With']
VBA = line_break.join(VBA)
modeler.add_to_history('define units',VBA)

#设置扫频范围
VBA = 'Solver.FrequencyRange "%f","%f"' %(Frq[0],Frq[1])
modeler.add_to_history('define frequency range',VBA)

#背景材料设置
VBA = ['With Background',
       '.ResetBackground',
       '.Type "Normal"',
       'End With']
VBA = line_break.join(VBA)
modeler.add_to_history('define background',VBA)

# 边界条件设置
VBA = ['With Boundary',
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
VBA = line_break.join(VBA)
modeler.add_to_history('define boundary',VBA)

# 新建所需介质材料：Rogers 5880
er1 = 2.2
VBA = ['With Material',
       '.Reset',
       '.Name "Rogers 5880"',
       '.FrqType "all"',
       '.Type "Normal"',
       '.Epsilon %f' %er1,
       '.Create',
       'End With']
VBA = line_break.join(VBA)
modeler.add_to_history('define material:Rogers 5880',VBA)
er2 = 2.08
VBA = ['With Material',
       '.Reset',
       '.Name "teflon"',
       '.FrqType "all"',
       '.Type "Normal"',
       '.Epsilon %f' %er2,
       '.Create',
       'End With']
VBA = line_break.join(VBA)
modeler.add_to_history('define material:teflon',VBA)

#3D建模
#创建介质基板
Str_Name = 'Substrate'
Str_Component = 'Substrate'
Str_Material = 'Rogers 5880'
VBA = ['With Brick',
       '.Reset',
       '.Name "%s"'%Str_Name,
       '.Component "%s"'%Str_Component,
       '.Material "%s"'%Str_Material,
       '.Xrange "-Wg/2","Wg/2"',
       '.Yrange "-Lg/2","Lg/2"',
       '.Zrange "0","h"',
       '.Create',
       'End With']
VBA = line_break.join(VBA)
modeler.add_to_history('define brick:%s:%s'%(Str_Component,Str_Name,),VBA)
#创建辐射贴片
Str_Name = 'Patch'
Str_Component = 'Patch'
Str_Material = 'PEC'
VBA = ['With Brick',
             '.Reset',
             '.Name "%s"' %Str_Name,
             '.Component "%S"' %Str_Component,
             '.Material "%s"' %Str_Material,
             '.Xrange "-W/2","W/2"',
             '.Yrange "-L/2","L/2"',
             '.Zrange "h","h"',
             '.Create',
             'End With']
VBA = line_break.join(VBA)
modeler.add_to_history('define brick:%s:%s' %(Str_Component,Str_Name,),VBA)
#创建地板