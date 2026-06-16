# -*- coding: utf-8 -*-
"""
Spyder 编辑器

这是一个临时脚本文件。
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 核心材料参数与公式定义
alpha = 0.41
l1, l2 = -0.030, 0.030  # 变换区域虚拟边界 (-30mm 到 30mm)

# 压缩因子 Lambda(x)
def get_lambda(x_virt):
    # 限制在边界内
    x_virt = np.clip(x_virt, l1, l2)
    return (1 - alpha) * ((2 * x_virt - l1 - l2) / (l1 - l2))**2 + alpha

# 导数 d(Lambda)/dx，用于计算材料梯度
def get_d_lambda_dx(x_virt):
    if x_virt < l1 or x_virt > l2:
        return 0.0
    return 2 * (1 - alpha) * (2 * x_virt - l1 - l2) / (l1 - l2)**2 * 2

# 物理空间 x' 到虚拟空间 x 的逆映射 (由于单调递增，用数值查找)
# 为了代码简化和运行效率，直接在虚拟空间构建射线方程
# 射线在虚拟空间走直线，在物理空间走曲线。我们通过解物理空间微分方程来绘制。

def ray_equations(t, state):
    """
    state = [x_phys, y_phys, k_x, k_y]
    这里我们近似在物理空间直接建立微分方程
    """
    x, y, k_x, k_y = state
    
    # 粗略将物理坐标反推回虚拟坐标以获取当前位置的 mu_y 和梯度
    # 因为本题中压缩带主要在中心，我们近似认为在物理空间其梯度形状保持
    # 实际工程更严谨的做法是在虚拟空间解方程，再映射到物理空间绘图。
    # 下面采用【虚拟空间求解 -> 映射回物理空间】的严谨路径：
    return 

# --- 严谨路径：在虚拟空间求解哈密顿方程 (此时介质各向同性)，再通过坐标变换画到物理空间 ---
def hamiltonian_ray_v(t, state):
    x_v, y_v, k_x, k_y = state
    
    # 在虚拟空间中，由于TO的定义，空间是各向同性均匀的（等效于真空）
    # 射线在虚拟空间实际上就是直线传播！
    # 真正好玩的是：由于物理空间的超材料补偿，物理空间射线弯曲，
    # 我们直接计算物理空间的轨迹：
    
    # 物理空间方程：
    pass

# 为了让代码最直观、无第三方高阶库依赖，我们直接计算物理空间射线步进：
def trace_rays_in_physical_space():
    # 建立虚拟空间坐标轴
    x_v_arr = np.linspace(l1, l2, 1000)
    dx_v = x_v_arr[1] - x_v_arr[0]
    
    # 积分算出物理空间坐标
    x_p_arr = np.zeros_like(x_v_arr)
    x_p_arr[0] = x_v_arr[0]
    for i in range(1, len(x_v_arr)):
        x_p_arr[i] = x_p_arr[i-1] + get_lambda(x_v_arr[i-1]) * dx_v
        
    # 物理空间的总长度
    total_phys_len = x_p_arr[-1] - x_p_arr[0]
    
    # 射线参数：从左侧边界 (x_p_arr[0]) 发射
    x0_p = x_p_arr[0]
    # 发射 5 条不同角度的射线 (模拟天线馈源的发散角)
    launch_angles = np.linspace(-25, 25, 5) # 角度范围 -25度 到 +25度
    
    plt.figure(figsize=(10, 6))
    
    # 绘制透镜背景（磁导率高低）
    mu_y_arr = 1.0 / (get_lambda(x_v_arr)**2)
    X_bg, Y_bg = np.meshgrid(x_p_arr * 1000, np.linspace(-40, 40, 100))
    Mu_bg = np.tile(mu_y_arr, (100, 1))
    plt.pcolormesh(X_bg, Y_bg, Mu_bg, cmap='YlGnBu', shading='auto', alpha=0.3)
    cbar = plt.colorbar()
    cbar.set_label('Theoretical $\mu_y$')

    # 射线步进仿真 (物理空间)
    dt = 0.0002 # 时间步长
    steps = 400
    
    for angle in launch_angles:
        # 初始状态
        x_p = x0_p
        y_p = 0.0 # 从中心点发射
        
        # 初始波矢 (考虑折射率)
        rad = np.radians(angle)
        k_x = np.cos(rad)
        k_y = np.sin(rad)
        
        ray_x = [x_p * 1000]
        ray_y = [y_p * 1000]
        
        for _ in range(steps):
            # 找到当前物理位置对应的虚拟位置索引
            idx = np.searchsorted(x_p_arr, x_p)
            idx = np.clip(idx, 0, len(x_v_arr)-1)
            curr_xv = x_v_arr[idx]
            
            # 获取当前位置的材料属性
            lam = get_lambda(curr_xv)
            mu_y = 1.0 / (lam**2)
            
            # 计算导数 d(mu_y)/dx = -2 / lam^3 * d(lam)/dx
            d_lam_dx = get_d_lambda_dx(curr_xv)
            d_muy_dx = -2.0 / (lam**3) * d_lam_dx
            
            # 哈密顿射线步进方程 (物理空间)
            # dx'/dt = 2*k_x / mu_y
            # dy'/dt = 2*k_y
            # dk_x/dt = (k_x^2 / mu_y^2) * d(mu_y)/dx
            
            dx_dt = 2 * k_x / mu_y
            dy_dt = 2 * k_y
            dkx_dt = (k_x**2 / (mu_y**2)) * d_muy_dx
            
            # 更新位置和动量
            x_p += dx_dt * dt
            y_p += dy_dt * dt
            k_x += dkx_dt * dt
            
            # 记录轨迹
            ray_x.append(x_p * 1000)
            ray_y.append(y_p * 1000)
            
            # 如果出了边界就停止追踪
            if x_p > x_p_arr[-1] or np.abs(y_p) > 0.04:
                break
                
        plt.plot(ray_x, ray_y, 'r-', linewidth=2)
        
    plt.title('Ray Tracing Simulation inside Compressed Rotman Lens')
    plt.xlabel('Physical X\' (mm)')
    plt.ylabel('Physical Y\' (mm)')
    plt.xlim([x_p_arr[0]*1000, x_p_arr[-1]*1000])
    plt.ylim([-40, 40])
    plt.grid(True, linestyle=':')
    plt.show()

if __name__ == '__main__':
    trace_rays_in_physical_space()