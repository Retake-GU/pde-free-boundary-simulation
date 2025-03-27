% 脚本名称：run_free_boundary_model.m
% 该脚本允许用户输入参数并运行自由边界模型

% 清除现有变量和图形
clear;
close all;

% 用户输入参数
d = input('请输入扩散率d：');
mu = input('请输入自由边界条件系数mu：');
h0 = input('请输入初始种群范围半径h0：');
umax = input('请输入初始种群的最大值umax：');
sigma = input('请输入初始种群分布的宽度参数sigma：');
T = input('请输入模拟的最终时间T：');
H = input('请输入大空间域边界H（需远大于h0）：');
dr = input('请输入空间步长dr：');
dt = input('请输入时间步长dt：');

% 调用函数进行模拟
[u, h, t, r] = free_boundary_model(d, mu, h0, umax, sigma, T, H, dr, dt);