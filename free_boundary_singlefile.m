function free_boundary_singlefile()
% 自由边界模型单文件整合版（含用户输入和图形保持）

% 清除工作区并关闭现有图形
clear;
close all;

% 用户输入参数
a = input('请输入内在增长率a：');
b = input('请输入种内竞争系数b：');
d = input('请输入扩散率d：');
mu = input('请输入自由边界条件系数mu：');
h0 = input('请输入初始种群范围半径h0：');
umax = input('请输入初始种群的最大值umax：');
sigma = input('请输入初始种群分布的宽度参数sigma：');
T = input('请输入模拟的最终时间T：');
H = input('请输入大空间域边界H（需远大于h0）：');
dr = input('请输入空间步长dr：');
dt = input('请输入时间步长dt：');

% 执行主计算函数
[u, h, t, r] = free_boundary_core(a, b, d, mu, h0, umax, sigma, T, H, dr, dt);

% 添加图形保持命令
disp('计算完成，图形窗口已保持打开状态。');
disp('请手动关闭图形窗口以结束程序。');
drawnow;  % 强制刷新图形
waitforbuttonpress;  % 等待用户操作后才退出
end

function [u, h, t, r] = free_boundary_core(a, b, d, mu, h0, umax, sigma, T, H, dr, dt)
% 核心计算函数
r = 0:dr:H;
t = 0:dt:T;
Nr = length(r);
Nt = length(t);

u = zeros(Nt, Nr);
h = zeros(Nt, 1);

% 初始化
h(1) = h0;
for j = 1:Nr
    if r(j) <= h0
        u(1,j) = umax * exp(-(r(j)^2)/(2*sigma^2));
    else
        u(1,j) = 0;
    end
end

[~, idx] = min(abs(r - h0));
h(1) = r(idx);

% 时间步进
for i = 1:Nt-1
    h_current = h(i);
    [~, hidx] = min(abs(r - h_current));
    
    u_new = u(i,:);
    
    % 空间离散
    for j = 2:hidx-1
        N = 2;  % 保持二维扩散
        if r(j) > 0
            laplacian = d*(u(i,j+1)-2*u(i,j)+u(i,j-1))/dr^2 + ...
                        d*(N-1)/r(j)*(u(i,j+1)-u(i,j-1))/(2*dr);
        else
            laplacian = d*2*(u(i,2)-u(i,1))/dr^2;
        end
        reaction = a*u(i,j) - b*u(i,j)^2;
        u_new(j) = u(i,j) + dt*(laplacian + reaction);
    end
    
    % 边界处理
    u_new(1) = u_new(2);
    u_new(hidx:end) = 0;
    
    % 更新自由边界
    if hidx > 1
        ur_at_boundary = (u_new(hidx-1)-u_new(hidx))/dr;
        if hidx > 3
            ur_smoothed = (3*ur_at_boundary + 2*(u_new(hidx-2)-u_new(hidx-1))/dr + ...
                          (u_new(hidx-3)-u_new(hidx-2))/dr)/6;
            ur_at_boundary = ur_smoothed;
        end
        h(i+1) = h(i) - mu*ur_at_boundary*dt;
        h(i+1) = max(h(i), h(i+1));
    else
        h(i+1) = h(i);
    end
    
    [~, new_hidx] = min(abs(r - h(i+1)));
    h(i+1) = r(new_hidx);
    u(i+1,:) = u_new;
    
    if h(i+1) >= H-dr
        u = u(1:i+1,:);
        h = h(1:i+1);
        t = t(1:i+1);
        break;
    end
end

% 可视化
figure(1)
plot(t, h, 'LineWidth', 2)
xlabel('Time t')
ylabel('Spreading radius h(t)')
title('Evolution of spreading front')
grid on

figure(2)
hold on
times_to_plot = unique(round(linspace(1, length(t), 5)));
colors = lines(length(times_to_plot));
for i = 1:length(times_to_plot)
    idx = times_to_plot(i);
    plot(r, u(idx,:), 'Color', colors(i,:), 'LineWidth', 1.5)
end
xlabel('Radius r')
ylabel('Population u(t,r)')
title('Population profiles at different times')
legend_str = arrayfun(@(x)sprintf('t = %.1f', t(x)), times_to_plot, 'UniformOutput', false);
legend(legend_str, 'Location', 'best')
grid on
set(gcf, 'Color', 'w')  % 设置白色背景

% 确保图形窗口保持打开
drawnow
end