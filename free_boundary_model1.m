function [u, h, t, r] = free_boundary_model1(d, mu, h0, umax, sigma, T, H, dr, dt)
% Numerical solution to the free boundary problem for invasive species spread
% 
% Inputs:
%   d - diffusion rate
%   mu - coefficient in the free boundary condition
%   h0 - initial radius of the population range
%   umax - maximum value of the initial population
%   sigma - parameter for the width of the initial population
%   T - final time
%   H - large spatial domain boundary (>> h0)
%   dr - spatial step size
%   dt - time step size
%
% Outputs:
%   u - solution matrix, u(i,j) is the population at time t(i) and position r(j)
%   h - vector of spreading front positions
%   t - time vector
%   r - spatial position vector

% Create spatial and time grids
r = 0:dr:H;
t = 0:dt:T;
Nr = length(r);
Nt = length(t);

% Initialize solution arrays
u = zeros(Nt, Nr);
h = zeros(Nt, 1);

% Set initial conditions
h(1) = h0;
for j = 1:Nr
    if r(j) <= h0
        u(1,j) = umax * exp(-(r(j)^2)/(2*sigma^2));
    else
        u(1,j) = 0;
    end
end

% Find grid point closest to h0 and adjust
[~, idx] = min(abs(r - h0));
h(1) = r(idx);

% Time stepping
for i = 1:Nt-1
    % Current position of the free boundary
    h_current = h(i);
    
    % Find the index corresponding to h_current
    [~, hidx] = min(abs(r - h_current));
    
    % Compute u at the next time step within the current domain
    u_new = u(i,:);
    
    % Interior points (using Crank-Nicolson scheme)
    for j = 2:hidx-1
        % Handle radial diffusion term with dimension N
        N = 2; % Can change this for different dimensions (1, 2 or 3)
        
        % Using a semi-implicit scheme for stability
        % Diffusion term: duxx + (N-1)/r * dur
        if r(j) > 0
            laplacian = d * (u(i,j+1) - 2*u(i,j) + u(i,j-1))/dr^2 + d*(N-1)/r(j) * (u(i,j+1) - u(i,j-1))/(2*dr);
        else
            % At r=0, use symmetry: dur(0) = 0
            laplacian = d * 2*(u(i,2) - u(i,1))/dr^2; % Factor 2 accounts for the dimension
        end
        
        % Reaction term
        reaction = u(i,j) - u(i,j)^2;
        
        % Update
        u_new(j) = u(i,j) + dt * (laplacian + reaction);
    end
    
    % Apply boundary conditions
    % At r=0: zero flux (symmetry)
    u_new(1) = u_new(2);
    
    % At r=h(t): u=0
    u_new(hidx:end) = 0;
    
    % Update the free boundary position
    % Compute ur at the boundary using one-sided difference
    if hidx > 1
        ur_at_boundary = (u_new(hidx-1) - u_new(hidx))/dr;
        
        % Smooth ur near the boundary to avoid oscillations
        if hidx > 3
            weighted_ur = (3*ur_at_boundary + 2*(u_new(hidx-2) - u_new(hidx-1))/dr + (u_new(hidx-3) - u_new(hidx-2))/dr)/6;
            ur_at_boundary = weighted_ur;
        end
        
        % Update h(t) using the free boundary condition h'(t) = -mu*ur
        h(i+1) = h(i) - mu * ur_at_boundary * dt;
        
        % Ensure h(i+1) doesn't go backward
        h(i+1) = max(h(i), h(i+1));
    else
        h(i+1) = h(i);
    end
    
    % Find the new grid point closest to h(i+1)
    [~, new_hidx] = min(abs(r - h(i+1)));
    h(i+1) = r(new_hidx);
    
    % Update u for the next time step
    u(i+1,:) = u_new;
    
    % Stop if h reaches the domain boundary
    if h(i+1) >= H-dr
        u = u(1:i+1,:);
        h = h(1:i+1);
        t = t(1:i+1);
        break;
    end
end

% Create a figure with two subplots
figure('Position', [100, 100, 1200, 500]);

% Subplot 1: Spreading front evolution
subplot(1, 2, 1);
plot(t, h, 'LineWidth', 2);
xlabel('Time t');
ylabel('Spreading radius h(t)');
title('Evolution of spreading front');
grid on;

% Subplot 2: Population profiles at different times
subplot(1, 2, 2);
times_to_plot = round(linspace(1, length(t), min(5, length(t))));
colors = lines(length(times_to_plot)); % Different colors for each time
for i = 1:length(times_to_plot)
    idx = times_to_plot(i);
    plot(r, u(idx,:), 'LineWidth', 1.5, 'Color', colors(i,:));
    hold on;
end
xlabel('Radius r');
ylabel('Population u(t,r)');
title('Population profiles at different times');

% Create legend with time values
legend_str = cell(length(times_to_plot), 1);
for i = 1:length(times_to_plot)
    legend_str{i} = sprintf('t = %.2f', t(times_to_plot(i)));
end
legend(legend_str, 'Location', 'northeast');
grid on;
hold off;
end