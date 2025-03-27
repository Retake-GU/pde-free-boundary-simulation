function [u, v, h, t, r] = free_boundary_model_uv(h0, umax, sigma, T, H, dr, dt, d, mu)
% Simulate the following coupled system:
%   u_t = u_{rr} + u*(1 - u - v),       0 < r < h(t)
%   v_t = d*v_{rr} + r*v*(1 - v - u),     0 < r < H (with H sufficiently large)
%   u_r(t,0) = v_r(t,0) = 0,              (symmetry boundary conditions)
%   h'(t)= -mu*u_r(t,h(t)),               (free boundary condition)
%   u(t,r) = 0,   for r >= h(t)
%
% Initial conditions:
%   u(0, r) = umax * exp(-r^2/(2*sigma^2)) for r <= h0, 0 elsewhere;
%   v(0, r) = v0(r) (here simply set to a constant, e.g. 0.1; can be modified as needed)
%
% Input parameters:
%   h0    - initial radius of the population range
%   umax  - maximum initial population density (used for the u initial distribution)
%   sigma - parameter controlling the width of the initial u distribution
%   T     - final simulation time
%   H     - spatial domain boundary (H >> h0)
%   dr    - spatial step size
%   dt    - time step size
%   d     - diffusion coefficient for v
%   mu    - parameter in the free boundary condition
%
% Output:
%   u - numerical solution matrix for u; u(i,j) is the value at time t(i) and position r(j)
%       (nonzero only for r < h(t))
%   v - numerical solution matrix for v; v(i,j) is the value at time t(i) and position r(j)
%   h - vector of the free boundary positions over time
%   t - time vector
%   r - spatial grid vector

% Create spatial and temporal grids
r = 0:dr:H;
t = 0:dt:T;
Nr = length(r);
Nt = length(t);

% Initialize solution matrices
u = zeros(Nt, Nr);
v = zeros(Nt, Nr);
h = zeros(Nt, 1);

% Set initial conditions
h(1) = h0;
for j = 1:Nr
    if r(j) <= h0
        u(1,j) = umax * exp( - (r(j)^2)/(2*sigma^2) );
    else
        u(1,j) = 0;
    end
    % Set initial distribution for v (here simply a constant, can be modified as needed)
    v(1,j) = 0.1;  
end

% Adjust the initial front to the grid point closest to h0
[~, idx] = min(abs(r - h0));
h(1) = r(idx);

% Time-stepping loop
for i = 1:Nt-1
    % Current free boundary position
    h_current = h(i);
    % Find the grid index corresponding to h_current
    [~, hidx] = min(abs(r - h_current));
    
    % Initialize temporary variables for u and v at the next time step
    u_new = u(i,:);
    v_new = v(i,:);
    
    % --- Update u ---
    % Update u only in the region 0 < r < h(t) (u remains 0 for r >= h(t))
    for j = 2:hidx-1
        % Compute the diffusion term for u (unit diffusion coefficient) including radial term (assuming 2D symmetry)
        if r(j) > 0
            lap_u = (u(i,j+1) - 2*u(i,j) + u(i,j-1))/dr^2 ...
                  + (u(i,j+1) - u(i,j-1))/(2*dr)/r(j);
        else
            % At r=0, use symmetry: u_r(0) = 0
            lap_u = 2*(u(i,2) - u(i,1))/dr^2;
        end
        % Reaction term for u: u*(1 - u - v)
        reaction_u = u(i,j) * (1 - u(i,j) - v(i,j));
        % Time stepping update
        u_new(j) = u(i,j) + dt * ( lap_u + reaction_u );
    end
    % At r=0, enforce zero-flux boundary condition
    u_new(1) = u_new(2);
    % For r >= h(t), u remains 0
    u_new(hidx:end) = 0;
    
    % --- Update v ---
    % v is updated over the entire spatial domain (0 <= r <= H)
    for j = 2:Nr-1
        % Compute the diffusion term for v (with diffusion coefficient d)
        if r(j) > 0
            lap_v = d * (v(i,j+1) - 2*v(i,j) + v(i,j-1))/dr^2 ...
                  + d * (v(i,j+1) - v(i,j-1))/(2*dr)/r(j);
        else
            lap_v = d * 2*(v(i,2) - v(i,1))/dr^2;
        end
        % Reaction term for v: r*v*(1 - v - u)
        reaction_v = r(j) * v(i,j) * (1 - v(i,j) - u(i,j));
        % Time stepping update
        v_new(j) = v(i,j) + dt * ( lap_v + reaction_v );
    end
    % At r=0, enforce zero-flux boundary condition for v
    v_new(1) = v_new(2);
    % At r=H, use zero-flux (or Dirichlet) condition; here we choose zero-flux
    v_new(end) = v_new(end-1);
    
    % --- Update the free boundary position h(t) ---
    % Use a one-sided difference to approximate u_r at r = h(t)
    if hidx > 1
        ur_boundary = (u_new(hidx-1) - u_new(hidx))/dr;
        % To reduce numerical oscillations, apply a simple smoothing
        if hidx > 3
            ur_smooth = (3*ur_boundary + 2*(u_new(hidx-2)-u_new(hidx-1))/dr + ...
                         (u_new(hidx-3)-u_new(hidx-2))/dr)/6;
            ur_boundary = ur_smooth;
        end
        % Free boundary condition: h'(t) = -mu * u_r(t, h(t))
        h_new = h(i) - mu * ur_boundary * dt;
        % Ensure that the free boundary does not retreat
        h(i+1) = max(h(i), h_new);
    else
        h(i+1) = h(i);
    end
    % Adjust the free boundary to the grid point closest to h(i+1)
    [~, new_hidx] = min(abs(r - h(i+1)));
    h(i+1) = r(new_hidx);
    
    % Update the solutions for u and v
    u(i+1,:) = u_new;
    v(i+1,:) = v_new;
    
    % Stop simulation if the free boundary reaches the domain boundary
    if h(i+1) >= H - dr
        u = u(1:i+1,:);
        v = v(1:i+1,:);
        h = h(1:i+1);
        t = t(1:i+1);
        break;
    end
end

% Plot the evolution of the free boundary over time
figure;
plot(t, h, 'LineWidth', 2);
xlabel('Time t');
ylabel('Front position h(t)');
title('Evolution of the spreading front');

% Plot u and v profiles at selected time instants
figure;
times_to_plot = round(linspace(1, length(t), min(5, length(t))));
for i = 1:length(times_to_plot)
    idx = times_to_plot(i);
    subplot(2,1,1);
    plot(r, u(idx,:), 'LineWidth', 1.5);
    hold on;
    title('Distribution of u(t,r)');
    xlabel('r');
    ylabel('u');
    
    subplot(2,1,2);
    plot(r, v(idx,:), 'LineWidth', 1.5);
    hold on;
    title('Distribution of v(t,r)');
    xlabel('r');
    ylabel('v');
end
legend(arrayfun(@(x) sprintf('t = %.2f', t(x)), times_to_plot, 'UniformOutput', false));
hold off;

end
