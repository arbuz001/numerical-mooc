% model parameters:
g = 9.81; % gravity in m s^{-2}
ms = 50.0; % weight of the rocket shell kg   
rho = 1.091; % average air density (assumed constant throughout flight) kg m^{-3}
r = 0.5; % maximum cross section linear dimention in m 
a1 = 20.0; % fuel burn rate in kg s{-1}
ve = 325.0; % exhaust speed in m s^{-1}
Cd = 0.15; % drag coefficient
A = pi*r^2; % maximum cross sectional area of the rocket in m^{2}

%%% set initial conditions %%%
h0 = 0.0; % initial height
v0 = 0.0; % start velocity
mp0 = 100.0; % initial weight of the rocket propellant in kg

T = 40.0;	% final time

dt = 0.1; % time increment
N = ceil(T/dt) + 1; % number of time-steps
t = linspace(0.0, T, N); % time discretization

u = zeros(N, 3);

u(1,:) = [h0, v0, mp0];
	
for n = 1:N
 u(n+1,:) = euler_step(u(n,:), dt, g, a1, ve, ms, rho, A, Cd);
end;

h_path = u(:,1);
v_path = u(:,2);
m_path = u(:,3);

% plot(h_path)

% t_x = 3.2;
% n_x = ceil(t_x/T*N);
% m_x = u(n_x,3)

% idx_max_v = find(v_path == max(v_path));
% t_max_v = idx_max_v/N*T
% h_max_v = h_path(idx_max_v) 

% max_h_path = max(h_path)
% idx_max_h = find(h_path == max(h_path))
% t_max_h = idx_max_h/N*T

idx_impact = min(find(h_path < 0.0)) - 1
t_impact = idx_impact/N*T
v_impact = v_path(idx_impact)

v = (h_path(idx_impact)- h_path(idx_impact-1))/dt
