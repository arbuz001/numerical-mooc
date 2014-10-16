% initial parameters
T = 30.0 / 365.0;
r = 0.05;
sigma = 0.40;
K = 100.12;
S0 = 100.0;

% set-up containers
nt = 300;
ns = 100;

% pricing region (Smin;Smax)
Smax = S0 + 3 * sigma * S0;
Smin = max(S0 - 3 * sigma * S0, 0.0);

dt = T / (nt - 1);
ds = (Smax - Smin) / (ns - 1);

v = zeros(ns,1);
b = zeros(ns,1);

% interm container to store solution v() for validation purposes
v_solution = zeros(nt, ns);

% set-up I.C. at time t = T : v(A,S,T) = max(A-K,0.0)
for i = 1:ns
    v(i) = max(ds * i - K, 0);
end;

% remove when validation is completed
v_solution(nt,:) = v;

% we are stepping backwards:
for n = nt-1:-1:1
	v_old = v;

    % prepare matrix A and vector b
    Amatrix = zeros(ns, ns);

    % i == 0 -> Smin
    % set BC-1: if Smin << K, then v(Smin,t) = max(Smin-K,0.0) = 0.0
    v(1) = 0.0;
    Amatrix(1,1) = 1.0;
    b(1) = 0.0;
	
    % i == ns - 1 -> Smax
    % BC-2: if Smax >> K, d^2V/dS^2 == 0, so dV/dt + rSdV/dS + 0.5*sigma^2*0.0 = rV
    % v(tn-1,i)*(1 + alpha1 - beta1) + v(tn-1,i-1)*beta1 = c1

    % where
    % alpha1 = 0.5*r*dt
    % beta1 = 0.5*r*S(i)*dt/dx
    % c1 = v(tn,i)*(1 - alpha1 + beta1) + v(tn,i-1)*beta
    alpha1 = 0.5 * r * dt;
    beta1 = 0.5 * r * (ns * ds) * dt / ds;

    % idx = (ns - 1);
    idx = ns;
	c1 = v_old(idx) * (1.0 - alpha1 + beta1) + v_old(idx - 1) * beta1;

    Amatrix(idx, idx) = (1 - alpha1 + beta1);
    Amatrix(idx, idx-1) = -beta1;
    b(idx) = c1;

    % %% 0 < i < ns-1
    % implement Black-Scholes discretization:
    % v(tn+1,i+1) + (1 -alpha + 2*gamma)/(beta + gamma)*v(tn+1,i) + v(tn+1,i-1) = b
    % where
    % alpha2 = 1/2*r*dt
    % beta2 = 1/4*r*S(i)*dt/dx
    % gamma2 = 1/8*sigma^2*S(i)^2*dt/(dx)^2
    % c2 = v_old(i+1)*(-1) + v_old(i)*(alpha2 + 2*gamma2 - 1.0)/(beta2 + gamma2) + v_old(i-1)*(gamma2 - beta2)/(beta2 + gamma2)
    for i = 1+1:ns-1
        alpha2 = 1 / 2. * r * dt;
        beta2 = 1 / 4. * r * (ds * i) * dt / ds;
        gamma2 = 1 / 8. * sigma^2 * (ds * i)^2 * dt / ds^2;
        c2 = v_old(i + 1) * (-1) + v_old(i) * (alpha2 + 2 * gamma2 - 1.0) / (beta2 + gamma2) + v_old(i - 1) * (beta2 - gamma2) / (beta2 + gamma2);

		Amatrix(i, i - 1) = 1.0;
        Amatrix(i, i) = -(1 + alpha2 + 2 * gamma2) / (beta2 + gamma2);
        Amatrix(i,i + 1) = (gamma2 - beta2) / (beta2 + gamma2);
        b(i) = c2;
	end;
	v = linsolve(Amatrix,b);
	% v = (Amatrix)\b;
	
	% remove when validation is completed
	v_solution(n,:) = v;	
end;

option_price = get_option_price(v,S0,ds,0.0,dt)