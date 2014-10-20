% initial parameters
T = 365.0 / 365.0;
r = 0.10;
sigma = 0.40;
K = 100.0;
S0 = 105.0;
Nsamples = 365;

% set-up containers
nt = 30;
ns = 10;
ng = 14;

% pricing region (Smin;Smax)
Smax = S0 + 3 * sigma * T * S0;
Smin = max(S0 - 3 * sigma * T * S0, 0.0);

% running sum region (Gmin;Gmax)
Gmax = K + Smax*Nsamples;
Gmin = Smin;

dt = T / (nt - 1);
ds = (Smax - Smin) / (ns - 1);
dg = (Gmax - Gmin) / (ng - 1);

vv = zeros(ng,ns);

% n = nt

% set-up I.C. at time t = T : v(S,g,T) = max(g/Nsample-K,0.0)
for j = 1:ng
	for i = 1:ns
		vv(j,i) = max(dg * j/ Nsamples - K, 0);
	end;
end;

% we are stepping backwards:
for n = nt-1:-1:1
	for j = 1:ng
		v_old = vv(j,:);

		% prepare matrix A and vector b
		Amatrix = zeros(ns, ns);
		b = zeros(ns,1);

		% idx == 0 -> Smin
		% BC-1: if Smin << K, d^2V/dS^2 == 0, so dV/dt + rSdV/dS + 0.5*sigma^2*0.0 = rV
		% v(tn-1,i) = v(tn,i)*(1 - alpha0 - beta0) + v(tn-1,i+1)*beta0

		% where
		% alpha0 = r*dt
		% beta0 = r*S(i)*dt/dx
		% c0 = v(tn,i)*(1 - alpha0 - beta0) + v(tn-1,i+1)*beta0
		idx = 1;
		alpha0 = r * dt;
		beta0 = r * (idx * ds) * dt / ds;

		c0 = v_old(idx) * (1.0 - alpha0 - beta0) + v_old(idx + 1) * beta0;

		Amatrix(idx, idx) = 1.0;
		b(idx) = c0;
		
		% idx == ns -> Smax
		% BC-2: if Smax >> K, d^2V/dS^2 == 0, so dV/dt + rSdV/dS + 0.5*sigma^2*0.0 = rV
		% v(tn-1,i) = v(tn,i)*(1 - alpha1 + beta1) + v(tn-1,i-1)*beta1
		
		% where
		% alpha1 = r*dt
		% beta1 = r*S(i)*dt/dx
		% c1 = v(tn,i)*(1 - alpha1 + beta1) + v(tn-1,i-1)*beta1
		
		idx = ns;
		alpha1 = r * dt;
		beta1 = r * (idx * ds) * dt / ds;

		c1 = v_old(idx) * (1.0 - alpha1 + beta1) - v_old(idx - 1) * beta1;

		Amatrix(idx, idx) = 1.0;
		b(idx) = c1;

		% %% 0 < idx < ns-1
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
		vv(j,:) = Amatrix\b;
	end;

	% implement jump condition V(S,g,ti^-) = V(S,g + S,ti)
	vv_copy = vv;

	for j = 1:ng
		for i = 1:ns
			jj = floor(j + i*ds/dg);
			% if for cases when g+S > Gmax set value of v same as for Gmax 
			if(jj+1 > ng)		
				vv(j,i) = vv_copy(ng,i);
			else
				vv(j,i) = vv_copy(jj,i) + (vv_copy(jj+1,i) - vv_copy(jj,i))*(j +i*ds/dg -jj) ;
			end;	
		end;
	end;	
end;

% validation step
csvwrite('vv.csv',vv)
mesh(vv)

Gx = 0.0;
Sx = 105.0;

% x = Smin + (1:ns)*ds;
% y = vv(1,:);

% option_price = spline(x,y,Sx);

% plot(x,y);
% hold on
% plot(Sx,option_price,'o');
% hold off

% option_price = get_option_price(vv,[Gx,Sx],[dg,ds],[Gmin,Smin])