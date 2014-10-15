import numpy

def lin_interp(u, ax, n1, n2, da, nsize):
    """ 
	function implements linear interpolation:
	fx = f1 + (f2-f1)/(a2-a1)*(ax-a1)
    
    parameters
    ----------
    u : FD solution at the time-step t(i) (i.e. just after)
	ax : grid point just before a[x]
	n1 : grid point just before a[x]
	n2 : grid point just after a[x]
    da : increment in A dimention
	nsize: size of the spatial grid in S dimention
    
    returns
    -------
    vn : interpolated value at the time-step t(i)^- (i.e. just before)
    """
    for i in range(1, nsize):
		ux[i, k] = u[i, n1] + (u[i, n2] - u[i, n1]) / (n2 - n1) * (ax / da - n1)

	return ux
	
def solve_lin_system(Amatrix, b):
    """ 
	function implements solution of linear system Ax  = b:
    
    parameters
    ----------
	Amatrix : 
	b : 
  
    returns
    -------
    vn : interpolated value at the time-step t(i)^- (i.e. just before)
    """
	y = Amatrix*x + b 
	
	return x	
	
# initial parameters
T = 30 / 365
r = 0.05
sigma = 0.40
K = 100.12
S0 = 100.0

# set-up containers
dt = 0.01  # T - time
da = 0.02  # A - average price
ds = 0.03  # S - spot price

# pricing region (Smin;Smax)
Smax = K + 3*sigma*S0
Smin = max(K - 3*sigma*S0,0.0)

nt = int(T / dt) + 1  # T - time dimension
na = int((Smax - Smin) / da) + 1  # A = (Smin;Smax) - average price dimension
ns = int((Smax - Smin) / ds) + 1  # S = (Smin;Smax) - spot price dimension

v = numpy.zeros((ns, na))
Amatrix = numpy.zeros((ns-1, ns-1))
b = numpy.zeros((ns-1))

# set-up I.C. at time t = T : v(A,S,T) = max(A-K,0.0)
for k in range(1, na):
    for i in range(1, ns):
        v[i, k] = max(na * k - K, 0.0)

		# we are stepping backwards:
for n in range(nt, 1, -1):
	
	# # set BC-1: if Smin << K, then A << K, so v(A,Smin,t) = max(A-K,0.0) = 0.0
    # for k in range(1, na):
		# v[0, k] = 0.0 
	
	# # BC-2: if Smax >> K, d^2V/dS^2 == 0, so dV/dt + rSdV/dS + 0.5*sigma^2*0.0 = rV	
	# for k in range(1, na):
		# v[ns,k] = 0.0 
	
	v_old = v.copy()
    for k in range(1, na):
        # a[x] = a[k] + (s[j]-a[k])/t(n)
        ax = da * k + (ds * j - da * k) / (dt * n)

        # determine grid point just before and just after a[x]
        nx = ax / da
        n1 = int(nx)
        n2 = n1 + 1

        # implement jump condition by using linear interpolation:
        # fx = f1 + (f2-f1)/(a2-a1)*(ax-a1)
		vn = lin_interp(v_old, ax, n1, n2, da, ns)

    # implement Black-Scholes equation backward time stepping
	# prepare matrix A and vector b

	### i == 0
	# set BC-1: if Smin << K, then A << K, so v(A,Smin,t) = max(A-K,0.0) = 0.0
	v[0] = 0.0 
        
	### 0 < i < ns
	# implement Black-Scholes discretization:
	# v(tn+1,i+1) + (1 -alpha + 2*gamma)/(beta + gamma)*v(tn+1,i) + v(tn+1,i-1) = b
		
	# where		
	# alpha = 0.5*r*dt
	# beta = 0.5*r*S(i)*dt/dx
	# gamma = 0.5*sigma^2*S(i)^2*dt/(dx)^2
	# b = v(tn,i+1)*(-1) + v(tn,i)*(1 + alpha + 2*gamma)/(beta + gamma) + v(tn,i-1)*(beta - gamma)/(beta + gamma) 
	for i in range(1, ns-1):
		Amatrix[i-1,i] = 1.0
		Amatrix[i,i] = (1 -alpha + 2*gamma)/(beta + gamma)
		Amatrix[i+1,i] = 1.0
		b[i] = v[i+1]*(-1) + v[i]*(1 + alpha + 2*gamma)/(beta + gamma) + v[i-1]*(beta - gamma)/(beta + gamma)
		
	### i == ns
	# BC-2: if Smax >> K, d^2V/dS^2 == 0, so dV/dt + rSdV/dS + 0.5*sigma^2*0.0 = rV			
	# v(tn+1,i) + beta*v(tn+1,i+1) - beta*v(tn+1,i) = c
		
	# where
	# beta = 0.5*r*S(i)*dt/dx
	# c = r*v(tn,i)*dt - beta*v(tn,i+1) +beta*v(tn,i)
	Amatrix[ns-1,ns] = beta
	Amatrix[ns,ns] = -beta
	b[ns] = c
		
	v_out = solve_lin_system(Amatrix, b)
	
print 'run everything'