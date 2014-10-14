import numpy

# initial parameters
T = 30 / 365
r = 0.05
sigma = 0.15
K = 100.12
S0 = 100.0

# set-up containers
dt = 0.01  # T - time
da = 0.02  # A - average price
ds = 0.03  # S - spot price

nt = int(T / dt) + 1  # T - time dimension
na = int(2 * K / da) + 1  # A = (-K;K) - average price dimension
ns = int(2 * K / ds) + 1  # S = (-K;K) - spot price dimension

v = numpy.zeros(ns, na)

# set-up I.C. at time t = T : v(A.S,T) = max(A-K,0.0)
for k in range(1, na):
    for i in range(1, ns):
        v[i, k] = max(na * k - K, 0.0)

for n in range(1, nt):
    # implement jump condition by using interpolation:
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
        for i in range(1, ns):
            vn[i, k] = v_old[i, n1] + (v_old[i, n2] - v_old[i, n1]) / (n2 - n1) * (ax / da - n1)

        # implement Black-Scholes equation backward time stepping
        for k in range(1, na):
            for i in range(1, ns):
                v[i, k] = vn[i, k] + 1.0 * vn[i - 1, k - 1]
                v[0] = 0.0

print 'run everything'