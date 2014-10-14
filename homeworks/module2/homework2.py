import numpy
import matplotlib.pyplot as plt

from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

L = 11.0
nx = 51
dx = L / (nx - 1)
dt = .001
n_minutes = 6.0
nt = int(n_minutes / 60.0 / dt)
Vmax = 80.
Rhomax = 250.

# as per our I.C.s
x = numpy.linspace(0, L, nx)
rho = numpy.ones(nx) * 10
rho[10:20] = 50.0

# print(rho)
# plt.plot(x, rho, color='#003366', ls='--', lw=3)
# plt.show()

# rho1 = numpy.ones(nx) * 10
# rho1[10:20] = 50.0  # as per our I.C.s

# at time t = 0
v = Vmax * (1 - rho / Rhomax)
v_min = min(v)
v_min_ms = v_min * 1000.0 / 3600.0
print 'minimum speed at time t = 0 in m/s:', v_min_ms

for n in range(1, nt):
    rhon = rho.copy()
    for i in range(1, nx):
        rho[i] = rhon[i] + Vmax/Rhomax * dt / dx * (rhon[i] - rhon[i - 1])*(2*rhon[i] - Rhomax)
        rho[0] = 10.0
		
# for n in range(1, nt):
	# rhon1 = rho.copy()
	# rho1[1:] = rhon1[1:] + Vmax/Rhomax * dt / dx * (rhon1[1:] - rhon1[:-1])*(2*rhon1[1:] - Rhomax)
	# rho1[0] = 10.0
	
# average speed
# v = Vmax * (1 - rho / Rhomax)
# v_avg = numpy.mean(v)
# v_avg_ms = v_avg * 1000.0 / 3600.0
# print 'average speed at time t =', n_minutes, 'in m/s:', v_avg_ms

# minimum speed
v = Vmax * (1 - rho / Rhomax)
v_min = min(v)
v_min_ms = v_min * 1000.0 / 3600.0
print 'minimum speed at time t =', n_minutes, 'in m/s:', v_min_ms

plt.plot(x, v, color='#003366', ls='--', lw=3)
plt.show()

print 'run everything'