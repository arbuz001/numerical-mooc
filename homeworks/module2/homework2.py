import numpy
import matplotlib.pyplot as plt

from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

nx = 51  # try changing this number from 41 to 81 and Run All ... what happens?
L = 11.0
dx = L / (nx - 1)
dt = .001
nt = int(20 / 60.0 / dt)
Vmax = 80.
Rhomax = 250.

x = numpy.linspace(0, L, nx)
rho = numpy.ones(nx) * 10
rho[10:20] = 50.0  # as per our I.C.s

# print(rho)
# plt.plot(x, rho, color='#003366', ls='--', lw=3)
# plt.show()

for n in range(1, nt):
    rhon = rho.copy()
    rho[0] = 10.0
    for i in range(1, nx):
        rho[i] = rhon[i - 1] - dt / dx * Vmax * (rhon[i] - rhon[i - 1]) * (1 - rhon[i] / Rhomax)

plt.plot(x, rho, color='#003366', ls='--', lw=3)
# plt.ylim(0, 2.5);
plt.show()

print 'run everything'