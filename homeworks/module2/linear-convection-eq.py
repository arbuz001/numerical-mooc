import numpy
import matplotlib.pyplot as plt

from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

nx = 41  # try changing this number from 41 to 81 and Run All ... what happens?
dx = 2. / (nx - 1)
nt = 25
dt = .02
c = 1.  # assume wavespeed of c = 1

u = numpy.ones(nx)  # numpy function ones()
u[.5 / dx: 1 / dx + 1] = 2  # setting u = 2 between 0.5 and 1 as per our I.C.s
print(u)

# plt.plot(numpy.linspace(0,2,nx), u, color='#003366', ls='--', lw=3)
# plt.ylim(0,2.5);
# plt.show()

for n in range(1, nt):
    un = u.copy()
    for i in range(1, nx):
        u[i] = un[i] - c * dt / dx * (un[i] - un[i - 1])

plt.plot(numpy.linspace(0, 2, nx), u, color='#003366', ls='--', lw=3)
plt.ylim(0, 2.5);
plt.show()

print 'run everything'