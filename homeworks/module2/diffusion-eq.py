import numpy
import matplotlib.pyplot as plt

from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16


def linearconv(nx):
    """Solve the linear convection equation.

    Solves the equation d_t u + c d_x u = 0 where
    * the wavespeed c is set to 1
    * the domain is x \in [0, 2]
    * 20 timesteps are taken, with \Delta t computed using the CFL 0.5
    * the initial data is the hat function

    Produces a plot of the results

    Parameters
    ----------

    nx : integer
        number of internal grid points

    Returns
    -------

    None : none
    """
    dx = 2. / (nx - 1)
    nt = 20
    c = 1
    sigma = .5

    dt = sigma * dx

    u = numpy.ones(nx)
    u[.5 / dx: 1 / dx + 1] = 2

    plt.plot(numpy.linspace(0, 2, nx), u, color='#003366', ls='--', lw=3)
    plt.ylim(0, 2.5);
    plt.show()

    un = numpy.ones(nx)

    for n in range(nt):
        un = u.copy()
        u[1:] = un[1:] - c * dt / dx * (un[1:] - un[0:-1])
        u[0] = 1.0

    plt.plot(numpy.linspace(0, 2, nx), u, color='#003366', ls='--', lw=3)
    plt.ylim(0, 2.5);
    plt.show()


nx = 41
linearconv(nx)

# nx = 41  # try changing this number from 41 to 81 and Run All ... what happens?
# dx = 2. / (nx - 1)
# nt = 25
# dt = .02
# c = 1.  # assume wave speed of c = 1

# u = numpy.ones(nx)  # numpy function ones()
# u[.5 / dx: 1 / dx + 1] = 2  # setting u = 2 between 0.5 and 1 as per our I.C.s
# # print(u)
#
# plt.plot(numpy.linspace(0,2,nx), u, color='#003366', ls='--', lw=3)
# plt.ylim(0,2.5);
# plt.show()

# for n in range(1, nt):
# un = u.copy()
#     for i in range(1, nx):
#         u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1])

# for n in range(1, nt):
# un = u.copy()
#     u[1:] = un[1:]-un[1:]*dt/dx*(un[1:]-un[0:-1])
#     u[0] = 1.0

# plt.plot(numpy.linspace(0, 2, nx), u, color='#003366', ls='--', lw=3)
# plt.ylim(0, 2.5);
# plt.show()

print 'run everything'