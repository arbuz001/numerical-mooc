from math import sin, cos, log, ceil, pi
import numpy
import matplotlib.pyplot as plt
from matplotlib import rcParams
from numpy import int
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

def f(u):
    """Returns the right-hand side of the phugoid system of equations.
    
    Parameters
    ----------
    u : array of float
        array containing the solution at time n.
        
    Returns
    -------
    dudt : array of float
        array containing the RHS given u.
    """
    
    h	= u[0]
    v	= u[1]
    mp	= u[2]
    
    if  mp > 0 :
        return numpy.array([v,
						-g + a1*ve/(ms+mp) - 0.5*rho*v*abs(v)*A*Cd/(ms+mp),
						-a1])
    else:
        return numpy.array([v,
                        -g - 0.5*rho*v*abs(v)*A*Cd/(ms+mp),
                        0.0])

def euler_step(u, f, dt):
    """Returns the solution at the next time-step using Euler's method.
    
    Parameters
    ----------
    u : array of float
        solution at the previous time-step.
    f : function
        function to compute the right hand-side of the system of equation.
    dt : float
        time-increment.
    
    Returns
    -------
    u_n_plus_1 : array of float
        approximate solution at the next time step.
    """
    
    return u + dt * f(u)
	
def get_diffgrid(u_current, u_fine, dt):
    """Returns the difference between one grid and the fine one using L-1 norm.
    
    Parameters
    ----------
    u_current : array of float
        solution on the current grid.
    u_finest : array of float
        solution on the fine grid.
    dt : float
        time-increment on the current grid.
    
    Returns
    -------
    diffgrid : float
        difference computed in the L-1 norm.
    """
    
    N_current = len(u_current[:,0])
    N_fine = len(u_fine[:,0])
   
    grid_size_ratio = ceil(N_fine/float(N_current))
    
    diffgrid = dt * numpy.sum( numpy.abs(\
            u_current[:,2]- u_fine[::grid_size_ratio,2])) 
    
    return diffgrid
    
# model parameters:
g	= 9.81		# gravity in m s^{-2}
ms	= 50.0		# weight of the rocket shell kg   
rho	= 1.091		# average air density (assumed constant throughout flight) kg m^{-3}
r	= 0.5		# maximum cross section linear dimention in m 
a1	= 20.0		# fuel burn rate in kg s{-1}
ve	= 325.0		# exhaust speed in m s^{-1}
Cd	= 0.15		# drag coefficient
A	= pi*r**2   # maximum cross sectional area of the rocket in m^{2}

### set initial conditions ###
h0	= 0.0    # initial height
v0	= 0.0    # start velocity
mp0	= 100.0 # initial weight of the rocket propellant in kg

T   = 40.0  # final time

dt_values	=	numpy.array([0.1])
u_values	=	numpy.empty_like(dt_values, dtype=numpy.ndarray)

for i, dt in enumerate(dt_values):
    
    N = int(T/dt) + 1    # number of time-steps
    
    ### discretize the time t ###
    t = numpy.linspace(0.0, T, N)
    
    # initialize the array containing the solution for each time-step
    u = numpy.empty((N, 3))
    u[0] = numpy.array([h0, v0, mp0])

    # time loop
    for n in range(N-1):
        u[n+1] = euler_step(u[n], f, dt)   ### call euler_step() ###
    
    # store the value of u related to one grid
    u_values[i] = u
	
# get the rocket height with respect to the time
h_path = u[:,0];
v_path = u[:,1];
m_path = u[:,2];
      
# plt.plot(v)
# plt.show()

t_x = 3.2
n_x = int(t_x/T*N)
m_x = m_path[n_x]
print 'remaining fuel at 3.2 s:', m_x

v_max = max(v_path)
print 'maximum speed of the rocket im m/s:', v_max
idx_max_v = numpy.where(v_path==v_max)[0][0]
t_max_v = idx_max_v/(N+0.0)*(T+0.0)
print 'maximum speed occurs at s:', t_max_v
h_max_v = h_path[idx_max_v]
print 'height at maximum speed in m:', h_max_v

idx_impact = numpy.where(h_path < 0.0)[0][0]
t_impact = idx_impact/(N+0.0)*(T+0.0)
v_impact = v_path[idx_impact]
print 'velocity at impact in m/s:', h_max_v

print 'run everything'