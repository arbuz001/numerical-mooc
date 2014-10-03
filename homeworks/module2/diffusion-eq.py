import numpy
import sympy

from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

from sympy import init_printing
init_printing()

from sympy.utilities.lambdify import lambdify

x, nu, t = sympy.symbols('x nu t')
phi = sympy.exp(-(x - 4 * t) ** 2 / (4 * nu * (t + 1))) + \
      sympy.exp(-(x - 4 * t - 2 * numpy.pi) ** 2 / (4 * nu * (t + 1)))
phi
print'phi: ', (phi), '\n'

phiprime = phi.diff(x)
phiprime
print'derivative of phi: ', (phiprime), '\n'

u = -2 * nu * (phiprime / phi) + 4
print(u)

ufunc = lambdify((t, x, nu), u)
print("The value of u at t=1, x=4, nu=3 is {}.".format(ufunc(1, 4, 3))), '\n'

f = sympy.Pow(sympy.cos(x), 2) * sympy.Pow(sympy.sin(x), 3) / (4 * sympy.Pow(x, 5) * sympy.exp(x))
print 'f: ', (f)
fprime = f.diff(x)
print 'derivative of f: ', (fprime)

ffunc = lambdify((x), fprime)
print("The value of derivative of f at x=2.2 is {}.".format(ffunc(2.2))), '\n'

print 'run everything'