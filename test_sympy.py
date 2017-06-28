from sympy import symbols, diff, sin, cos, exp, lambdify, Rational
from sympy.solvers import solve
import numpy as np
from scipy.optimize import root

#from mpmath import findroot
#from mpmath.calculus.optimization import MDNewton

x1, x2, z = symbols('x1 x2 z')
y1 = exp(-exp(-(x1+x2))) - x2*(1+x1**2) + x1**Rational(1/2)
y2 =  x1*cos(x2) + x2*sin(x1) - 0.5
f1 = lambdify([x1, x2], y1)
f2 = lambdify([x1, x2], y2)

d_f1_x1 = lambdify([x1, x2], y1.diff(x1))
d_f1_x2 = lambdify([x1, x2], y1.diff(x2))
d_f2_x1 = lambdify([x1, x2], y2.diff(x1))
d_f2_x2 = lambdify([x1, x2], y2.diff(x2))

def f(vec):
    par0, par1 = vec
    print(par0, par1)
    return [f1(par0, par1), f2(par0, par1)]
    
def jac(vec):
    par0, par1 = vec
    return [
        [d_f1_x1(par0, par1), d_f1_x2(par0, par1)],
        [d_f2_x1(par0, par1), d_f2_x2(par0, par1)]
    ]

root(f, [.5, .3], jac=jac)
#findroot(f, [0,0], J=jac, solver=MDNewton)
#print(y.subs({x1 : 1, x2: -2}).evalf())
#print(f(1,2.0))
