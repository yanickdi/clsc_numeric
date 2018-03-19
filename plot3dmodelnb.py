import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from solver import ImprovedGridSearch2D, ModelNBSolver, Parameter, MODEL_NB, _UNDEFINED_CASE

X=[]; Y=[]; Z=[]

def middleware_profit(par, wn, b, case=_UNDEFINED_CASE):
    man_profit, sol_tuple = ModelNBSolver.profit(par, wn, b)
    if man_profit is not None and man_profit > 0.025:
        X.append(wn); Y.append(b); Z.append(man_profit)
    return man_profit, sol_tuple

par = Parameter(MODEL_NB, tau=0.05, a=0.007, s=0, cn=0.5)

raster_size = 100; wn_range = [0, 1]; b_range = [0, 1]; iter=1
grid_search_func = middleware_profit

sol_tuple = ImprovedGridSearch2D.maximize(grid_search_func,
    par, wn_range, b_range, raster_size, iter)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X, Y, Z)
ax.set_xlabel('wn'); ax.set_ylabel('b'); ax.set_zlabel('man_profit')
plt.show()