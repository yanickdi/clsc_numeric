import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from solver import MODEL_2_QUAD, Parameter, ModelTwoQuadSolver


def profitm(par, wn, pr, pn, rho):
    qn = 1 - (pn - pr)/(1 - par.delta)
    qr = (pn - pr)/(1 - par.delta) - pr/par.delta
    if rho < 1 or not (0 <= qr <= (par.tau/rho)*qn):
        return np.NaN
    prof = qn*(wn - par.cn - par.tau/rho*wn) + qr*(pr - par.cr) + (par.tau/rho*qn - qr)*par.s
    if prof > 0.1:
        return prof
    else:
        return np.NaN

def pn_rho_case_two(par, wn, pr):
    pn = (1/2) * (1 - par.delta + pr + wn)
    rho = 1
    return pn, rho
    
def pn_rho_case_one(par, wn, pr):
    wn, pr = float(wn), float(pr)
    pn = (1/2) * (1 - par.delta + pr + wn)
    #rho = -((par.tau**(1/3) *(-1 + par.delta - pr + wn)**(2/3))/(2 * (par.a * (-1 + par.delta))**(1/3)))
    #rho = ((-1)**(1/3) * par.tau**(1/3) * (-1 + par.delta - pr + wn)**(2/3))/(2 * (par.a *(-1 + par.delta))**(1/3))
    rho = -(((-1)**(2/3) * par.tau**(1/3) *(-1 + par.delta - pr + wn)**(2/3))/( 2 *(par.a *(-1 + par.delta))**(1/3)))
    assert rho.imag < 0.00001
    rho = rho.real
    return pn, rho
    
    
par = Parameter(MODEL_2_QUAD, tau=0.15, cn=0.1, s=0.05, cr=0.1*0.1, delta=0.85, a=0.001)

all_wn = np.linspace(0, 1, 100)
all_pr = np.linspace(0, 1, 100)
X = np.zeros(len(all_wn)*len(all_pr))
Y = np.zeros(len(X))
Z_CASE_ONE = np.zeros(len(X))
Z_CASE_TWO = np.zeros(len(X))
pos = 0
for i, j in [(wn, pr) for wn in range(len(all_wn)) for pr in range(len(all_pr))]:
    wn, pr = all_wn[i], all_pr[j]
    pn, rho = pn_rho_case_one(par, wn, pr)
    Z_CASE_ONE[pos] = profitm(par, wn, pr, pn, rho)
    
    pn, rho = pn_rho_case_two(par, wn, pr)
    Z_CASE_TWO[pos] = profitm(par, wn, pr, pn, rho)
    
    X[pos] = wn
    Y[pos] = pr
    pos += 1
    
    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X, Y, Z_CASE_ONE, label='case one')
ax.scatter(X, Y, Z_CASE_TWO, label='case two')
ax.set_xlabel('wn'); ax.set_ylabel('pr'); ax.set_zlabel('man_profit')

print('max case one: {}'.format(np.nanmax(Z_CASE_ONE)))
print('max case two: {}'.format(np.nanmax(Z_CASE_TWO)))
# solver solution
sol = ModelTwoQuadSolver.solve(par)
ax.scatter(sol.dec.wn, sol.dec.pr, sol.profit_man)
print('solver sol  : {} (case {})'.format(sol.profit_man, sol.case))

plt.legend()
plt.show()