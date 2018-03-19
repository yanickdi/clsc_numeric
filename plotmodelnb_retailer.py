import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from solver import Parameter, MODEL_NB
import utils


def profitr(par, wn, b, pn, rho):
    qn = 1 - pn
    profit = qn * (pn*(1 - par.tau/rho) - wn + par.tau/ rho * b) - par.a*(rho - 1)
    if profit > 0:
        return profit
    else:
        return np.NaN
    
    
par = Parameter(MODEL_NB, tau=0.15, a=0.001, s=0, cn=0.1)
wn, b = 0.5325294961734692, 0

all_pn = np.linspace(0, 1, 50)
all_rho = np.linspace(1, 15, 100)
X = np.zeros(len(all_pn)*len(all_rho))
Y = np.zeros(len(X))
Z = np.zeros(len(X))
pos = 0
for i, j in [(pn, rho) for pn in range(len(all_pn)) for rho in range(len(all_rho))]:
    pn, rho = all_pn[i], all_rho[j]
    Z[pos] = profitr(par, wn, b, pn, rho)
    
    X[pos] = pn
    Y[pos] = rho
    pos += 1
    
    
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X, Y, Z, label='all points in case one', alpha=0.1)
ax.set_xlabel('pn'); ax.set_ylabel('rho'); ax.set_zlabel('ret_profit')

index = np.nanargmax(Z)
print('here in this case, max would be here: ', X[index], Y[index], Z[index])

pn, rho = utils.model_nb_pn_rho_case_one(par, wn, b)
prof = profitr(par, wn, b, pn, rho)
print(X[index], Y[index], Z[index])
#ax.scatter(X[index], Y[index], Z[index], color='red')
ax.scatter(pn, rho, prof, color='red', marker='*', label='optimal retailer dec')

assert wn < b + ( (1-b)**2 - 4*(par.a/par.tau))**(1/2)

rho = (((1+b)**2 -4*b) / ((wn - b)**2 + 4*(par.a/par.tau)))**(1/2)
pn  = ((1+b) + ( (1+b)**2 - 4*(b+rho**2 * (par.a/par.tau)))**(1/2) ) / 2
prof = profitr(par, wn, b, pn, rho)
print(pn, rho, profitr(par, wn, b, pn, rho))
ax.scatter(pn, rho, prof, color='green', label='laut paper')

plt.legend()
plt.show()