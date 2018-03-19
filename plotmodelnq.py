import matplotlib.pyplot as plt
import numpy as np

from solver import MODEL_1_QUAD, Parameter, ModelOneQuadGridSearch

def profitm(par, wn, pn, rho):
    qn = 1 - pn
    return qn*(wn*(1 - par.tau/rho) - par.cn + (par.tau/rho)*par.s)
    
def profitr(par, wn, pn, rho):
    qn = 1 - pn
    return qn*(pn - wn)*(1 - par.tau/rho) - par.a*(rho**2 - 1)

def solvenq_case_2(par):
    # solution rho = 1:
    wn = (-1 - par.cn + par.tau + par.s * par.tau)/(2 * (-1 + par.tau))
    pn = (1 + wn)/2
    rho = 1
    return wn, profitm(par, wn, pn, rho)
  
  
par = Parameter(MODEL_1_QUAD, tau=0.15, cn=0.1, s=0.05, a=0.01)
    
all_wn = np.linspace(0, 1, 100)
# case rho=1:
profit = np.zeros(len(all_wn))
for i, wn in enumerate(all_wn):
    pn = (1 + wn)/2
    rho = 1
    if profitr(par, wn, pn, rho) > 0:
        profit[i] = profitm(par, wn, pn, rho)
    else:
        profit[i] = np.NaN
    
plt.plot(all_wn, profit)
plt.plot(solvenq_case_2(par)[0], solvenq_case_2(par)[1], linestyle='', marker='o')


# case one, rho > 1:
profit = np.zeros(len(all_wn))
for i, wn in enumerate(all_wn):
    wn = float(wn)
    pn = (1 + wn)/2
    rho = -(((-1)**(1/3) * par.tau**(1/3) * (-1 + wn)**(2/3))/(2 *par.a**(1/3)))
    if rho.imag < 0.00001: rho = rho.real
    if rho >= 1 and profitr(par, wn, pn, rho) > 0:
        profit[i] = profitm(par, wn, pn, rho)
    else:
        profit[i] = np.NaN
plt.plot(all_wn, profit)


# show solver optimum:
search = ModelOneQuadGridSearch()
sol = search.search(par)
plt.plot(sol.dec.wn, sol.profit_man, linestyle='', marker='o')
plt.plot(sol.dec.wn, profitm(par, sol.dec.wn, sol.dec.pn, sol.dec.rho), marker='*')
print(sol.profit_man, sol.dec.rho)

plt.show()