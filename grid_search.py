import sys, os
if __name__ == '__main__': sys.path.append(os.path.abspath('..'))
import numpy as np
from cmath import sqrt
from matplotlib import pyplot as plt
from clsc_numeric.solver import Parameter, MODEL_2, MODEL_2_QUAD, DecisionVariables, Solution, _CASE_ONE, _CASE_TWO
from itertools import product

STEP_SIZE_WN = 0.01
STEP_SIZE_PR = 0.01

class ModelTwoGridSearch:
    def __init__(self):
        pass
        
    def _retailer_decision(self, par, wn, pr):
        rho_1 = .5 * (1-par.delta + pr - wn)*sqrt(par.tau/(par.a*(1-par.delta)))
        pn_1  = .5 * (1+wn-par.delta+pr)
        rho_2 = 1
        pn_2 = .5 * (1+wn-par.delta+pr)
        valids = []
        for pn, rho in ((pn_1, rho_1), (pn_2, rho_2)):
            qn = 1 - (pn - pr)/(1-par.delta)
            qr = (pn-pr)/(1-par.delta) - pr/par.delta
            if rho == 0: continue
            ret_profit = qn*(pn-wn)*(1-par.tau/rho)-par.a*(rho-1)
            if (0 <= qr <= (par.tau/rho)*qn) and rho >= 1 and ret_profit >= 0:
                valids.append([pn, rho, qn, qr, ret_profit])
        if len(valids) > 0:
            return max(valids, key=lambda k: k[4])
        else:
            return None
                
        
    def search(self, par):
        max_man_profit = -1
        #            0:wn   1:pr  2:pn 3:rho 4:qn  5:qr  6:ret_profit
        max_param = None
        for wn in np.arange(0, 1+STEP_SIZE_WN, STEP_SIZE_WN):
            for pr in np.arange(0, 1+STEP_SIZE_PR, STEP_SIZE_PR):
                ret_dec = self._retailer_decision(par, wn, pr)
                if ret_dec is None: continue
                pn, rho, qn, qr, ret_prof = ret_dec
                man_profit = qn*(wn*(1-par.tau/rho)-par.cn) + qr*(pr-par.cr)+((par.tau/rho)*qn - qr)*par.s
                if man_profit > max_man_profit:
                    max_man_profit = man_profit
                    max_param = [wn, pr, pn, rho, qn, qr, ret_prof]
        
        # if we found something, build a solution object and return
        if max_man_profit >= 0:
            dec = DecisionVariables(MODEL_2, pn=max_param[2], pr=max_param[1], wn=max_param[0],
                rho=max_param[3], qn=max_param[4], qr=max_param[5])
            case = None
            sol = Solution(dec, max_man_profit, max_param[6], None)
            return sol
        else:
            return None
    
class ModelTwoQuadGridSearch:
    def __init__(self):
        pass
        
    def _retailer_decision(self, par, wn, pr):
        tau, a, s, cr, cn, delta = par.tau, par.a, par.s, par.cr, par.cn, par.delta
        #rho_1 = (8*2**(1/3)*a**2*(-1 + delta)**2 + 4*a*(-1 + delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2) + 3 * sqrt(3) * sqrt(-a**4 *(-1 + delta)**4 *tau *(-1 + delta - pr + wn)**2 * (32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2)))**(1/3) + (2 * a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2) + 6 * sqrt(3) * sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr + wn)**2 * (32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2)))**(   2/3))/(12 * a * (-1 + delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2) + 3 * sqrt(3)* sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr + wn)**2 * (32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2)))**(1/3))
        rho_1 = (-16 * (-2)**(1/3) * a**2 * (-1 + delta)**2 +  8 * a * (-1 +  delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) -  27 * tau * (-1 + delta - pr + wn)**2) +  3 * sqrt(3) *sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr +  wn)**2 * (32 * a * (-1 + delta) -  27 * tau * (-1 + delta - pr + wn)**2)))**(1/3) +  1j * (1j + sqrt( 3)) * (2 * a**2 * (-1 + delta)**2 * (16 *a *(-1 + delta) -  27 * tau * (-1 + delta - pr + wn)**2) +  6 * sqrt(3)*  sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr +  wn)**2 * (32 * a * (-1 + delta) -  27 * tau * (-1 + delta - pr + wn)**2)))**( 2/3))/(24 * a * (-1 +  delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) -  27 * tau * (-1 + delta - pr + wn)**2) +  3 * sqrt(3)*  sqrt(-a**4 * (-1 + delta)**4*  tau * (-1 + delta - pr +  wn)**2 * (32 * a * (-1 + delta) -  27 * tau * (-1 + delta - pr + wn)**2)))**(1/3))
        #rho_1 = (8 * 1j *2**(1/3) * (1j + sqrt(3)) * a**2 * (-1 + delta)**2 +   8 * a * (-1 +   delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2) +   3 * sqrt(3)*   sqrt(-a**4 * (-1 + delta)**4 *tau * (-1 + delta - pr +   wn)**2 * (32 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2)))**(  1/3) + (-1 -   1j * sqrt(3)) * (2 * a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2) +   6 * sqrt(3)*   sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr +   wn)**2 * (32 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2)))**(  2/3))/(24 * a *(-1 +   delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2) +   3 * sqrt(3)*   sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr +   wn)**2 * (32 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2)))**(1/3))
        pn_1  = .5 * (1+wn-delta+pr)
        rho_2 = 1
        pn_2 = .5 * (1+wn-delta+pr)
        valids = []
        if type(rho_1) == complex or type(rho_1) == np.complex128:
            if rho_1.imag < .000001 and rho_1.real > 1:
                rho_1 = rho_1.real
            else:
                rho_1 = -1
        for pn, rho in ((pn_1, rho_1), (pn_2, rho_2)):
            qn = 1 - (pn - pr)/(1-delta)
            qr = (pn-pr)/(1-delta) - pr/delta
            if rho == 0: continue
            ret_profit = qn*(pn-wn)*(1-par.tau/rho)-par.a*(rho-1)**2
            if (0 <= qr <= (par.tau/rho)*qn) and rho >= 1 and ret_profit >= 0:
                valids.append([pn, rho, qn, qr, ret_profit])
        if len(valids) > 0:
            return max(valids, key=lambda k: k[4])
        else:
            return None
                
        
    def search(self, par):
        """ Returns a solution object or None if no solution found """
        max_man_profit = -1
        #            0:wn   1:pr  2:pn 3:rho 4:qn  5:qr  6:ret_profit
        max_param = None
        for wn in np.arange(0, 1+STEP_SIZE_WN, STEP_SIZE_WN):
            for pr in np.arange(0, 1+STEP_SIZE_PR, STEP_SIZE_PR):
                wn, pr = float(wn), float(pr)
                ret_dec = self._retailer_decision(par, wn, pr)
                if ret_dec is None: continue
                pn, rho, qn, qr, ret_prof = ret_dec
                man_profit = qn*(wn*(1-par.tau/rho)-par.cn) + qr*(pr-par.cr)+((par.tau/rho)*qn - qr)*par.s
                if man_profit > max_man_profit:
                    max_man_profit = man_profit
                    max_param = [wn, pr, pn, rho, qn, qr, ret_prof]
        
        # if we found something, build a solution object and return
        if max_man_profit >= 0:
            dec = DecisionVariables(MODEL_2_QUAD, pn=max_param[2], pr=max_param[1], wn=max_param[0],
                rho=max_param[3], qn=max_param[4], qr=max_param[5])
            case = _CASE_ONE if max_param[3] > 1 else _CASE_TWO
            sol = Solution(dec, max_man_profit, max_param[6], case)
            return sol
        else:
            return None

def combinations(setting):
    TAU, A, S, CR, CN, DELTA = 0, 1, 2, 3, 4, 5
    MIN, MAX, STEP_SIZE = 1, 2, 3
    tau_arr = np.arange(setting[TAU][MIN], setting[TAU][MAX]+setting[TAU][STEP_SIZE], setting[TAU][STEP_SIZE])
    a_arr   = np.arange(setting[A][MIN], setting[A][MAX]+setting[A][STEP_SIZE], setting[A][STEP_SIZE])
    s_arr   = np.arange(setting[S][MIN], setting[S][MAX]+setting[S][STEP_SIZE], setting[S][STEP_SIZE])
    for tau, a, s in product(tau_arr, a_arr, s_arr):
        for cr in np.arange(s, setting[CR][MAX]+setting[CR][STEP_SIZE], setting[CR][STEP_SIZE]):
            for cn in np.arange(cr, setting[CN][MAX]+setting[CN][STEP_SIZE], setting[CN][STEP_SIZE]):
                for delta in np.arange(cr, setting[DELTA][MAX], setting[DELTA][STEP_SIZE]):
                    yield tau, a, s, cr, cn, delta
                    
def parameter_combinations(setting):
    for tau, a, s, cr, cn, delta in combinations(setting):
        yield Parameter(MODEL_2_QUAD, tau=tau, a=a, s=s, cr=cr, cn=cn, delta=delta)
        
def plot_profits():
    solver = ModelTwoQuadGridSearch()
    
    num = 5
    all_a = np.linspace(0.0, 0.01, num=num)
    profit_man, profit_ret = np.zeros(num), np.zeros(num)
    cn = .1
    for i, a in enumerate(all_a):
        par = Parameter(MODEL_2_QUAD, tau=.09, cr=0.4*cn, s=0.4*cn, delta=.7956, cn=cn, a=a)
        sol = solver.search(par)
        if sol == None or a == 0:
            profit_man[i] = None
            profit_ret[i] = None
        else:
            profit_man[i] = sol.profit_man
            profit_ret[i] = sol.profit_ret
    
    print(profit_man[0])
    plt.plot(all_a, profit_man, label=r'$\pi_{M}^{Q}$')
    plt.plot(all_a, profit_ret, label=r'$\pi_{R}^{Q}$')
    plt.legend()
    plt.show()

def main():
    plot_profits()
    sys.exit()
    search = ModelTwoQuadGridSearch()
    par = Parameter(MODEL_2_QUAD, tau=0.4, a=.001, s=.1, cr=.1, cn=.2, delta=.7)
    sol = search.search(par)
    print('retailer prof: ', sol.profit_ret)
    print(sol.case)
    print('manufacturer prof: ', sol.profit_man)
    print(sol.dec.rho)
    sys.exit()
    setting = [['tau', 0.4, 0.4, 0.1], ['a', 0.001, 0.001, 0.01], ['s', 0.1, 0.1, 0.1], ['cr', 's', 1, 0.1], ['cn', 'cr', 0.2, 0.1], ['delta', 'cr', 1, 0.3]]
    i = 0
    for par in parameter_combinations(setting):
        sol = search.search(par)
        if sol is not None:
            print(par, sol, sol.case, sol.dec.rho)
            #print('hello')
    
if __name__ == '__main__':
    #import sys, os
    #sys.path.append(os.path.abspath('..'))
    #from clsc_numeric import grid_search
    main()