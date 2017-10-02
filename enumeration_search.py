import numpy as np
from math import sqrt
from solver import Parameter, MODEL_2, DecisionVariables, Solution

STEP_SIZE_WN = 0.01
STEP_SIZE_PR = 0.01

class EnumerationSearch:
    def __init__(self, model):
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
        
def main():
    search = EnumerationSearch(MODEL_2)
    par = Parameter(MODEL_2, tau=0.3, a=.01, s=.1, cr=.02, cn=.04, delta=.7)
    sol = search.search(par)
    print(sol)
    
if __name__ == '__main__':
    main()