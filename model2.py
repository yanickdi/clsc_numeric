from math import sqrt
import numpy as np
import sys
from scipy.optimize import minimize

from solver import ImprovedGridSearch2D, _CASE_ONE, _CASE_TWO, _UNDEFINED_CASE, MODEL_2, DecisionVariables, Solution, Parameter, SolverProxy
from history import HillClimber

class ModelTwoSolver:
    """ This class solves Model Two by first using the GridSearch - then a local search algorithm (nelder-mead) """

    @staticmethod
    def solve(par):
        if par.a == 0: return None
        case_solutions = []
        # we have to check both cases:
        for case in (_CASE_ONE, _CASE_TWO):
            startP = ModelTwoGridSearch.search(par, iter=1, raster_size=51, case=case)
            start_vec = [startP.dec.wn, startP.dec.pr]
            # improve wn, pr
            result = minimize(ModelTwoSolver._minimize_func, args={'par': par, 'case': case}, x0=start_vec,
                              method='Nelder-Mead', options={'xatol': 0.00000000001, 'fatol': 0.000000001})
            #result = HillClimber.minimize(ModelTwoSolver._minimize_func, args={'par': par, 'case': case}, x0=start_vec,
            #                  iterations=20, step_sizes=0.001)
            # build a new solution using wn, pr
            sol = ModelTwoSolver.getSolution(par, float(result.x[0]), float(result.x[1]))
            # assert that wn/pr combination lies really in case `case`
            if sol.case != case:
                print('here')
            assert sol.case == case
            case_solutions.append(sol)

        if len(case_solutions) >= 1:
            return max(case_solutions, key=lambda k: k.profit_man)
        else:
            return None

    @staticmethod
    def _minimize_func(x, args):
        wn, pr = float(x[0]), float(x[1])
        par = args['par']
        case = args['case']
        man_profit, data = ModelTwoSolver.profit(par, wn, pr, case)
        if man_profit is None:
            return sys.maxsize
        return -man_profit

    @staticmethod
    def getSolution(par, wn, pr, case=_UNDEFINED_CASE):
        man_profit, sol_tuple = ModelTwoSolver.profit(par, wn, pr, case)
        if sol_tuple is None: return None
        dec = DecisionVariables(MODEL_2, pn=sol_tuple[2], pr=sol_tuple[1],
                                wn=sol_tuple[0], rho=sol_tuple[3], qn=sol_tuple[4], qr=sol_tuple[5])
        case = _CASE_ONE if sol_tuple[3] > 1 else _CASE_TWO
        sol = Solution(dec, sol_tuple[6], sol_tuple[7], case)
        return sol

    @staticmethod
    def _retailer_decision(par, wn, pr, case=_UNDEFINED_CASE):
        assert type(wn) == float and type(pr) == float
        tau, a, s, cr, cn, delta = par.tau, par.a, par.s, par.cr, par.cn, par.delta
        rho_1 = .5 * (1 - delta + pr - wn)* sqrt(tau/(a*(1 - delta)))
        pn_1  = .5 * (1+wn-delta+pr)
        rho_2 = 1
        pn_2 = .5 * (1+wn-delta+pr)
        valids = []

        for pn, rho in ((pn_1, rho_1), (pn_2, rho_2)):
            qn = 1 - (pn - pr) / (1 - delta)
            qr = (pn - pr) / (1 - delta) - pr / delta
            if rho == 0: continue
            ret_profit = qn * (pn - wn) * (1 - par.tau / rho) - par.a * (rho - 1)
            if (0 <= qr <= (par.tau / rho) * qn) and rho >= 1 and ret_profit >= 0:
                valids.append([pn, rho, qn, qr, ret_profit])
        if len(valids) > 0:
            best_valid = max(valids, key=lambda k: k[4])
            # only return solution if retailer decision is in case `case`
            best_rho = best_valid[1]
            if case == _CASE_ONE:
                return best_valid if best_rho > 1 else None
            if case == _CASE_TWO:
                return best_valid if best_rho == 1 else None
            if case == _UNDEFINED_CASE:
                return best_valid
        else:
            return None

    @staticmethod
    def profit(par, wn, pr, case=_UNDEFINED_CASE):
        ret_dec = ModelTwoSolver._retailer_decision(par, wn, pr, case=case)
        if ret_dec is None: return None, None
        pn, rho, qn, qr, ret_prof = ret_dec
        man_profit = qn * (wn * (1 - par.tau / rho) - par.cn) + qr * (pr - par.cr) + ((par.tau / rho) * qn - qr) * par.s
        return man_profit, (wn, pr, pn, rho, qn, qr, man_profit, ret_prof)


class ModelTwoGridSearch:
    @staticmethod
    def search(par, iter=1, raster_size=51, case=_UNDEFINED_CASE):
        if par.a == 0: return None
        wn_range = [0, 1]
        pr_range = [0, 1]
        raster_start_size = raster_size
        sol_tuple = ImprovedGridSearch2D.maximize(ModelTwoGridSearch._grid_search_func,
                                                  {'par': par, 'case': case}, wn_range, pr_range, raster_size, iter,
                                                  raster_start_size=raster_start_size)
        if sol_tuple is None: return None
        dec = DecisionVariables(MODEL_2, pn=sol_tuple[2], pr=sol_tuple[1],
                                wn=sol_tuple[0], rho=sol_tuple[3], qn=sol_tuple[4], qr=sol_tuple[5])
        case = _CASE_ONE if sol_tuple[3] > 1 else _CASE_TWO
        sol = Solution(dec, sol_tuple[6], sol_tuple[7], case)
        return sol

    @staticmethod
    def _grid_search_func(options, wn, pr):
        par, case = options['par'], options['case']
        return ModelTwoSolver.profit(par, wn, pr, case)


def plottiplotti():
    proxy = SolverProxy()
    a = [float(a) for a in np.linspace(0.0, 0.01, num=48 + 2)]
    x = np.zeros((len(a), 1))
    y = np.zeros((len(a), 1))
    dist = 0
    for i, act_a in enumerate(a):
        if act_a == 0:
            x[0] = None
            y[0] = None
            continue
        par = Parameter(MODEL_2, tau=0.09, a=act_a, s=0.04000000000000001, cn=0.1, cr=0.04000000000000001, delta=0.7956)
        solP = ModelTwoSolver.solve(par)
        solA = proxy.calculate(par)
        x[i] = solP.dec.wn
        y[i] = solA.dec.wn
        dist += abs(solP.dec.wn - solA.dec.wn)
    from matplotlib import pyplot as plt
    plt.plot(a, x, label='search')
    plt.plot(a, y, label='analy')
    plt.legend()
    plt.show()
    print(dist)


if __name__ == '__main__':
    plottiplotti()