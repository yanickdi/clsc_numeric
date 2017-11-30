import sys

from solver import ModelTwoQuadGridSearch, Parameter, MODEL_2_QUAD, ModelTwoSolver
from scipy.optimize import minimize, fmin_tnc

def profit(par, wn, pr, case):
    ret_dec = ModelTwoQuadGridSearch._retailer_decision(par, wn, pr, case=case)
    print(ret_dec)
    #sys.exit()
    if ret_dec is None: return None, None
    pn, rho, qn, qr, ret_prof = ret_dec
    man_profit = qn*(wn*(1-par.tau/rho)-par.cn) + qr*(pr-par.cr)+((par.tau/rho)*qn - qr)*par.s
    return man_profit, (wn, pr, pn, rho, qn, qr, man_profit, ret_prof)
    
def minimize_func(x, args):
    wn, pr = x[0], x[1]
    par = args['par']
    args['search_points'].append((wn, pr))
    man_profit, data = profit(par, wn, pr, 0)
    print(wn, pr, man_profit)
    if man_profit is None:
        return sys.maxsize
    return -man_profit
    
if __name__ == '__main__':
    par = Parameter(MODEL_2_QUAD, tau=0.09, a=0.0006122448979591836, s=0.04000000000000001, cn=0.1, cr=0.04000000000000001, delta=0.7956)
    gridsearch = ModelTwoQuadGridSearch()
    #start_sol = gridsearch.search(par, raster_size=201)
    #start_vec = [start_sol.dec.wn, start_sol.dec.pr]
    start_vec = [0.5, 0.4]
    #start_vec = [0.5436073773054079, 0.4889027821194123] #besserer startpunkt
    args = {'par': par, 'search_points': []}
    #sol_opt = ModelTwoSolver.solve_analytical(par)
    #opt_vec = [sol_opt.dec.wn, sol_opt.dec.pr]
    result = minimize(minimize_func, args=args, x0=start_vec, method='nelder-mead')#,
    #    options={'fatol': 0.000000000001, 'xatol': 0.00000000001 })
    from matplotlib import pyplot as plt
    x, y = [], []
    for p1, p2 in args['search_points']:
        x.append(p1)
        y.append(p2)
    plt.plot(x, y, 'o', markersize=4, label='search points')
    #plt.plot(opt_vec[0], opt_vec[1], label='global opt', marker='8', color='magenta', linestyle='')
    plt.plot(start_vec[0], start_vec[1], marker='o', linestyle='', color='green', label='starting point')
    plt.plot(result.x[0], result.x[1], marker='x', color='red', label='nelder-mead result')
    plt.legend()
    plt.xlabel('wn')
    plt.ylabel('pr')
    plt.show()
    print(result)