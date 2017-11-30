from solver import ModelTwoQuadSolver, Parameter, MODEL_2_QUAD, ModelTwoQuadGridSearch
from scipy.optimize import minimize
import numpy as np
import sys


def main():
    par = Parameter(MODEL_2_QUAD, tau=0.09, a=0.0038775510204081634, s=0.04000000000000001, cn=0.1, cr=0.04000000000000001, delta=0.7956)
    # find any start vector of case 1
    startP = ModelTwoQuadGridSearch.search(par, iter=1, case=1)
    start_vec = [startP.dec.wn, startP.dec.pr]
    # improve wn, pr
    result = minimize(ModelTwoQuadSolver._minimize_func, args={'par': par, 'case': 1}, x0=start_vec, method='Nelder-Mead', options={'xatol': 0.000000001, 'fatol': 0.000000001})
    # build a new solution using wn, pr
    solP = ModelTwoQuadSolver.getSolution(par, float(result.x[0]), float(result.x[1]))
    print(startP.profit_man, startP.dec.wn, startP.dec.pr, startP.dec.rho)
    print(solP.profit_man, solP.dec.wn, solP.dec.pr, solP.dec.rho)


def main2():
    par = Parameter(MODEL_2_QUAD, tau=0.09, a=0.0038775510204081634, s=0.04000000000000001, cn=0.1, cr=0.04000000000000001, delta=0.7956)
    sol = ModelTwoQuadSolver.solve(par)
    print(sol.profit_man)


def test_ob_gleich():
    par_21 = Parameter(MODEL_2_QUAD, tau=0.09, a=0.004081632653061225, s=0.04000000000000001, cn=0.1, cr=0.04000000000000001, delta=0.7956)
    wn_pr_21 = (0.5429821819318537, 0.48640200062519534)
    wn_pr_20 = (0.5461081587996248, 0.4885901844326352)
    # now test with wn,pr combination of par
    sol_with_21 = ModelTwoQuadSolver.getSolution(par_21, wn_pr_21[0], wn_pr_21[1])
    sol_with_20 = ModelTwoQuadSolver.getSolution(par_21, wn_pr_20[0], wn_pr_20[0])
    assert sol_with_21.profit_man == 0.15756762319737758
    print(sol_with_21.profit_man)
    print(sol_with_20.profit_man)


def plottiplotti():
    from matplotlib import pyplot as plt
    a = [float(a) for a in np.linspace(0.0, 0.01, num=48 + 2)]
    x = np.zeros((len(a), 1))
    y = np.zeros((len(a), 1))
    for i, act_a in enumerate(a):
        if act_a == 0:
            x[0] = None
            y[0] = None
            continue
        par = Parameter(MODEL_2_QUAD, tau=0.09, a=act_a, s=0.04000000000000001, cn=0.1, cr=0.04000000000000001, delta=0.7956)
        solP = ModelTwoQuadSolver.solve(par)
        x[i] = solP.dec.wn
        y[i] = solP.profit_ret
        # plot actual a on x[i]
        #plt.text(act_a, x[i], str(act_a))

    plt.plot(a, x, label='x', marker='o')
    plt.plot(a, y, label='y')
    plt.legend()
    plt.show()

def test_db_and_plot():
    from matplotlib import pyplot as plt
    from solver import SolverProxy
    proxy = SolverProxy()
    a = [float(a) for a in np.linspace(0.0, 0.01, num=48 + 2)]
    x = np.zeros((len(a), 1))
    y = np.zeros((len(a), 1))
    for i, act_a in enumerate(a):
        if act_a == 0:
            x[0] = None
            y[0] = None
            continue
        par = Parameter(MODEL_2_QUAD, tau=0.09, a=act_a, s=0.04000000000000001, cn=0.1, cr=0.04000000000000001, delta=0.7956)
        #solP = ModelTwoQuadSolver.solve(par)
        solP = proxy.read_or_calc_and_write(par, comment='hello')
        x[i] = solP.dec.wn
        y[i] = solP.profit_ret
        # plot actual a on x[i]
        #plt.text(act_a, x[i], str(act_a))
    proxy.commit()
    plt.plot(a, x, label='x', marker='o')
    plt.plot(a, y, label='y')
    plt.legend()
    plt.savefig('figure_quad.pdf')


#main()
#plottiplotti()
test_db_and_plot()
#plottiplotti()
#test_ob_gleich()