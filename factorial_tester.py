from matplotlib import pyplot as plt
import numpy as np

from solver import Parameter, SolverProxy, MODEL_1, MODEL_2_QUAD, ModelTwoQuadSolver, ModelNBSolver


proxy = SolverProxy()

def some_shit():
    x = []
    y = []

    for a in np.linspace(0.001, 0.1, num=100):
        a = float(a)
        #a = 0.0069387755102040816
        tau = 0.15
        s = 0.04000000000000001
        cn = 0.1
        par = Parameter(MODEL_1, tau=tau, a=0.0069387755102040816, s=0.04000000000000001, cn=0.1)
        sol = proxy.calculate(par)

        if sol.case == '2': continue
        x.append(a)
        y.append(sol.dec.qn)

    plt.plot(x,y, linestyle='', marker='s')
    plt.show()


def case_quad_fixed_a():
    cn = np.linspace(0.0, 0.9, num=100)
    y = []
    for act_cn in np.linspace(0.0, 0.9, num=100):
        act_cn = float(act_cn)
        par = Parameter(MODEL_2_QUAD, tau=0.09, a=.4, s=.1, cr=.2, cn=act_cn, delta=.4)
        sol = ModelTwoQuadSolver.solve(par, 'low')
        y.append(sol.dec.qr)
    plt.plot(cn, y)
    print(y)
    plt.show()

def case_model_nb_fixed_cn():
    a = np.linspace(0.001, 0.02, num=10)
    y = []
    for act_a in a:
        act_a = float(act_a)
        par = Parameter(MODEL_2_QUAD, tau=0.09, a=act_a, s=0.4*0.1, cr=0.4*0.1, cn=0.1, delta=.7956)
        sol = ModelNBSolver.solve(par, 'low')
        print(act_a)
        y.append(sol.dec.b)
        y.append()
    plt.plot(a, y)
    plt.show()


case_model_nb_fixed_cn()