from sympy import symbols, init_printing, solve, Rational, lambdify
from scipy.optimize import root
import numpy as np

from solver import Parameter, MODEL_2

init_printing()
par = Parameter(MODEL_2, tau=0.1, a=0.05, s=0.1, cr=0.2, cn=0.3, delta=0.5)

tau, a, s, cr, cn, delta = symbols('tau a s cr cn delta')
qn, qr = symbols('qn qr')
pn, rho = symbols('pn rho')
wn, pr = symbols('wn pr')
lam, mu1, mu2 = symbols('lam mu1 mu2')

qn_ = 1 - (pn - pr)/(1 - delta)
qr_ = (pn - pr)/(1 - delta) - pr/delta

profit_man = (wn - cn - tau/rho*wn)*qn_ + qr_*(pr - cr) + (tau/rho*qn_ - qr_)*s

# case 1:
pn_one = Rational(1,2)*(1 - delta + pr + wn)
rho_one = -((tau**Rational(1,3) * (-1 + delta - pr + wn)**Rational(2,3))/(2*(a*(-1 + delta))**Rational(1,3)))
f_pn_one = lambdify([delta, wn, pr], pn_one)
f_rho_one = lambdify([tau, a, delta, wn, pr], rho_one)

profit_man_one = profit_man.subs([(pn, pn_one), (rho, rho_one)])
lagr_man_one = (profit_man_one -mu1*(-qr_) - mu2*(qr_ - tau/rho*qn_)).subs([(pn, pn_one), (rho, rho_one)])
lag_diff_wn = lagr_man_one.diff(wn)
f_lag_diff_wn = lambdify([tau, a, s, cr, cn, delta, wn, pr, mu1, mu2], lagr_man_one.diff(wn))
lag_diff_pr = lagr_man_one.diff(pr)
lag_diff_mu1 = lagr_man_one.diff(mu1)


# case 1a:  diff_wn == 0,   diff_pr == 0,   diff_mu1 == 0,  mu2 == 0
def f_case_1a(vector, args):
    """
        vector must be [wn, pr, mu1]
        args must be a tuple (ParameterObject)
        returns a 3 dimensional scalar vector [lag_diff_wn, lag_diff_pr, lag_diff_mu1]
    """
    par = args
    repl_static = {tau : par.tau, a : par.a, s: par.s, cr : par.cr, cn : par.cn, delta : par.delta, mu2 : 0}
    repl_vars = {wn : vector[0], pr: vector[1], mu1: vector[2]}
    y = np.zeros(3)
    wn_val = lag_diff_wn.subs(repl_static).subs(repl_vars).evalf()
    print(wn_val)
    print(f_lag_diff_wn(par.tau, par.a, par.s, par.cr, par.cn, par.delta, vector[0], vector[1], vector[2], 0))
    pr_val = lag_diff_pr.subs(repl_static).subs(repl_vars).evalf()
    mu1_val = lag_diff_mu1.subs(repl_static).subs(repl_vars).evalf()
    #print(wn_val, pr_val, mu1_val)
    return y
    
result = root(f_case_1a, x0=[0.1, 0.1, 0.1], args=par)
print(result)