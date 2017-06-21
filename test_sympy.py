from sympy import symbols, diff
from sympy.solvers import solve
#Parameter(MODEL_1, tau=0.1, a=0.005, s=0.0005, cn=0.01)

#tau, a, s, cn = 0.1, 0.005, 0.0005, 0.01 
tau, a, s, cn = symbols('tau a s cn')
qn, pn, wn, roh = symbols('qn pn wn roh')
lam = symbols('lamda')


retailer_profit = qn * (pn - wn) * (1 - tau/roh) - a * roh
retailer_profit = retailer_profit.subs(qn, 1 - qn)
lag_ret = retailer_profit + lam * ( 1 - roh)

print(solve([lag_ret.diff(roh).subs(lam, 0), lag_ret.diff(pn).subs(lam, 0)], roh, pn))
#ret_der_roh = diff(lag_ret, roh)
#ret_der_pn = diff(lag_ret, pn)
#x = solve([diff(lag_ret, roh).subs(lam, 0), diff(lag_ret, pn).subs(lam, 0)], roh, pn, set=True)
#print(x)