from cmath import sqrt
from solver import Parameter, ModelTwoQuadSolver, MODEL_2_QUAD, _UNDEFINED_CASE


def valid_rhos(tau, a, s, cr, cn, delta):
    rho_1 = (8*2**(1/3)*a**2*(-1 + delta)**2 + 4*a*(-1 + delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2) + 3 * sqrt(3) * sqrt(-a**4 *(-1 + delta)**4 *tau *(-1 + delta - pr + wn)**2 * (32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2)))**(1/3) + (2 * a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2) + 6 * sqrt(3) * sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr + wn)**2 * (32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2)))**(   2/3))/(12 * a * (-1 + delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2) + 3 * sqrt(3)* sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr + wn)**2 * (32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2)))**(1/3))
    rho_2 = (-16 * (-2) ** (1 / 3) * a ** 2 * (-1 + delta) ** 2 + 8 * a * (-1 + delta) * (
    a ** 2 * (-1 + delta) ** 2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2) + 3 * sqrt(3) * sqrt(
        -a ** 4 * (-1 + delta) ** 4 * tau * (-1 + delta - pr + wn) ** 2 * (
        32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2))) ** (1 / 3) + 1j * (1j + sqrt(3)) * (
             2 * a ** 2 * (-1 + delta) ** 2 * (
             16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2) + 6 * sqrt(3) * sqrt(
                 -a ** 4 * (-1 + delta) ** 4 * tau * (-1 + delta - pr + wn) ** 2 * (
                 32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2))) ** (2 / 3)) / (
            24 * a * (-1 + delta) * (
            a ** 2 * (-1 + delta) ** 2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2) + 3 * sqrt(
                3) * sqrt(-a ** 4 * (-1 + delta) ** 4 * tau * (-1 + delta - pr + wn) ** 2 * (
            32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2))) ** (1 / 3))
    rho_3 = (8 * 1j *2**(1/3) * (1j + sqrt(3)) * a**2 * (-1 + delta)**2 +   8 * a * (-1 +   delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2) +   3 * sqrt(3)*   sqrt(-a**4 * (-1 + delta)**4 *tau * (-1 + delta - pr +   wn)**2 * (32 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2)))**(  1/3) + (-1 -   1j * sqrt(3)) * (2 * a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2) +   6 * sqrt(3)*   sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr +   wn)**2 * (32 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2)))**(  2/3))/(24 * a *(-1 +   delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2) +   3 * sqrt(3)*   sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr +   wn)**2 * (32 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2)))**(1/3))
    valids = []
    for rho in (rho_1, rho_2, rho_3):
        if abs(rho.imag) < .000001 and rho.real > 1:
            valids.append(rho.real)
    return valids


par = Parameter(MODEL_2_QUAD, tau=0.09, a=0.7, s=0.04000000000000001, cn=0.1, cr=0.04000000000000001, delta=0.7956)
wn=0.5479085657786293
pr=0.4883607806652195

rhos = valid_rhos(par.tau, par.a, par.s, par.cr, par.cn, par.delta)

# case 1:
print(ModelTwoQuadSolver.retailer_profit(par, wn, pr, rhos[0]))





def _retailer_decision(par, wn, pr, case=_UNDEFINED_CASE):
    assert type(wn) == float and type(pr) == float
    tau, a, s, cr, cn, delta = par.tau, par.a, par.s, par.cr, par.cn, par.delta
    # rho_1 = (8*2**(1/3)*a**2*(-1 + delta)**2 + 4*a*(-1 + delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2) + 3 * sqrt(3) * sqrt(-a**4 *(-1 + delta)**4 *tau *(-1 + delta - pr + wn)**2 * (32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2)))**(1/3) + (2 * a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2) + 6 * sqrt(3) * sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr + wn)**2 * (32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2)))**(   2/3))/(12 * a * (-1 + delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2) + 3 * sqrt(3)* sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr + wn)**2 * (32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn)**2)))**(1/3))
    rho_1 = (-16 * (-2) ** (1 / 3) * a ** 2 * (-1 + delta) ** 2 + 8 * a * (-1 + delta) * (
    a ** 2 * (-1 + delta) ** 2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2) + 3 * sqrt(3) * sqrt(
        -a ** 4 * (-1 + delta) ** 4 * tau * (-1 + delta - pr + wn) ** 2 * (
        32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2))) ** (1 / 3) + 1j * (1j + sqrt(3)) * (
             2 * a ** 2 * (-1 + delta) ** 2 * (
             16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2) + 6 * sqrt(3) * sqrt(
                 -a ** 4 * (-1 + delta) ** 4 * tau * (-1 + delta - pr + wn) ** 2 * (
                 32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2))) ** (2 / 3)) / (
            24 * a * (-1 + delta) * (
            a ** 2 * (-1 + delta) ** 2 * (16 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2) + 3 * sqrt(
                3) * sqrt(-a ** 4 * (-1 + delta) ** 4 * tau * (-1 + delta - pr + wn) ** 2 * (
            32 * a * (-1 + delta) - 27 * tau * (-1 + delta - pr + wn) ** 2))) ** (1 / 3))
    # rho_1 = (8 * 1j *2**(1/3) * (1j + sqrt(3)) * a**2 * (-1 + delta)**2 +   8 * a * (-1 +   delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2) +   3 * sqrt(3)*   sqrt(-a**4 * (-1 + delta)**4 *tau * (-1 + delta - pr +   wn)**2 * (32 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2)))**(  1/3) + (-1 -   1j * sqrt(3)) * (2 * a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2) +   6 * sqrt(3)*   sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr +   wn)**2 * (32 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2)))**(  2/3))/(24 * a *(-1 +   delta) * (a**2 * (-1 + delta)**2 * (16 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2) +   3 * sqrt(3)*   sqrt(-a**4 * (-1 + delta)**4 * tau * (-1 + delta - pr +   wn)**2 * (32 * a * (-1 + delta) -   27 * tau * (-1 + delta - pr + wn)**2)))**(1/3))
    # pn_1  = .5 * (1+wn-delta+pr)
    rho_2 = 1
    # pn_2 = .5 * (1+wn-delta+pr)
    valids = []

    if type(rho_1) == complex or type(rho_1) == np.complex128:
        if rho_1.imag < .000001 and rho_1.real > 1:
            rho_1 = rho_1.real
        else:
            rho_1 = -1

    for pn, rho in ((pn_1, rho_1), (pn_2, rho_2)):
        qn = 1 - (pn - pr) / (1 - delta)
        qr = (pn - pr) / (1 - delta) - pr / delta
        if rho == 0: continue
        ret_profit = qn * (pn - wn) * (1 - par.tau / rho) - par.a * (rho - 1) ** 2
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
