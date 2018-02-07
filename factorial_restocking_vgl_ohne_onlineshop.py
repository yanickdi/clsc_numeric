import sys, os.path, re
from statistics import mean

from jinja2 import Template
from solver import SolverProxy, Parameter, MODEL_1, MODEL_NB, ModelNBSolver

LATEX_SUBS = (
    (re.compile(r'\\'), r'\\textbackslash'),
    (re.compile(r'([{}_#%&$])'), r'\\\1'),
    (re.compile(r'~'), r'\~{}'),
    (re.compile(r'\^'), r'\^{}'),
    (re.compile(r'"'), r"''"),
    (re.compile(r'\.\.\.+'), r'\\ldots'),
)

def escape_tex(value):
    newval = value
    for pattern, replacement in LATEX_SUBS:
        newval = pattern.sub(replacement, newval)
    return newval

class Factorial:
    tau_lo   = lambda : .005
    tau_hi   = lambda : 0.35
    tau_lh   = lambda : (.005, .35)
    
    s_lo     = lambda : 0
    s_hi     = lambda cn: cn/2
    s_lh     = lambda cn: (0, cn/2)
    
    cr_lo    = lambda : .1
    cr_hi    = lambda cn, delta : delta*(cn/2)
    cr_lh    = lambda cn, delta : (.1, delta*(cn/2))
    
    delta_lo = lambda : .5
    delta_hi = lambda : .85
    delta_lh = lambda : (.5, .85)
    
    cn_lo    = lambda : .1
    cn_hi    = lambda : .5
    cn_lh    = lambda : (.1, .5)
    
    a_lo     = lambda : .001
    a_hi     = lambda : .1
    a_lh     = lambda : (.001, .1)

    LOW, HIGH = 'LOW', 'HIGH'

    def __init__(self):
        self.__dict = {}
        self.__proxy = SolverProxy()
        self.__calc__()
        self.__fullFactorialTable()

    def __key(self, model, tau, a, s, cr, cn, delta):
        assert model in (MODEL_NB, MODEL_1)
        for val in (tau, a, s, cr, cn): assert val in (Factorial.LOW, Factorial.HIGH)
        return '{},tau={},a={},s={},cr={},cn={},delta={}'.format(model, tau, a, s, cr, cn,delta)

    def _add(self, key, par, sol):
        assert key not in self.__dict
        self.__dict[key] = {'par' : par, 'sol': sol}
        
    def _get(self, key):
        return self.__dict[key]

    def __calc__(self):
        p = Factorial
        iter = self.__state_iter
        proxy = SolverProxy()
        num = 0
        for tau, tau_state in iter([p.tau_lo(), p.tau_hi()]):
            for cn, cn_state in iter([p.cn_lo(), p.cn_hi()]):
                for s, s_state in iter([p.s_lo(), p.s_hi(cn)]):
                    for a, a_state in iter([p.a_lo(), p.a_hi()]):
                        # calculate model with restocking_fee
                        par_nb = Parameter(MODEL_NB, tau=tau, a=a, s=s, cn=cn)
                        sol_nb = proxy.read_or_calc_and_write(par_nb, resolution='high')
                        num += 1; print(num)
                        #sol_nb = ModelNBSolver.solve(par_nb, 'low')
                        key_nb = self.__key(MODEL_NB, tau_state, a_state, s_state, Factorial.HIGH, cn_state, Factorial.HIGH)
                        self._add(key_nb, par_nb, sol_nb)

                        # calculate model without online store
                        par_no = Parameter(MODEL_1, tau=tau, a=a, s=s, cn=cn)
                        sol_no = proxy.calculate(par_no)
                        key_no = self.__key(MODEL_1, tau_state, a_state, s_state, Factorial.HIGH, cn_state, Factorial.HIGH)
                        self._add(key_no, par_no, sol_no)
        proxy.commit()
                                
    def __fullFactorialTable(self):
        self.__ff_table = [[None for j in range(4)] for i in range(4)]
        lohi = (Factorial.LOW, Factorial.HIGH)
        for i, (s, a) in enumerate([(s, a) for s in lohi for a in lohi]):
            for j, (cn, tau) in enumerate([(cn, tau) for cn in lohi for tau in lohi]):
                key_nb = self.__key(MODEL_NB, tau, a, s, Factorial.HIGH, cn, Factorial.HIGH)
                key_no = self.__key(MODEL_1, tau, a, s, Factorial.HIGH, cn, Factorial.HIGH)
                val_nb = self._get(key_nb)
                val_no = self._get(key_no)
                self.__ff_table[i][j] = {'no' : val_no, 'nb' : val_nb}
                
    def __anal_line_parameters(self, table, par):
        f = Factorial
        pars = []
        if table == 'delta':
            delta = par
            for tau, cn in [(tau,cn) for tau in f.tau_lh() for cn in f.cn_lh()]:
                for a, s, cr in [(a,s,cr) for a in f.a_lh() for s in f.s_lh(cn) for cr in f.cr_lh(cn,delta)]:
                    pars.append([tau, a, s, cr, cn, delta])
        if table == 'tau':
            tau = par
            for delta, cn in [(delta,cn) for delta in f.delta_lh() for cn in f.cn_lh()]:
                for a, s, cr in [(a,s,cr) for a in f.a_lh() for s in f.s_lh(cn) for cr in f.cr_lh(cn,delta)]:
                    pars.append([tau, a, s, cr, cn, delta])
        if table == 'cn':
            cn = par
            for tau, delta in [(tau, delta) for tau in f.tau_lh() for delta in f.delta_lh()]:
                for a, s, cr in [(a,s,cr) for a in f.a_lh() for s in f.s_lh(cn) for cr in f.cr_lh(cn,delta)]:
                    pars.append([tau, a, s, cr, cn, delta])
        if table == 'a':
            a = par
            for tau, delta, cn in [(tau,delta,cn) for tau in f.tau_lh() for delta in f.delta_lh() for cn in f.cn_lh()]:
                for s, cr in [(s,cr) for s in f.s_lh(cn) for cr in f.cr_lh(cn,delta)]:
                    pars.append([tau, a, s, cr, cn, delta])
        assert len(pars) == int(2**5)
        return pars
                
    def getAnalysisLine(self, table, par):
        pars = self.__anal_line_parameters(table, par)
        profits, rhos, pns, wns, prs = [], [], [], [], []
        # calc all vals
        for par in pars:
            tau, a, s, cr, cn, delta = par
            sol_n, sol_o = self.__both_solutions(tau, a, s, cr, cn, delta)
            profits.append(sol_o.profit_man / sol_n.profit_man)
            rhos.append(sol_o.dec.rho / sol_n.dec.rho)
            pns.append(sol_o.dec.pn / sol_n.dec.pn)
            wns.append(sol_o.dec.wn / sol_n.dec.wn)
            prs.append(sol_o.dec.pr)
        mean_profits = mean(profits)*100
        mean_rhos = mean(rhos)*100
        mean_pns = mean(pns)*100
        mean_wns = mean(wns)*100
        mean_prs = mean(prs)
        return r'{:.2f}\%&{:.2f}\%&{:.2f}\%&{:.2f}\%&{:.3f}'.format(
            mean_profits, mean_rhos, mean_pns, mean_wns, mean_prs)
        
    def __both_solutions(self, tau, a, s, cr, cn, delta):
        par_n = Parameter(MODEL_1, tau=tau, a=a, s=s, cn=cn)
        sol_n = self.__proxy.calculate(par_n)
        par_o = Parameter(MODEL_2, tau=tau, a=a, s=s, cr=cr, cn=cn, delta=delta)
        sol_o = self.__proxy.calculate(par_o)
        return (sol_n, sol_o)

    def __state_iter(self, iterable):
        state = Factorial.LOW
        for val in iterable:
            yield (val, state)
            state = Factorial.HIGH if state == Factorial.LOW else Factorial.LOW
            
    def getTableValue(self, table, i, j):
        no_val, nb_val = self.__ff_table[i][j]['no'], self.__ff_table[i][j]['nb']
        no_sol, nb_sol = no_val['sol'], nb_val['sol']
        if table == 'case':
            return '{}/{}'.format(no_sol.case, nb_sol.case)
        elif table == 'profits':
            prof = (nb_sol.profit_man / no_sol.profit_man)*100
            return '{:.2f}\\%'.format(prof)
        elif table == 'retailerprof':
            prof = (nb_sol.profit_ret / no_sol.profit_ret)*100
            return '{:.2f}\\%'.format(prof)
        elif table == 'rho':
            rho_dec = (nb_sol.dec.rho / no_sol.dec.rho) *100
            return '{:.2f}\\%'.format(rho_dec)
        elif table == 'price_new':
            price_dec =  (nb_sol.dec.pn / no_sol.dec.pn) * 100
            return '{:.2f}\\%'.format(price_dec)
        elif table == 'wholesale_price':
            wn_dec = nb_sol.dec.wn
            return '{:.2f}\\%'.format(wn_dec)
        elif table == 'restocking_price':
            rest_price = nb_sol.dec.b
            #'{:.2f}\\%'.format(b_rel)
            return '{:.2f}'.format(rest_price)
        elif table.startswith('par'):
            par = table.split('_')[1]
            if par == 'tau':
                return o_val['par'].tau
            elif par == 'cn':
                return o_val['par'].cn
            elif par == 'cr':
                return o_val['par'].cr
            elif par == 'a':
                return o_val['par'].a
            elif par == 's':
                return o_val['par'].s
            elif par == 'delta':
                return o_val['par'].delta

def main():
    if len(sys.argv) != 2:
        print('usage {} template.tex'.format(os.path.basename(sys.argv[0])))
        sys.exit()

    fullFactorial = Factorial()
    tplfile = sys.argv[1]
    # read template file
    with open(tplfile, 'r') as f:
        template = Template(f.read())

    # open output file
    with open('output.tex', 'w', newline='\n') as f:
        renderedString  = template.render({
            'full': fullFactorial,
            'ft' : fullFactorial.getTableValue,
            'at' : fullFactorial.getAnalysisLine,
            'esc' : escape_tex
        })
        f.write(renderedString)



if __name__ == '__main__':
    main()