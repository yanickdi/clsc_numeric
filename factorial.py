import sys, os.path, re
from statistics import mean

from jinja2 import Template
from solver import SolverProxy, Parameter, MODEL_2, MODEL_1

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
    tau_lh   = lambda : (.15, .5)
    s_lh     = lambda cn: (0, cn/2)
    cr_lh    = lambda cn, delta : (.1*cn, .4*cn)
    delta_lh = lambda : (.5, .85)
    cn_lh    = lambda : (.1, .5)
    a_lh     = lambda : (.001, .01) 

    LOW, HIGH = 'LOW', 'HIGH'

    def __init__(self):
        self.__dict = {}
        self.__proxy = SolverProxy()
        self.__calc__()
        self.__fullFactorialTable()

    def __key(self, model, tau, a, s, cr, cn, delta):
        assert model in (MODEL_1, MODEL_2)
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

        for tau, tau_state in iter(p.tau_lh()):
            for delta, delta_state in iter(p.delta_lh()):
                for cn, cn_state in iter(p.cn_lh()):
                    for cr, cr_state in iter(p.cr_lh(cn, delta)):
                        for s, s_state in iter(p.s_lh(cn)):
                            for a, a_state in iter(p.a_lh()):
                                # calculate model one
                                par_n = Parameter(MODEL_1, tau=tau, a=a, s=s, cn=cn)
                                sol_n = proxy.calculate(par_n)
                                key_n = self.__key(MODEL_1, tau_state, a_state, s_state, cr_state, cn_state, delta_state)
                                self._add(key_n, par_n, sol_n)

                                # calculate model two
                                par_o = Parameter(MODEL_2, tau=tau, a=a, s=s, cr=cr, cn=cn, delta=delta)
                                sol_o = proxy.calculate(par_o)
                                key_o = self.__key(MODEL_2, tau_state, a_state, s_state, cr_state, cn_state, delta_state)
                                self._add(key_o, par_o, sol_o)
                                
    def __fullFactorialTable(self):
        self.__ff_table = [[None for j in range(8)] for i in range(8)]
        lohi = (Factorial.LOW, Factorial.HIGH)
        for i, (s, cr, a) in enumerate([(s, cr, a) for s in lohi for cr in lohi for a in lohi]):
            for j, (cn, delta, tau) in enumerate([(cn, delta, tau) for cn in lohi for delta in lohi for tau in lohi]):
                key_n = self.__key(MODEL_1, tau, a, s, cr, cn, delta)
                key_o = self.__key(MODEL_2, tau, a, s, cr, cn, delta)
                val_n = self._get(key_n)
                val_o = self._get(key_o)
                self.__ff_table[i][j] = {'o' : val_o, 'n' : val_n}
                
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

    def getAnalysisLineAbsoluteMN(self, table, par):
        pars = self.__anal_line_parameters(table, par)
        profits, rhos, pns, wns, prs = [], [], [], [], []
        sol_ids = []
        # calc all vals
        for i, par in enumerate(pars):
            tau, a, s, cr, cn, delta = par
            sol_n, sol_o = self.__both_solutions(tau, a, s, cr, cn, delta)
            if sol_o.case in ('1a', '1b', '1c'):
                sol_ids.append(i)
                profits.append( sol_o.profit_man  )
                rhos.append( sol_o.dec.rho )
                pns.append( sol_o.dec.pn )
                wns.append( sol_o.dec.wn )
                prs.append( sol_o.dec.qn )

        if len(profits) == 0:
            return r'-&-&-&-&-'
        mean_profits = mean(profits)
        mean_rhos = mean(rhos)
        mean_pns = mean(pns)
        mean_wns = mean(wns)
        mean_prs = mean(prs)
        id_string = sum(sol_ids)
        return r'{:.5f}&{:.5f}&{:.5f}&{:.5f}&{:.5f} (n={},id={})'.format(
            mean_profits, mean_rhos, mean_pns, mean_wns, mean_prs, len(profits), id_string)

    def getAnalysisLineAbsoluteMO(self, table, par):
        pars = self.__anal_line_parameters(table, par)
        profits, rhos, pns, wns, prs, qrs, real_prs = [], [], [], [], [], [], []
        sol_ids = []
        # calc all vals
        for i, par in enumerate(pars):
            tau, a, s, cr, cn, delta = par
            sol_n, sol_o = self.__both_solutions(tau, a, s, cr, cn, delta)
            if sol_o.case in ('2a', '2b', '2c'):
            #if sol_o.case in ('1c'):
                sol_ids.append(i)
                profits.append( sol_o.profit_man  )
                rhos.append( sol_o.dec.rho )
                pns.append( sol_o.dec.pn )
                wns.append( sol_o.dec.wn )
                prs.append( sol_o.dec.qn )
                qrs.append (sol_o.dec.qr )
                real_prs.append( sol_o.dec.pr)

        if len(profits) == 0:
            return r'-&-&-&-&-'
        mean_profits = mean(profits)
        mean_rhos = mean(rhos)
        mean_pns = mean(pns)
        mean_wns = mean(wns)
        mean_prs = mean(prs)
        mean_qrs = mean(qrs)
        mean_real_prs = mean(real_prs)
        mean_real_prs = mean(real_prs)
        id_string = sum(sol_ids)
        return r'{:.5f}&{:.5f}&{:.5f}&{:.5f}&{:.5f} qr={:.5f} pr={:.5f}, (n={},id={})'.format(
            mean_profits, mean_rhos, mean_pns, mean_wns, mean_prs, mean_qrs, mean_real_prs, len(profits), id_string)
                
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
        o_val, n_val = self.__ff_table[i][j]['o'], self.__ff_table[i][j]['n']
        o_sol, n_sol = o_val['sol'], n_val['sol']
        
        if table == 'case':
            o_case = '-' if o_sol is None else o_sol.case
            n_case = '-' if n_sol is None else n_sol.case
            return '{}/{}'.format(o_case, n_case)
        
        if None in (o_sol, n_sol): return '-'
        if table == 'profits':
            prof = (o_val['sol'].profit_man / n_val['sol'].profit_man)*100
            return '{:.2f}\\%'.format(prof)
        elif table == 'retailerprof':
            prof = (o_val['sol'].profit_ret / n_val['sol'].profit_ret)*100
            return '{:.2f}\\%'.format(prof)
        elif table == 'rho':
            rho_dec = (o_val['sol'].dec.rho / n_val['sol'].dec.rho) *100
            return '{:.2f}\\%'.format(rho_dec)
        elif table == 'price_new':
            price_dec =  (o_val['sol'].dec.pn / n_val['sol'].dec.pn) * 100
            return '{:.2f}\\%'.format(price_dec)
        elif table == 'wholesale_price':
            wn_dec =  (o_val['sol'].dec.wn / n_val['sol'].dec.wn) * 100
            return '{:.2f}\\%'.format(wn_dec)
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
    output_file = 'output.tex' if tplfile.endswith('.tex') else 'output.csv'
    with open(output_file, 'w', newline='\n') as f:
        renderedString  = template.render({
            'full': fullFactorial,
            'ft' : fullFactorial.getTableValue,
            'at' : fullFactorial.getAnalysisLine,
            'esc' : escape_tex
        })
        f.write(renderedString)



if __name__ == '__main__':
    main()