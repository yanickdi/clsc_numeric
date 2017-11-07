import argparse
from math import log
import sys, os
if __name__ == '__main__': sys.path.append(os.path.abspath('..'))

import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib import cm

from clsc_numeric import solver
from clsc_numeric.solver import drange, Parameter, MODEL_1, MODEL_2, MODEL_NB,  \
    MODEL_1_QUAD, MODEL_2_QUAD, \
    ModelOneNumericalSolver, ModelTwoNumericalSolver, SolverProxy

# my color palette:
RED_DARK = '#a70000'
RED_MEDIUM = '#ec5300'
RED_LIGHT = '#db1414'
BLUE_DARK = '#0d0548'
BLUE_MEDIUM = '#1b51a6'
BLUE_LIGHT = '#1b51a6'

PLOT_PROFIT_DIFFERENCE_MAN = 'profit-difference-man'
PLOT_PROFIT_DIFFERENCE_RET = 'profit-difference-ret'
PLOT_RHO_DIFFERENCE = 'rho-difference'
PLOT_CASES_MODEL_ONE = 'cases-model-one'
PLOT_CASES_MODEL_TWO = 'cases-model-two'
PLOT_COMPARE_CASES = 'compare-cases'
PLOT_FIXED_PLOT = 'fixed-plot'
PLOT_SPONT_PLOT = 'spont-plot'
ALLOWED_PLOTS = PLOT_PROFIT_DIFFERENCE_MAN, PLOT_PROFIT_DIFFERENCE_RET, PLOT_RHO_DIFFERENCE, PLOT_CASES_MODEL_ONE, PLOT_CASES_MODEL_TWO, PLOT_COMPARE_CASES, PLOT_FIXED_PLOT, PLOT_SPONT_PLOT
AUTOMATED_PLOTS = 'automated-plots'

_ALL_CASES_MODEL_1 = [solver._CASE_ONE, solver._CASE_TWO]
_ALL_CASES_MODEL_2 = [solver._CASE_ONE_A, solver._CASE_ONE_B, solver._CASE_ONE_C,
                       solver._CASE_TWO_A, solver._CASE_TWO_B, solver._CASE_TWO_C]

class CountourPlotter:
    def __init__(self, type, params):
        self.proxy = SolverProxy()
        self.type = type
        self.tau, self.s, self.cr, self.delta = params['tau'], params['s'], params['cr'], params['delta']
        self.step_size_a, self.lower_bound_a, self.upper_bound_a = params['step_size_a'], params['lower_bound_a'], params['upper_bound_a']
        self.step_size_cn, self.lower_bound_cn, self.upper_bound_cn = params['step_size_cn'], params['lower_bound_cn'], params['upper_bound_cn']
        self._calc_func = {
            PLOT_RHO_DIFFERENCE : self.__rho_calc_func,
            PLOT_PROFIT_DIFFERENCE_MAN : self.__profit_calc_func_man,
            PLOT_PROFIT_DIFFERENCE_RET : self.__profit_calc_func_ret,
            PLOT_CASES_MODEL_ONE : self.__case_calc_func_model_one,
            PLOT_CASES_MODEL_TWO : self.__case_calc_func_model_two,
            PLOT_COMPARE_CASES : self.__cases_compare_calc_func
        }[type]
        self.absolute = params['absolute']
        self.gray = params['gray']
        self.output = params['output']
        self.nolegend = params['nolegend']
        if type == PLOT_PROFIT_DIFFERENCE_MAN and self.absolute: self._calc_func = self.__profit_calc_func_man_absolute
        elif type == PLOT_PROFIT_DIFFERENCE_RET and self.absolute: self._calc_func = self.__profit_calc_func_ret_absolute
        elif type == PLOT_RHO_DIFFERENCE  and self.absolute: self._calc_func = self.__rho_calc_func_absolute
        self.matrix = None
    
    def _a_cn_generator(self):
        """
            Helper generator to generate parameters for plotting
            it varies cn from [cr, 1] and a from [0, upper_bound_a]
            Returns a tuple of (par_model_1, par_model_2)
        """
        def __cr(cn):
            #if self.cr == 'delta*cn/2':
            #    return self.delta * cn / 2
            if self.cr == '0.4*cn':
                return 0.4 * cn
            else:
                return self.cr
                
        def __s(cn):
            if self.s == 'cn/2':
                return cn / 2
            else:
                return self.s
        
        if self.s == '0.4*cn' and self.cr == '0.4*cn':
            for line, cn in enumerate(drange(self.lower_bound_cn, self.upper_bound_cn, self.step_size_cn)):
                for col, a in enumerate(drange(self.lower_bound_a, self.upper_bound_a, self.step_size_a)):
                    par_model_1 = Parameter(MODEL_1, tau=self.tau, a=a, s=0.4 * cn, cn=cn)
                    par_model_2 = Parameter(MODEL_2, tau=self.tau, a=a, s=0.4 * cn,  cr=0.4 * cn, cn=cn, delta=self.delta)
                    yield (line, col, par_model_1, par_model_2)
        else:
            for line, cn in enumerate(drange(self.lower_bound_cn, self.upper_bound_cn, self.step_size_cn)):
                for col, a in enumerate(drange(self.lower_bound_a, self.upper_bound_a, self.step_size_a)):
                    par_model_1 = Parameter(MODEL_1, tau=self.tau, a=a, s=__s(cn), cn=cn)
                    par_model_2 = Parameter(MODEL_2, tau=self.tau, a=a, s=__s(cn),  cr=__cr(cn), cn=cn, delta=self.delta)
                    yield (line, col, par_model_1, par_model_2)
    
    def calc(self):
        self.nr_cols = int((self.upper_bound_a-self.lower_bound_a)/self.step_size_a) + 1
        self.nr_lines = int((self.upper_bound_cn-self.lower_bound_cn)/self.step_size_cn) + 1
        self.matrix = np.zeros([self.nr_lines, self.nr_cols])
        solver_m1, solver_m2 = ModelOneNumericalSolver(), ModelTwoNumericalSolver()
        i = 0
        numall = self.nr_cols * self.nr_lines
        print(numall)
        #self.proxy.beginWrite()
        print(self.nr_cols*self.nr_lines)
        for line, col, par_model_1, par_model_2 in self._a_cn_generator():
            # calc solutions
            if par_model_1.a == .0:
                sol_model_1, sol_model_2 = None, None
            else:
                #sol_model_1 = solver_m1.optimize(par_model_1)
                #sol_model_2 = solver_m2.optimize(par_model_2)
                sol_model_1 = self.proxy.read_or_calc_and_write(par_model_1)
                sol_model_2 = self.proxy.read_or_calc_and_write(par_model_2)
                if i % 1000 == 0:
                    self.proxy.commit()
                    print(i)
            self.matrix[line, col] = self._calc_func(sol_model_1, sol_model_2, par_model_2)
            i += 1
        #self.proxy.endWrite()
        
    def __rho_calc_func(self, sol_model_1, sol_model_2, par):
        if sol_model_1 is None or sol_model_2 is None:
            return np.nan
        return (sol_model_2.dec.rho - sol_model_1.dec.rho) / sol_model_1.dec.rho
        
    def __rho_calc_func_absolute(self, sol_model_1, sol_model_2, par):
        if sol_model_1 is None or sol_model_2 is None:
            return np.nan
        return (sol_model_2.dec.rho - sol_model_1.dec.rho)
        
        
    def __cases_compare_calc_func(self, sol_model_1, sol_model_2, par):
        #0 equal cases
        #1 model 1 case 1, model 2 case 2
        #2 model 1 case 2, model 2 case 1
        #3 solution in 1, not in 2
        #4 solution in 2, not in 1
        if sol_model_1 is None or sol_model_2 is None:
            if sol_model_1:
                return 3
            elif sol_model_2:
                return 4
            else:
                return np.nan
        case_1, case_2 = sol_model_1.case, sol_model_2.case
        switch = 0
        if case_1 == solver._CASE_ONE and _ALL_CASES_MODEL_2.index(case_2) not in [0, 1, 2]:
            switch = 1
        elif case_1 == solver._CASE_TWO and _ALL_CASES_MODEL_2.index(case_2) not in [3,4,5]:
            switch = 2
        
        return switch
        
    def __profit_calc_func_man(self, sol_model_1, sol_model_2, par):
        if sol_model_1 is None or sol_model_2 is None:
            return np.nan
        rel =  (sol_model_2.profit_man - sol_model_1.profit_man) / sol_model_1.profit_man
        #if rel <= 0 or rel >= .3:
        #if rel < .3:
        #    return np.nan
        #return max(0, min(.5, rel))
        return rel
        
        
    def __profit_calc_func_ret(self, sol_model_1, sol_model_2, par):
        if sol_model_1 is None or sol_model_2 is None:
            return np.nan
        rel =  (sol_model_2.profit_ret - sol_model_1.profit_ret) / sol_model_1.profit_ret
        return rel
        
    def __profit_calc_func_ret_absolute(self, sol_model_1, sol_model_2, par):
        if sol_model_1 is None or sol_model_2 is None:
            return np.nan
        rel =  sol_model_2.profit_ret - sol_model_1.profit_ret
        return rel
        
    def __profit_calc_func_man_absolute(self, sol_model_1, sol_model_2, par):
        if sol_model_1 is None and sol_model_2 is None or par.cn <= par.cr:
            return np.nan
        prof_1 = sol_model_1.profit_man if sol_model_1 is not None else 0
        prof_2 = sol_model_2.profit_man if sol_model_2 is not None else 0
        
        rel =  (prof_2 - prof_1)
        return rel
        
    def __case_calc_func_model_one(self, sol_model_1, sol_model_2, par):
        if sol_model_1 == None:
            return np.nan
        else:
            return _ALL_CASES_MODEL_1.index(sol_model_1.case) + 1
            
    def __case_calc_func_model_two(self, sol_model_1, sol_model_2, par):
        if sol_model_2 == None:
            return np.nan
        #elif round(sol_model_2.profit_man, 1) == 0 or round(sol_model_2.profit_ret, 1) == 0:
        #    return np.nan
        else:
            return _ALL_CASES_MODEL_2.index(sol_model_2.case) + 1
        
    def plot_contourf(self):
        cmap = cm.gray if self.gray else None         
        
        fig = plt.figure()
        ax = plt.subplot(111)
        
        levels, colors, extend = None, None, 'both'
        yticklabels = None
        
        if self.type == PLOT_PROFIT_DIFFERENCE_MAN :
            title = 'Increase of Profits with vs. without Online Store'
            side_title = r'relative increase of $\pi_{man}$'
            if self.absolute:
                side_title = side_title.replace('relative', 'absolute')
            cbar_ticks = None
        elif self.type == PLOT_PROFIT_DIFFERENCE_RET :
            title = 'Increase of Profits with vs. without Online Store'
            side_title = r'relative increase of $\pi_{ret}$'
            if self.absolute:
                side_title = side_title.replace('relative', 'absolute')
            cbar_ticks = None
        elif self.type == PLOT_RHO_DIFFERENCE:
            title = r'Increase of Efforts $\rho$ with vs. without Online Store'
            side_title = r'relative increase of $\rho$'
            if self.absolute:
                side_title = side_title.replace('relative', 'absolute')
            cbar_ticks = None
        elif self.type == PLOT_CASES_MODEL_ONE:
            title = r'Which case will be active in the event of no Online Shop'
            side_title = r'0 = no solution, 1 = case 1,   2 = case 2'
            levels = [0, 1, 2]
            yticklabels = r'case $\rho\geq1$', r'case $\rho = 1$'
        elif self.type == PLOT_CASES_MODEL_TWO:
            title = r'Which case will be active in the event of having an Online Shop'
            #side_title = r'1=c1a, 2=c1b, 3=c1c, 4=c2a, 5=c2b, 6=c3b'
            side_title = ''
            levels = [0,1,2,3,4,5,6]
            colors = ['#ee4b4b', '#e12726', '#960000', '#44db6b', '#0ed040', '#009727']
            yticklabels = '1a', '1b', '1c', '2a', '2b', '2c'
        elif self.type == PLOT_COMPARE_CASES:
            title = r'Comparison of Cases Model 1 vs. 2'
            cbar_ticks = None
            side_title = ''
            yticklabels = ['same case', 'M1C1, M2C2', 'M1C2, M2C1', 'S1N2', 'S2N1']
            levels = [0, 1, 2, 3, 4]
        #fig.suptitle(title)
        x_vector = np.linspace(self.lower_bound_a, self.upper_bound_a, num=self.nr_cols)       # x is a
        y_vector = np.linspace(self.lower_bound_cn, self.upper_bound_cn, num=self.nr_lines)    # y is cn
        
        
        cont = ax.contourf(x_vector, y_vector, self.matrix, label='heatmap',
            levels=levels, colors=colors, origin='lower', extend=extend)
        cont.cmap.set_under('yellow')
        cont.cmap.set_over('cyan')
        cbar = fig.colorbar(cont)
        if yticklabels:
            cbar.ax.set_yticklabels(yticklabels)# vertically oriented colorbar
        cbar.ax.set_ylabel(side_title)
        ax.set_xlabel('a')
        ax.set_ylabel(r'$c_n$')
        
        # Put a text of paramaters below current axis
        if not self.nolegend:
            txt = self.__par_txt()
            fig.text(0.20, 0.83, txt, fontsize=12)#, bbox=dict(facecolor='white', alpha=0.5))
            plt.subplots_adjust(bottom=0.15)
        if self.output == None: plt.show()
        else:
            plt.savefig(self.output)
        plt.close(fig)
            
    def __par_txt(self):
        if self.cr == 'delta*cn/2':
            _cr = r'$\frac{\delta*cn}{2}$'
        elif self.cr == '0.4*cn':
            _cr = '40% of cn'
        else:
            _cr = '{:.2f}'.format(self.cr)
        if self.s == 'cn/2':
            _s = r'$\frac{1}{2}cn$'
        else:
            _s = '{:.2f}'.format(self.s)
        return r'par: $\tau$={:.2f}, s={}, cr={}, $\delta$={:.2f}'.format(self.tau, _s, _cr, self.delta)
        
    def plot(self):
        self.plot_contourf()

class FixedPlot:
    def __init__(self, filename=None):
        self.proxy = SolverProxy()
        self.filename = filename
        self.all_a = [float(a) for a in np.linspace(0.0, 0.01, num=50)]
        self.nr_elements = len(self.all_a)
        self.calc_model_nb()
        self.calc_model_o()
        self.calc_model_n()
        self.calc_model_oq()
        self.calc_model_nq()
        
    def calc_model_oq(self):
        # vectors for 'with online shop quadratic':
        self.oq_profit_man = np.zeros(self.nr_elements)
        self.oq_profit_ret = np.zeros(self.nr_elements)
        self.oq_rho = np.zeros(self.nr_elements)
        self.oq_qn = np.zeros(self.nr_elements)
        self.oq_qr = np.zeros(self.nr_elements)
        self.oq_pn = np.zeros(self.nr_elements)
        self.oq_pr = np.zeros(self.nr_elements)
        self.oq_wn = np.zeros(self.nr_elements)
        self.proxy.beginWrite()
        cn = 0.1
        for i, a in enumerate(self.all_a):
            par_oq = Parameter(MODEL_2_QUAD, tau=.09, cr=0.4*cn, s=0.4*cn, delta=.7956, cn=cn, a=a)
            sol_oq = self.proxy.read_or_calc_and_write(par_oq)
            if sol_oq == None or a == 0:
                self.oq_profit_ret[i] = None
                self.oq_profit_man[i] = None
                self.oq_rho[i] = None
                self.oq_qn[i] = None
                self.oq_qr[i] = None
                self.oq_pn[i] = None
                self.oq_pr[i] = None
                self.oq_wn[i] = None
            else:
                self.oq_profit_ret[i] = sol_oq.profit_ret
                self.oq_profit_man[i] = sol_oq.profit_man
                self.oq_rho[i] = sol_oq.dec.rho
                self.oq_qn[i] = sol_oq.dec.qn
                self.oq_qr[i] = sol_oq.dec.qr
                self.oq_pn[i] = sol_oq.dec.pn
                self.oq_pr[i] = sol_oq.dec.pr
                self.oq_wn[i] = sol_oq.dec.wn
        self.proxy.endWrite()
        self.oq_max_qr = self.oq_qn * (0.09/self.oq_rho)
        
    def calc_model_nq(self):
        # vectors for 'without online shop quadratic':
        self.nq_profit_man = np.zeros(self.nr_elements)
        self.nq_profit_ret = np.zeros(self.nr_elements)
        self.nq_rho = np.zeros(self.nr_elements)
        self.nq_qn = np.zeros(self.nr_elements)
        self.nq_qr = np.zeros(self.nr_elements)
        self.nq_pn = np.zeros(self.nr_elements)
        self.nq_wn = np.zeros(self.nr_elements)
        cn = 0.1
        self.proxy.beginWrite()
        for i, a in enumerate(self.all_a):
            par_nq = Parameter(MODEL_1_QUAD, tau=.09, a=a, s=0.4*cn, cn=cn)
            sol_nq = self.proxy.read_or_calc_and_write(par_nq)
            if sol_nq == None or a == 0:
                self.nq_profit_ret[i] = None
                self.nq_profit_man[i] = None
                self.nq_rho[i] = None
                self.nq_qn[i] = None
                self.nq_pn[i] = None
                self.nq_wn[i] = None
            else:
                self.nq_profit_ret[i] = sol_nq.profit_ret
                self.nq_profit_man[i] = sol_nq.profit_man
                self.nq_rho[i] = sol_nq.dec.rho
                self.nq_qn[i] = sol_nq.dec.qn
                self.nq_pn[i] = sol_nq.dec.pn
                self.nq_wn[i] = sol_nq.dec.wn
        self.proxy.endWrite()
        
    def calc_model_o(self):
        # vectors for 'with online shop':
        self.on_profit_man = np.zeros(self.nr_elements)
        self.on_profit_ret = np.zeros(self.nr_elements)
        self.on_rho = np.zeros(self.nr_elements)
        self.on_qn = np.zeros(self.nr_elements)
        self.on_qr = np.zeros(self.nr_elements)
        self.on_pn = np.zeros(self.nr_elements)
        self.on_pr = np.zeros(self.nr_elements)
        self.on_wn = np.zeros(self.nr_elements)
        self.proxy.beginWrite()
        cn = 0.1
        for i, a in enumerate(self.all_a):
            par_o = Parameter(MODEL_2, tau=.09, cr=0.4*cn, s=0.4*cn, delta=.7956, cn=cn, a=a)
            sol_o = self.proxy.read_or_calc_and_write(par_o)
            if sol_o == None or a == 0:
                self.on_profit_ret[i] = None
                self.on_profit_man[i] = None
                self.on_rho[i] = None
                self.on_qn[i] = None
                self.on_qr[i] = None
                self.on_pn[i] = None
                self.on_pr[i] = None
                self.on_wn[i] = None
            else:
                self.on_profit_ret[i] = sol_o.profit_ret
                self.on_profit_man[i] = sol_o.profit_man
                self.on_rho[i] = sol_o.dec.rho
                self.on_qn[i] = sol_o.dec.qn
                self.on_qr[i] = sol_o.dec.qr
                self.on_pn[i] = sol_o.dec.pn
                self.on_pr[i] = sol_o.dec.pr
                self.on_wn[i] = sol_o.dec.wn
        self.proxy.endWrite()
        self.on_max_qr = self.on_qn * (0.09/self.on_rho)

    def calc_model_n(self):
        # vectors for 'without online shop':
        self.no_profit_man = np.zeros(self.nr_elements)
        self.no_profit_ret = np.zeros(self.nr_elements)
        self.no_rho = np.zeros(self.nr_elements)
        self.no_qn = np.zeros(self.nr_elements)
        self.no_qr = np.zeros(self.nr_elements)
        self.no_pn = np.zeros(self.nr_elements)
        self.no_wn = np.zeros(self.nr_elements)
        cn = 0.1
        self.proxy.beginWrite()
        for i, a in enumerate(self.all_a):
            par_n = Parameter(MODEL_1, tau=.09, a=a, s=0.4*cn, cn=cn)
            sol_n = self.proxy.read_or_calc_and_write(par_n)
            if sol_n == None or a == 0:
                self.no_profit_ret[i] = None
                self.no_profit_man[i] = None
                self.no_rho[i] = None
                self.no_qn[i] = None
                self.no_pn[i] = None
                self.no_wn[i] = None
            else:
                self.no_profit_ret[i] = sol_n.profit_ret
                self.no_profit_man[i] = sol_n.profit_man
                self.no_rho[i] = sol_n.dec.rho
                self.no_qn[i] = sol_n.dec.qn
                self.no_pn[i] = sol_n.dec.pn
                self.no_wn[i] = sol_n.dec.wn
        self.proxy.endWrite()
    
    def calc_model_nb(self):
        # vectors for 'without online shop':
        self.nb_profit_man = np.zeros(self.nr_elements)
        self.nb_profit_ret = np.zeros(self.nr_elements)
        self.nb_rho = np.zeros(self.nr_elements)
        self.nb_qn = np.zeros(self.nr_elements)
        self.nb_b = np.zeros(self.nr_elements)
        self.nb_pn = np.zeros(self.nr_elements)
        self.nb_wn = np.zeros(self.nr_elements)
        cn = 0.1
        self.proxy.beginWrite()
        for i, a in enumerate(self.all_a):
            par_nb = Parameter(MODEL_NB, tau=.09, a=a, s=0.4*cn, cn=cn)
            sol_nb = self.proxy.read_or_calc_and_write(par_nb)
            if sol_nb == None or a == 0:
                self.nb_profit_ret[i] = None
                self.nb_profit_man[i] = None
                self.nb_rho[i] = None
                self.nb_qn[i] = None
                self.nb_b[i] = None
                self.nb_pn[i] = None
                self.nb_wn[i] = None
            else:
                self.nb_profit_ret[i] = sol_nb.profit_ret
                self.nb_profit_man[i] = sol_nb.profit_man
                self.nb_rho[i] = sol_nb.dec.rho
                self.nb_b[i] = sol_nb.dec.b
                self.nb_qn[i] = sol_nb.dec.qn
                self.nb_pn[i] = sol_nb.dec.pn
                self.nb_wn[i] = sol_nb.dec.wn
        self.proxy.endWrite()
        

    def plot(self):
        ## ORIGINAL MODELS:
        #self.plot_profits(relative=True)
        self.plot_rhos()
        #self.plot_quantities()
        #self.plot_prices()
        
        ## MODEL NB VS O:
        #self.plot_profits_nb_vs_o()
        #self.plot_rhos_nb_vs_o()
        #self.plot_prices_nb_vs_o()
        
        ## QUADRATIC MODEL:
        #self.plot_profits_nq_vs_oq()
        #self.plot_rhos_nq_vs_oq()
        #self.plot_prices_nq_vs_oq()
        #self.plot_quantities_nq_vs_oq()
        pass

    def plot_profits_nb_vs_o(self):
        fig, ax = plt.subplots()
        
        o_profit_man, o_profit_ret, nb_profit_ret, nb_profit_man = self.on_profit_man, self.on_profit_ret, self.nb_profit_ret, self.nb_profit_man
        
        # with online store:
        pl1,  = ax.plot(self.all_a, o_profit_ret, color=RED_DARK)
        ax.text(self.all_a[-1], o_profit_ret[-1]*1.2, r'$\pi_{R}^{O}$', color=RED_DARK)
        pl2, = ax.plot(self.all_a, o_profit_man, color=RED_MEDIUM)
        ax.text(self.all_a[-1], o_profit_man[-1], r'$\pi_{M}^{O}$', color=RED_MEDIUM)
        
        # with model nb:
        ax.plot(self.all_a, nb_profit_ret, color=BLUE_DARK)
        ax.text(self.all_a[-1], nb_profit_ret[-1], r'$\pi_{R}^{NB}$', color=BLUE_DARK)
        ax.plot(self.all_a, nb_profit_man, color=BLUE_MEDIUM)
        ax.text(self.all_a[-1], nb_profit_man[-1], r'$\pi_{M}^{NB}$', color=BLUE_MEDIUM)
        
        # with model n:
        ax.plot(self.all_a, self.no_profit_ret, color=BLUE_DARK, linestyle='dashed')
        ax.text(self.all_a[-1]*.95, self.no_profit_ret[-1]+.003, r'$\pi_{R}^{N}$', color=BLUE_DARK)
        ax.plot(self.all_a, self.no_profit_man, color=BLUE_MEDIUM, linestyle='dashed')
        ax.text(self.all_a[-1]*.95, self.no_profit_man[-1]-.01, r'$\pi_{M}^{N}$', color=BLUE_MEDIUM)
        ax.set_xlabel('a')
        plt.show()

    def plot_rhos_nb_vs_o(self):
        fig, ax = plt.subplots()
        
        # with online store:
        ax.plot(self.all_a, self.on_rho, color=RED_DARK, label=r'$\rho_{*}^{O}$')
        #ax.text(self.all_a[-1], self.on_rho[-1]*1.2, r'$\rho_{*}^{O}$', color=RED_DARK)
        
        # with model nb:
        ax.plot(self.all_a, self.nb_rho, color=BLUE_DARK, label=r'$\rho_{*}^{NB}$')
        #ax.text(self.all_a[-1], self.nb_rho[-1], r'$\rho_{*}^{NB}$', color=BLUE_DARK)
        
        # with model n:
        ax.plot(self.all_a, self.no_rho, color=BLUE_DARK, linestyle='dashed', label=r'$\rho_{*}^{N}$')
        #ax.text(self.all_a[-1]*.95, self.no_rho[-1]+.003, r'$\rho_{*}^{N}$', color=BLUE_DARK)
        ax.set_xlabel('a')
        ax.set_ylabel(r'Effort ($\rho$)')
        ax.legend()
        plt.show()
        
    def plot_prices_nb_vs_o(self):
        fig, ax = plt.subplots()
        
        # with online store:
        ax.plot(self.all_a, self.on_pn, color=RED_DARK, label=r'$pn_{*}^{O}$')
        ax.text(self.all_a[-1], self.on_pn[-1], r'$pn_{*}^{O}$', color=RED_DARK)
        
        ax.plot(self.all_a, self.on_wn, color=RED_LIGHT, label=r'$wn_{*}^{O}$')
        ax.text(self.all_a[-1], self.on_wn[-1], r'$wn_{*}^{O}$', color=RED_LIGHT)
        
        
        # with model nb:
        ax.plot(self.all_a, self.nb_pn, color=BLUE_DARK, label=r'$pn_{*}^{NB}$')
        ax.text(self.all_a[-1], self.nb_pn[-1], r'$pn_{*}^{NB}$', color=BLUE_DARK)
        
        ax.plot(self.all_a, self.nb_b, color=BLUE_MEDIUM, label=r'$b_{*}^{NB}$')
        ax.text(self.all_a[-1], self.nb_b[-1], r'$b_{*}^{NB}$', color=BLUE_MEDIUM)
        
        ax.plot(self.all_a, self.nb_wn, color=BLUE_LIGHT, label=r'$wn_{*}^{NB}$')
        ax.text(self.all_a[-1], self.nb_wn[-1], r'$wn_{*}^{NB}$', color=BLUE_LIGHT)
        
        # with model n:
        #ax.plot(self.all_a, self.no_rho, color=BLUE_DARK, linestyle='dashed', label=r'$\rho_{*}^{N}$')
        #ax.text(self.all_a[-1]*.95, self.no_rho[-1]+.003, r'$\rho_{*}^{N}$', color=BLUE_DARK)
        ax.set_xlabel('a')
        ax.set_xlim([0, 0.011])
        #ax.legend()
        plt.show()
        
    def plot_profits(self, relative=False):
        fig, ax = plt.subplots()
        on_profit_sc = self.on_profit_man + self.on_profit_ret
        no_profit_sc = self.no_profit_man + self.no_profit_ret
        
        if relative:
            no_profit_man = self.no_profit_man * 0 + 1# / self.no_profit_man
            no_profit_ret = self.no_profit_ret / self.no_profit_man
            on_profit_ret = self.on_profit_ret / self.no_profit_man
            on_profit_man = self.on_profit_man / self.no_profit_man
        else:
            no_profit_man, no_profit_ret, on_profit_ret, on_profit_man = self.no_profit_man, self.no_profit_ret, self.on_profit_ret, self.on_profit_man
        
        # with online store:
        ax.plot(self.all_a, on_profit_ret, color=RED_DARK)
        ax.text(self.all_a[-1], on_profit_ret[-1]*1.2, r'$\pi_{R}^{O}$', color=RED_DARK)
        ax.plot(self.all_a, on_profit_man, color=RED_LIGHT)
        ax.text(self.all_a[-1], on_profit_man[-1], r'$\pi_{M}^{O}$', color=RED_LIGHT)
        
        # without online store:
        ax.plot(self.all_a, no_profit_ret, color=BLUE_DARK)
        ax.text(self.all_a[-1], no_profit_ret[-1], r'$\pi_{R}^{N}$', color=BLUE_DARK)
        ax.plot(self.all_a, no_profit_man, color=BLUE_LIGHT)
        ax.text(self.all_a[-1], no_profit_man[-1], r'$\pi_{M}^{N}$', color=BLUE_LIGHT)
        #pl6, = ax.plot(self.all_a, no_profit_sc, color='#a70000')
        #ax.text(self.all_a[-1], no_profit_sc[-1], r'$\pi_{SC}^{N}$', color=pl6.get_c())
        ax.set_xlabel('a')
        plt.show()
        
    def plot_prices(self):
        fig, ax = plt.subplots()
        # with online store:
        ax.plot(self.all_a, self.on_pn, color=BLUE_MEDIUM)
        ax.text(self.all_a[-1]*.95, self.on_pn[-1]+.01, r'$pn_{O}^{*}$', color=BLUE_MEDIUM)
        ax.plot(self.all_a, self.on_pr, color=BLUE_LIGHT)
        ax.text(self.all_a[-1]*.95, self.on_pr[-1]-.015, r'$pr_{O}^{*}$', color=BLUE_LIGHT)
        # without online store:
        ax.plot(self.all_a, self.no_pn, color=RED_MEDIUM)
        ax.text(self.all_a[-1]*.95, self.no_pn[-1]-.015, r'$pn_{N}^{*}$', color=RED_MEDIUM)
        #plot wholesale prices:
        ax.plot(self.all_a, self.on_wn, color=RED_DARK)
        ax.text(self.all_a[-1]*.95, self.on_wn[-1]+.01, r'$wn_{O}^{*}$', color=RED_DARK)
        ax.plot(self.all_a, self.no_wn, color=BLUE_DARK)
        ax.text(self.all_a[-1]*.95, self.no_wn[-1]+.015, r'$wn_{N}^{*}$', color=BLUE_DARK)
        ax.set_xlabel('a')
        plt.show()
        
    def plot_quantities(self):
        fig, ax = plt.subplots()
        # with online store:
        pl1,  = ax.plot(self.all_a, self.on_qn, color='#1b51a6')
        ax.text(self.all_a[-1]*.95, self.on_qn[-1]+.01, r'$qn_{O}^{*}$', color=pl1.get_c())
        pl2, = ax.plot(self.all_a, self.on_qr, color='#1b51a6')
        ax.text(self.all_a[-1]*.95, self.on_qr[-1]+.01, r'$qr_{O}^{*}$', color=pl2.get_c())
        ax.plot(self.all_a, self.on_max_qr, linestyle='dashed')
        
        # without online store:
        pl3, = ax.plot(self.all_a, self.no_qn, color='#ec5300')
        ax.text(self.all_a[-1]*.95, self.no_qn[-1]+.01, r'$qn_{N}^{*}$', color=pl3.get_c())
        ax.set_xlabel('a')
        plt.show()
        
    def plot_rhos(self):
        fig, ax = plt.subplots()
        pl1,  = ax.plot(self.all_a, self.on_rho, color='#1b51a6')
        ax.text(self.all_a[-1], self.on_rho[-1]*.7,  r'$\rho_{O}^{*}$', color=pl1.get_c())
        pl2, = ax.plot(self.all_a, self.no_rho, color='#a70000')
        ax.text(self.all_a[-1], self.no_rho[-1]*1.3, r'$\rho_{N}^{*}$', color=pl2.get_c())
        ax.set_ylim([0, 5])
        ax.set_xlabel('a')
        ax.set_ylabel(r'Effort ($\rho$)')
        plt.show()
        
    def plot_rhos_nq_vs_oq(self):
        """ no onlineshop quadratic vs. onlineshop quadratic """
        fig, ax = plt.subplots()
        pl1,  = ax.plot(self.all_a, self.oq_rho, color='#1b51a6')
        ax.text(self.all_a[-1], self.oq_rho[-1]*.7,  r'$\rho_{OQ}^{*}$', color=pl1.get_c())
        
        pl2, = ax.plot(self.all_a, self.nq_rho, color='#a70000')
        ax.text(self.all_a[-1], self.nq_rho[-1]*1.3, r'$\rho_{NQ}^{*}$', color=pl2.get_c())
        ax.set_ylim([0, 5])
        ax.set_xlabel('a')
        ax.set_ylabel(r'Effort ($\rho$)')
        plt.show()
        
    def plot_prices_nq_vs_oq(self):
        fig, ax = plt.subplots()
        # with online store:
        ax.plot(self.all_a, self.oq_pn, color=BLUE_MEDIUM)
        ax.text(self.all_a[-1]*.95, self.on_pn[-1]+.01, r'$pn_{OQ}^{*}$', color=BLUE_MEDIUM)
        ax.plot(self.all_a, self.oq_pr, color=BLUE_LIGHT)
        ax.text(self.all_a[-1]*.95, self.on_pr[-1]-.015, r'$pr_{OQ}^{*}$', color=BLUE_LIGHT)
        # without online store:
        ax.plot(self.all_a, self.nq_pn, color=RED_MEDIUM)
        ax.text(self.all_a[-1]*.95, self.no_pn[-1]-.015, r'$pn_{NQ}^{*}$', color=RED_MEDIUM)
        #plot wholesale prices:
        ax.plot(self.all_a, self.oq_wn, color=RED_DARK)
        ax.text(self.all_a[-1]*.95, self.on_wn[-1]+.01, r'$wn_{OQ}^{*}$', color=RED_DARK)
        ax.plot(self.all_a, self.nq_wn, color=BLUE_DARK)
        ax.text(self.all_a[-1]*.95, self.no_wn[-1]+.015, r'$wn_{NQ}^{*}$', color=BLUE_DARK)
        ax.set_xlabel('a')
        plt.show()
        
    def plot_quantities_nq_vs_oq(self):
        fig, ax = plt.subplots()
        # with online store quadratic:
        ax.plot(self.all_a, self.oq_qn, color=RED_DARK)
        ax.text(self.all_a[-1]*.95, self.oq_qn[-1]+.01, r'$qn_{OQ}^{*}$', color=RED_DARK)
        
        ax.plot(self.all_a, self.oq_qr, color=RED_LIGHT)
        ax.text(self.all_a[-1]*.95, self.oq_qr[-1]+.01, r'$qr_{OQ}^{*}$', color=RED_LIGHT)
        
        ax.plot(self.all_a, self.oq_max_qr, linestyle='', marker='x')
        
        
        # without online store quadratic:
        ax.plot(self.all_a, self.no_qn, color=BLUE_DARK)
        ax.text(self.all_a[-1]*.95, self.no_qn[-1]+.01, r'$qn_{NQ}^{*}$', color=BLUE_DARK)
        ax.set_xlabel('a')
        plt.show()
        
    def plot_profits_nq_vs_oq(self, relative=False):
        fig, ax = plt.subplots()
        
        if relative:
            nq_profit_man = self.nq_profit_man * 0 + 1# / self.no_profit_man
            nq_profit_ret = self.nq_profit_ret / self.nq_profit_man
            oq_profit_ret = self.oq_profit_ret / self.nq_profit_man
            oq_profit_man = self.oq_profit_man / self.nq_profit_man
        else:
            nq_profit_man, nq_profit_ret, oq_profit_ret, oq_profit_man = self.nq_profit_man, self.nq_profit_ret, self.oq_profit_ret, self.oq_profit_man
        
        # with online store:
        ax.plot(self.all_a, oq_profit_ret, color=RED_DARK)
        ax.text(self.all_a[-1], oq_profit_ret[-1]*1.2, r'$\pi_{R}^{OQ}$', color=RED_DARK)
        ax.plot(self.all_a, oq_profit_man, color=RED_LIGHT)
        ax.text(self.all_a[-1], oq_profit_man[-1], r'$\pi_{M}^{OQ}$', color=RED_LIGHT)
        
        # without online store:
        ax.plot(self.all_a, nq_profit_ret, color=BLUE_DARK)
        ax.text(self.all_a[-1], nq_profit_ret[-1], r'$\pi_{R}^{N}$', color=BLUE_DARK)
        ax.plot(self.all_a, nq_profit_man, color=BLUE_LIGHT)
        ax.text(self.all_a[-1], nq_profit_man[-1], r'$\pi_{M}^{N}$', color=BLUE_LIGHT)
        
        ax.set_xlabel('a')
        plt.show()
        
class SpontPlot:
    def plot(self):
        print('hello world')
    
    def plot_model_nb_vs_o_profits(self):
        pass

def __parser_output_file(string):
    if string == 'stdout':
        return string
    if len(string.split('.')) <= 1:
        raise argparse.ArgumentTypeError('The output file has to have a suffix.')
    allowed_suffixes = ('png', 'pdf')
    suffix = string.split('.')[-1]
    if suffix not in allowed_suffixes:
        raise argparse.ArgumentTypeError('Supported output types are: .png or .pdf')
    else:
        for all in allowed_suffixes:
            if all == string: return all
    return string
    
def __parser_plot_name(string):
    allowed = ALLOWED_PLOTS
    if string not in allowed:
        raise argparse.ArgumentTypeError('Implemented plots are: ' + ' or '.join(allowed))
    return string
    

def __gen_parameter_automated_plots():
    def __gen():
        for plot in ('profit-difference', 'rho-difference', 'cases-model-one', 'cases-model-two'):
            for absolute in (True, False):
                for delta in (.4, .8):
                    for s in (0, 'cn/2'):
                        for tau in (.02, .05):
                            for cr in (.1, 'delta*cn/2'):
                                yield {'plot' : plot, 'absolute' : absolute, 'delta' : delta, 's' : s, 'tau' : tau, 'cr' : cr}
    for elem in __gen():
        if elem['plot'] in ('cases-model-one', 'cases-model-two') and elem['absolute'] == True:
            pass
        else:
            yield elem

def automated_plots(step_size_a, step_size_cn, gray):
    for nr, parms in enumerate(__gen_parameter_automated_plots()):
        filename = _automated_plots_filename(nr, parms)
        plotter = CountourPlotter(parms['plot'], params={
            'tau': parms['tau'], 's': parms['s'], 'cr': parms['cr'], 'delta' : parms['delta'],
            'step_size_a' : step_size_a, 'lower_bound_a' : .0, 'upper_bound_a' : .04,
            'step_size_cn' : step_size_cn, 'lower_bound_cn' : .0, 'upper_bound_cn' : 1.0,
            'absolute' : parms['absolute'],
            'gray'   : gray,
            'output' : filename
        })
        print('plotting {} ... '.format(filename), end='')
        sys.stdout.flush()
        plotter.calc()
        plotter.plot()
        print('done')
        sys.stdout.flush()
        
def _automated_plots_filename(nr, parms):
    plot_name = parms['plot']
    if plot_name in ('profit-difference', 'rho-difference'):
        plot_name += '_abs_' if parms['absolute'] else '_rel_'
        
    par_keys = ['delta', 's', 'tau', 'cr']
    par_strings = []
    for key in par_keys:
        if type(parms[key]) == str:
            str_expr = parms[key].replace('/', '_over_').replace('*', '_times_')
            par_strings.append('{}_{}'.format(key, str_expr))
        else:
            par_strings.append('{}_{:.3f}'.format(key, parms[key])) #par should be a float
    par_string = '_'.join(par_strings)
    filename = 'plot_{}_{}_{}.pdf'.format(nr+1, plot_name, par_string)
    return filename
    
    
if __name__ == '__main__':
    # parse the command line
    parser = argparse.ArgumentParser(description='Plotting module of clsc solver')
    parser.add_argument('-plot', type=__parser_plot_name, nargs=1, required=True)
    parser.add_argument('-output', type=__parser_output_file, nargs=1, required=False)
    parser.add_argument('--absolute', action='store_true')
    parser.add_argument('--low-qual', action='store_true')
    parser.add_argument('--gray', action='store_true')
    args = parser.parse_args()
    
    quality = 'low' if args.low_qual else 'high'
    if quality == 'high':
        step_size_a = .0001
        step_size_cn = .001
    elif quality == 'low':
        step_size_a = .001
        step_size_cn = .01
    gray = args.gray
    
    plot = args.plot[0]
    if plot == AUTOMATED_PLOTS:
        automated_plots(step_size_a, step_size_cn, gray)
    elif plot == PLOT_FIXED_PLOT:
        fixedPlot = FixedPlot()
        fixedPlot.plot()
    elif plot == PLOT_SPONT_PLOT:
        spontPlot = SpontPlot()
        spontPlot.plot()
    else:
        absolute = args.absolute
        if args.output:
            output = args.output[0]
        else:
            output = None
         
        plotter = CountourPlotter(args.plot[0], params={
            'tau': .09, 's': '0.4*cn', 'cr': '0.4*cn', 'delta' : .7956,
            'step_size_a' : step_size_a, 'lower_bound_a' : .0, 'upper_bound_a' : .025,
            'step_size_cn' : step_size_cn, 'lower_bound_cn' : .0, 'upper_bound_cn' : .9,
            'absolute' : absolute,
            'gray'   : gray,
            'nolegend': True,
            'output' : output
        })
        plotter.calc()
        plotter.plot()
