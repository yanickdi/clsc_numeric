import argparse
from math import log
import sys, os
if __name__ == '__main__': sys.path.append(os.path.abspath('..'))

import numpy as np
import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib import cm

import solver
from solver import drange, Parameter, MODEL_1, MODEL_2, MODEL_NB,  \
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
        self.count_a, self.lower_bound_a, self.upper_bound_a = params['count_a'], params['lower_bound_a'], params['upper_bound_a']
        self.count_cn, self.lower_bound_cn, self.upper_bound_cn = params['count_cn'], params['lower_bound_cn'], params['upper_bound_cn']
        self._calc_func = {
            PLOT_CASES_MODEL_TWO : self.__case_calc_func_model_two
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
            if self.cr == '0.5*cn':
                return 0.5 * cn
            else:
                return self.cr
                
        def __s(cn):
            if self.s == 'cn/2':
                return cn / 2
            else:
                return self.s
        
        for line, cn in enumerate(np.linspace(self.lower_bound_cn, self.upper_bound_cn, self.count_cn)):
            for col, a in enumerate(np.linspace(self.lower_bound_a, self.upper_bound_a, self.count_a)):
                cn, a = float(cn), float(a)
                par_model_1 = None
                par_model_2 = Parameter(MODEL_2_QUAD, tau=self.tau, a=a, s=__s(cn),  cr=__cr(cn), cn=cn, delta=self.delta)
                yield (line, col, par_model_1, par_model_2)
    
    def calc(self):
        self.nr_cols = self.count_a
        self.nr_lines = self.count_cn
        self.matrix = np.zeros([self.nr_lines, self.nr_cols])
        i = 0
        numall = self.nr_cols * self.nr_lines
        print(numall)
        #self.proxy.beginWrite()
        for line, col, par_model_1, par_model_2 in self._a_cn_generator():
            # calc solutions
            if par_model_2.a == .0:
                sol_model_1, sol_model_2 = None, None
            else:
                #sol_model_1 = solver_m1.optimize(par_model_1)
                #sol_model_2 = solver_m2.optimize(par_model_2)
                #sol_model_1 = self.proxy.read_or_calc_and_write(par_model_1, resolution='low')
                sol_model_1 = None #todo hier fuer model 1 aufschalten
                #sol_model_2 = self.proxy.read_or_calc_and_write(par_model_2, resolution='low')
                sol_model_2 = self.proxy.calculate(par_model_2, resolution=RESOLUTION)
            self.matrix[line, col] = self._calc_func(sol_model_1, sol_model_2, par_model_2)
            if i % 100 == 0:
                #self.proxy.commit()
                print('{} ({:.2f} %)    at {}'.format(i, (i/numall)*100, str(datetime.datetime.now())))
            i += 1
        self.proxy.commit()
        #self.proxy.endWrite()
        
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
        return NotImplementedError()
            
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
        lower_bound_a = .0
        upper_bound_a = .025
        count_x = 400
        lower_bound_cn = .0
        upper_bound_cn = .9
        count_y = 400
        step_size_cn = (upper_bound_cn - lower_bound_cn) / count_y
        RESOLUTION = 'super high'
    elif quality == 'low':
        lower_bound_a = .0
        upper_bound_a = .025
        count_x = 40
        lower_bound_cn = .0
        upper_bound_cn = .9
        count_y = 40
        RESOLUTION = 'super high'
    gray = args.gray
    
    plot = args.plot[0]
    plotter = CountourPlotter(args.plot[0], params={
        'tau': .3, 's': 0.07, 'cr': 0.1, 'delta' : 0.3,
        'count_a' : count_x, 'lower_bound_a' : lower_bound_a, 'upper_bound_a' : upper_bound_a,
        'count_cn' : count_y, 'lower_bound_cn' : lower_bound_cn, 'upper_bound_cn' : upper_bound_cn,
        'absolute' : False,
        'gray'   : gray,
        'nolegend': True,
        'output' : None
    })
    plotter.calc()
    plotter.plot()