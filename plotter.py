import argparse
from math import log

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from solver import drange, Parameter, MODEL_1, MODEL_2, ModelOneNumericalSolver, ModelTwoNumericalSolver
import solver

PLOT_PROFIT_DIFFERENCE = 'profit-difference'
PLOT_RHO_DIFFERENCE = 'rho-difference'
PLOT_CASES_MODEL_ONE = 'cases-model-one'
PLOT_CASES_MODEL_TWO = 'cases-model-two'

_ALL_CASES_MODEL_1 = [solver._CASE_ONE, solver._CASE_TWO]
_ALL_CASES_MODEL_2 = [solver._CASE_ONE_A, solver._CASE_ONE_B, solver._CASE_ONE_C,
                       solver._CASE_TWO_A, solver._CASE_TWO_B, solver._CASE_TWO_C]

class CountourPlotter:
    def __init__(self, type, params):
        self.type = type
        self.tau, self.s, self.cr, self.delta = params['tau'], params['s'], params['cr'], params['delta']
        self.step_size_a, self.lower_bound_a, self.upper_bound_a = params['step_size_a'], params['lower_bound_a'], params['upper_bound_a']
        self.step_size_cn, self.lower_bound_cn, self.upper_bound_cn = params['step_size_cn'], params['lower_bound_cn'], params['upper_bound_cn']
        self._calc_func = {
            PLOT_RHO_DIFFERENCE : self.__rho_calc_func,
            PLOT_PROFIT_DIFFERENCE : self.__profit_calc_func,
            PLOT_CASES_MODEL_ONE : self.__case_calc_func_model_one,
            PLOT_CASES_MODEL_TWO : self.__case_calc_func_model_two
        }[type]
        self.absolute = params['absolute']
        self.gray = params['gray']
        if type == PLOT_PROFIT_DIFFERENCE and self.absolute: self._calc_func = self.__profit_calc_func_absolute
        elif type == PLOT_RHO_DIFFERENCE  and self.absolute: self._calc_func = self.__rho_calc_func_absolute
        self.matrix = None
    
    def _a_cn_generator(self):
        """
            Helper generator to generate parameters for plotting
            it varies cn from [cr, 1] and a from [0, upper_bound_a]
            Returns a tuple of (par_model_1, par_model_2)
        """
        #assert self.lower_bound_cn >= self.cr
        for line, cn in enumerate(drange(self.lower_bound_cn, self.upper_bound_cn, self.step_size_cn)):
            for col, a in enumerate(drange(self.lower_bound_a, self.upper_bound_a, self.step_size_a)):
                par_model_1 = Parameter(MODEL_1, tau=self.tau, a=a, s=self.s, cn=cn)
                par_model_2 = Parameter(MODEL_2, tau=self.tau, a=a, s=self.s,  cr=self.cr, cn=cn, delta=self.delta)
                yield (line, col, par_model_1, par_model_2)
    
    def calc(self):
        self.nr_cols = int((self.upper_bound_a-self.lower_bound_a)/self.step_size_a) + 1
        self.nr_lines = int((self.upper_bound_cn-self.lower_bound_cn)/self.step_size_cn) + 1
        self.matrix = np.zeros([self.nr_lines, self.nr_cols])
        solver_m1, solver_m2 = ModelOneNumericalSolver(), ModelTwoNumericalSolver()
        i = 0
        for line, col, par_model_1, par_model_2 in self._a_cn_generator():
            # calc solutions
            if par_model_1.a == .0:
                sol_model_1, sol_model_2 = None, None
            else:
                sol_model_1 = solver_m1.optimize(par_model_1)
                sol_model_2 = solver_m2.optimize(par_model_2)
            self.matrix[line, col] = self._calc_func(sol_model_1, sol_model_2, par_model_2)
            i += 1
        print(i)
        
    def __rho_calc_func(self, sol_model_1, sol_model_2, par):
        if sol_model_1 is None or sol_model_2 is None:
            return np.nan
        return (sol_model_2.dec.rho - sol_model_1.dec.rho) / sol_model_1.dec.rho
        
    def __rho_calc_func_absolute(self, sol_model_1, sol_model_2, par):
        if sol_model_1 is None or sol_model_2 is None:
            return np.nan
        return (sol_model_2.dec.rho - sol_model_1.dec.rho)
        
    def __profit_calc_func(self, sol_model_1, sol_model_2, par):
        if sol_model_1 is None or sol_model_2 is None or sol_model_1.profit_man < 0.01:
            return np.nan
        rel =  (sol_model_2.profit_man - sol_model_1.profit_man) / sol_model_1.profit_man
        return rel
        
    def __profit_calc_func_absolute(self, sol_model_1, sol_model_2, par):
        if sol_model_1 is None or sol_model_2 is None:
            return np.nan
        rel =  (sol_model_2.profit_man - sol_model_1.profit_man)
        return rel
        
    def __case_calc_func_model_one(self, sol_model_1, sol_model_2, par):
        if sol_model_1 == None:
            return np.nan
        else:
            return _ALL_CASES_MODEL_1.index(sol_model_1.case) + 1
            
    def __case_calc_func_model_two(self, sol_model_1, sol_model_2, par):
        if sol_model_2 == None:
            return np.nan
        else:
            return _ALL_CASES_MODEL_2.index(sol_model_2.case) + 1
            
    def __plot_cases(self):
        """ not used any more"""
        if self.type == PLOT_CASES_MODEL_ONE:
            title = r'Which case will be active in the event of no Online Shop'
            side_title = r'0 = no solution, 1 = case 1,   2 = case 2'
        elif self.type == PLOT_CASES_MODEL_TWO:
            title = r'Which case will be active in the event of having an Online Shop'
            side_title = r'0=no sol, 1=c1a, 2=c1b, 3=c1c, 4=c2a, 5=c2b, 6=c3b'
            
        fig = plt.figure()
        ax = plt.subplot(111)
        
        imgplot = ax.matshow(self.matrix)
        fig.suptitle(title)
        cbar = plt.colorbar(imgplot)
        # Put a text of paramaters below current axis
        txt = r'par: $\tau$={:.2f}, s={:.2f}, cr={:.2f}, $\delta$={:.2f}'.format(self.tau, self.s, self.cr, self.delta)
        fig.text(0.6, 0.05, txt, fontsize=8, bbox=dict(facecolor='white', alpha=0.5))
        plt.subplots_adjust(bottom=0.15)
        plt.show()
        
    def plot_contourf(self):
        if self.type == PLOT_PROFIT_DIFFERENCE:
            title = 'Increase of Profits with vs. without Online Store'
            side_title = r'relative increase of $\pi_{man}$'
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
            cbar_ticks = range(2)
        elif self.type == PLOT_CASES_MODEL_TWO:
            title = r'Which case will be active in the event of having an Online Shop'
            side_title = r'0=no sol, 1=c1a, 2=c1b, 3=c1c, 4=c2a, 5=c2b, 6=c3b'
            cbar_ticks = range(2,7)
        cmap = cm.gray if self.gray else None            
        fig = plt.figure()
        ax = plt.subplot(111)
        
        x_vector = np.linspace(self.lower_bound_a, self.upper_bound_a, num=self.nr_cols)       # x is a
        y_vector = np.linspace(self.lower_bound_cn, self.upper_bound_cn, num=self.nr_lines)    # y is cn
        cont = ax.contourf(x_vector, y_vector, self.matrix, label='heatmap', cmap=cmap)
        fig.suptitle(title)
        ax.set_xlabel('a')
        ax.set_ylabel(r'$c_n$')
        
        cbar = plt.colorbar(cont, ticks=cbar_ticks)
        cbar.ax.set_ylabel(side_title)
        
        # Put a text of paramaters below current axis
        txt = r'par: $\tau$={:.2f}, s={:.2f}, cr={:.2f}, $\delta$={:.2f}'.format(self.tau, self.s, self.cr, self.delta)
        fig.text(0.6, 0.05, txt, fontsize=8, bbox=dict(facecolor='white', alpha=0.5))
        plt.subplots_adjust(bottom=0.15)
        plt.show()
        
    def plot(self):
        self.plot_contourf()
            
def __parser_output_file(string):
    if string == 'stdout':
        return string
    if len(string.split('.')) <= 1:
        raise argparse.ArgumentTypeError('The output file has to have a suffix.')
    suffix = string.split('.')[-1]
    if suffix not in ('png'):
        raise argparse.ArgumentTypeError('Supported output types are: .png')
    else:
        for all in allowed:
            if all == string: return all
    
def __parser_plot_name(string):
    allowed = (PLOT_PROFIT_DIFFERENCE, PLOT_RHO_DIFFERENCE, PLOT_CASES_MODEL_ONE, PLOT_CASES_MODEL_TWO)
    if string not in allowed:
        raise argparse.ArgumentTypeError('Implemented plots are: ' + ' or '.join(allowed))
    return string

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
    absolute = args.absolute
    gray = args.gray
    
    if quality == 'high':
        step_size_a = .0001
        step_size_cn = .001
    elif quality == 'low':
        step_size_a = .001
        step_size_cn = .01
    
    plotter = CountourPlotter(args.plot[0], params={
        'tau': .3, 's': .07, 'cr': .1, 'delta' : .3,
        'step_size_a' : step_size_a, 'lower_bound_a' : .0, 'upper_bound_a' : .04,
        'step_size_cn' : step_size_cn, 'lower_bound_cn' : .0, 'upper_bound_cn' : .7,
        'absolute' : absolute,
        'gray' : gray
    })
    plotter.calc()
    plotter.plot()
