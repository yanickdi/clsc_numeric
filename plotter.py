import argparse
from math import log

import numpy as np
import matplotlib.pyplot as plt

from solver import drange, Parameter, MODEL_1, MODEL_2, ModelOneNumericalSolver, ModelTwoNumericalSolver

PLOT_PROFIT_DIFFERENCE = 'profit-difference'
PLOT_RHO_DIFFERENCE = 'rho-difference'

def __output_file(string):
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
    
def __plot_name(string):
    allowed = (PLOT_PROFIT_DIFFERENCE, PLOT_RHO_DIFFERENCE)
    if string not in allowed:
        raise argparse.ArgumentTypeError('Implemented plots are: ' + ' or '.join(allowed))
    return string

class CountourPlotter:
    def __init__(self, type, params):
        self.type = type
        assert type in (PLOT_PROFIT_DIFFERENCE, PLOT_RHO_DIFFERENCE)
        self.tau, self.s, self.cr, self.delta = params['tau'], params['s'], params['cr'], params['delta']
        self.step_size_a, self.lower_bound_a, self.upper_bound_a = params['step_size_a'], params['lower_bound_a'], params['upper_bound_a']
        self.step_size_cn, self.lower_bound_cn, self.upper_bound_cn = params['step_size_cn'], params['lower_bound_cn'], params['upper_bound_cn']
        self._calc_func = self.__rho_calc_func if type == PLOT_RHO_DIFFERENCE else self.__profit_calc_func
        self.matrix = None
    
    def _a_cn_generator(self):
        """
            Helper generator to generate parameters for plotting
            it varies cn from [cr, 1] and a from [0, upper_bound_a]
            Returns a tuple of (par_model_1, par_model_2)
        """
        assert self.lower_bound_cn >= self.cr
        for line, cn in enumerate(drange(self.lower_bound_cn, self.upper_bound_cn, self.step_size_cn)):
            for col, a in enumerate(drange(0.01, self.upper_bound_a, self.step_size_a)):
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
            sol_model_1 = solver_m1.optimize(par_model_1)
            sol_model_2 = solver_m2.optimize(par_model_2)
            # both possible?
            if None not in (sol_model_1, sol_model_2):
                self.matrix[line, col] = self._calc_func(sol_model_1, sol_model_2)
            i += 1
        print(i)
        
    def __rho_calc_func(self, sol_model_1, sol_model_2):
        return (sol_model_2.dec.rho - sol_model_1.dec.rho) / sol_model_1.dec.rho
        
    def __profit_calc_func(self, sol_model_1, sol_model_2):
        return (sol_model_2.profit_man - sol_model_1.profit_man) / sol_model_1.profit_man
        
    def plot(self):
        if self.type == PLOT_PROFIT_DIFFERENCE:
            title = 'Increase of Profits with vs. without Online Store'
            side_title = r'increase of $\pi_{man}$ in percentage'
        elif self.type == PLOT_RHO_DIFFERENCE:
            title = r'Increase of Efforts $\rho$ with vs. without Online Store'
            side_title = r'increase of $\rho$ in percentage'
        
        x_vector = np.linspace(self.lower_bound_a, self.upper_bound_a, num=self.nr_cols)       # x is a
        y_vector = np.linspace(self.lower_bound_cn, self.upper_bound_cn, num=self.nr_lines)    # y is cn
        cont = plt.contourf(x_vector, y_vector, self.matrix)
        plt.legend()
        plt.title(title)
        plt.xlabel('a')
        plt.ylabel('cn')
        
        cbar = plt.colorbar(cont)
        cbar.ax.set_ylabel(side_title)
        plt.show()

if __name__ == '__main__':
    # parse the command line
    parser = argparse.ArgumentParser(description='Plotting module of clsc solver')
    parser.add_argument('-plot', type=__plot_name, nargs=1, required=True)
    parser.add_argument('-output', type=__output_file, nargs=1, required=False)
    args = parser.parse_args()
    
    plotter = CountourPlotter(args.plot[0], params={
        'tau': .3, 's': .1, 'cr': .1, 'delta' : .4,
        'step_size_a' : .0001, 'lower_bound_a' : .01, 'upper_bound_a' : .04,
        'step_size_cn' : .001, 'lower_bound_cn' : .1, 'upper_bound_cn' : .4})
    plotter.calc()
    plotter.plot()