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
    
def plot_heatmap_profit_diff():
    _heatmap_profit_diff_data()

def _heatmap_profit_diff_data():
    """ returns a matrice """
    tau, s, cr, delta = 0.3, 0.1, 0.1, 0.4
    step_size_a, lower_bound_a, upper_bound_a = 0.001, .01, .04
    step_size_cn, lower_bound_cn, upper_bound_cn = 0.01, cr, 1
    
    gen = _a_cn_generator(tau, s, cr, delta, step_size_a, step_size_cn, lower_bound_a, upper_bound_a, lower_bound_cn, upper_bound_cn)
    nr_cols = int((upper_bound_a-lower_bound_a)/step_size_a) + 1
    nr_lines = int((upper_bound_cn-lower_bound_cn)/step_size_cn) + 1
    mat = np.empty([nr_lines, nr_cols])
    mat[:,:] = 0
    solver_m1, solver_m2 = ModelOneNumericalSolver(), ModelTwoNumericalSolver()
    i = 0
    for line, col, par_model_1, par_model_2 in gen:
        # calc solutions
        sol_model_1 = solver_m1.optimize(par_model_1)
        sol_model_2 = solver_m2.optimize(par_model_2)
        # both possible?
        if None not in (sol_model_1, sol_model_2):
            mat[line, col] = (sol_model_2.profit_man - sol_model_1.profit_man) / sol_model_1.profit_man
        i += 1
    #mat = mat.transpose()
    print(i)
    plt.imshow(mat, cmap='hot', interpolation='nearest')
    plt.show()
        
def _a_cn_generator(
    tau, s, cr, delta, step_size_a=0.01, step_size_cn=0.1, lower_bound_a = 0.01,
    upper_bound_a=.04, lower_bound_cn=0, upper_bound_cn=1):
    """
        Helper generator for generate parameters for plotting
        it varies cn from [cr, 1] and a from [0, upper_bound_a]
        Returns a tuple of (par_model_1, par_model_2)
    """
    assert lower_bound_cn >= cr
    round_digits_a = int(log(1/step_size_a, 10))
    round_digits_cn = int(log(1/step_size_cn, 10))
    for line, cn in enumerate(drange(lower_bound_cn, upper_bound_cn, step_size_cn)):
        for col, a in enumerate(drange(0.01, upper_bound_a, step_size_a)):
            par_model_2 = Parameter(
                MODEL_2, tau=round(tau, 1), a=round(a, round_digits_a), s=round(s, 1),
                cr=round(cr, 1), cn=round(cn, round_digits_cn), delta=round(delta, 1))
            par_model_1 = Parameter(
                MODEL_1, tau=round(tau, 1), a=round(a, round_digits_a), s=round(s, 1),
                cr=round(cr, 1), cn=round(cn, round_digits_cn), delta=round(delta, 1))
            yield (line, col, par_model_1, par_model_2)

if __name__ == '__main__':
    # parse the command line
    parser = argparse.ArgumentParser(description='Plotting module of clsc solver')
    parser.add_argument('-plot', type=__plot_name, nargs=1, required=True)
    parser.add_argument('-output', type=__output_file, nargs=1, required=False)
    args = parser.parse_args()
    
    if args.plot[0] == PLOT_PROFIT_DIFFERENCE:
        plot_heatmap_profit_diff()