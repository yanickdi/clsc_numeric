import sys
from collections import namedtuple

class HillClimber:
    """
        Finds an local optimum
        
        The function to optimize `func` must accept 2 arguments:
        (x, args) - where args is passed through and x is a vector including the actual variables
        x and start_vec have the same length
        
        The function must return a floating point number
    """
    
    resObj = namedtuple('HillClimberResult', ['x', 'func'])

    @staticmethod
    def minimize(func, args, x0, iterations=200, step_sizes=.1):
        x, xsize = x0[:], len(x0)
        step_size = [step_sizes for i in range(xsize)]
        acceleration = 1.2
        directions = [-1, 0, 1]
        act_it = 0
        while act_it < iterations:
            before = func(x, args)
            # each variable once:
            for i in range(xsize):
                best_direction = -1
                best_direction_fval = sys.maxsize-10
                # try each direction
                for dir_j in range(len(directions)):
                    # update current point
                    x[i] += step_size[i] * directions[dir_j]
                    # evaluate
                    tmp = func(x, args)
                    # undo update current point
                    x[i] -= step_size[i] * directions[dir_j]
                    # direction was good?
                    if tmp < best_direction_fval:
                        best_direction_fval = tmp
                        best_direction = dir_j
                if directions[best_direction] == 0:
                    # current point was already best, decrease step_size
                    step_size[i] = step_size[i] / acceleration
                    pass
                else:
                    # found a better point:
                    x[i] += step_size[i] * directions[best_direction]
                    # increase step_size
                    step_size[i] = step_size[i] * acceleration
            if func(x, args) == before:
                # we didnt improve, we are (i) exact at optimum (ii) have to zoom in, i.e.
                # decrease *all* step size s
                #for i in range(len(step_size)): step_size[i] = step_size[i] / acceleration
                pass
            else:
                pass
            act_it += 1
            assert func(x, args) != sys.maxsize
        return HillClimber.resObj(x, 0)
        #return {'x': currentP, 'fun': best_score}
        
if __name__ == '__main__':
    from model2 import ModelTwoSolver
    from solver import SolverProxy, Parameter, MODEL_2

    proxy = SolverProxy()
    minimize_func = ModelTwoSolver._minimize_func

    par = Parameter(MODEL_2, tau=0.09, a=0.00040816326530612246, s=0.04000000000000001, cn=0.1, cr=0.04000000000000001, delta=0.7956)

    analySol = proxy.calculate(par)
    hcSol = ModelTwoSolver.solve(par)

    print('anal sol', analySol, analySol.profit_man, analySol.case)
    print('found sol', hcSol, hcSol.profit_man, hcSol.case)

    #res = minimize(func, args=args, x0=[1], method='Nelder-Mead')
    #myres = HillClimber.minimimize(func, args=None, x0=[5,5])
    #print(myres)