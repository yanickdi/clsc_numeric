import scipy.optimize
import sys

_CASE_ONE, _CASE_TWO = 1, 2

class ModelOneNumericalSolver:
    """
        This class offers methods to solve Model 1 (Without Online Store of the Manufacturer) numerically
    """
    def __init__(self):
        pass
    
    def calc_profits(self, const_args, dec_vars):
        """
            Returns the numeric result of the profit of the manufacturer and the retailer (a tuple containing first manufacturer, second retailer)
            having set all decision variables
        """
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        wn, pn, roh, qn = dec_vars['wn'], dec_vars['pn'], dec_vars['roh'], dec_vars['qn']
        manu_profit = qn * (wn * (1- tau/roh) - cn + (tau/roh) * s)
        retailer_profit = qn * (pn - wn) * (1 - tau/roh) - a * roh
        return manu_profit, retailer_profit
    
    def optimize(self, const_args):
        """
        This is the core method of this class. It will return all four
        decision variables to maximize the retailer's profit (with
        respect to the profit maximization condition of the retailer)
        """
        ## test two cases:
        #       case 1 - roh is >= 1
        #       case 2 - roh is == 1
        
        dec_vars_case_1 = self._optimize_case_one(const_args)
        prof_man_case_1, prof_man_ret_1 = self.calc_profits(const_args, dec_vars_case_1)
        dec_vars_case_2 = self._optimize_case_two(const_args)
        prof_man_case_2, prof_man_ret_2 = self.calc_profits(const_args, dec_vars_case_2)
        
        
        case = None
        
        if dec_vars_case_1['roh'] < 1:
            case = _CASE_TWO
        else:
            # roh is greater than 1, we can decide between both
            # the manufacturer decides:
            case = _CASE_ONE if prof_man_case_1 >= prof_man_case_2 else _CASE_TWO
        
        return dec_vars_case_1 if case == _CASE_ONE else dec_vars_case_2
    
    def _optimize_case_one(self, const_args):
        """
        Returns a dec_vars dict having wn,pn,roh and qn stored
        """
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        
        # defining our helper function:
        def __manufacturer_derivation_case_1(wn):
            return (-1/2) - (cn/2) + wn + a * (tau/a)**(1/2)
        
        # let scipy do the job:
        opt = scipy.optimize.fsolve(__manufacturer_derivation_case_1, x0=0.5, full_output=True)
        wn = opt[0][0]
        pn = (1 + wn) / 2
        roh = (1-wn) * (1/2) * (tau/a)**(1/2)
        return {'wn' : wn, 'pn' : pn, 'roh' : roh, 'qn' : 1 - pn}
        
    def _optimize_case_two(self, const_args):
        """
        Returns a dec_vars dict having wn,pn,roh and qn stored
        """
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        
        #helper function
        def __manufacturer_derivation_case_2(wn):
            return (1/2) * (-1 -cn -2*wn * (-1 + tau) + tau + s*tau)
        
        # hello scipy:
        opt = scipy.optimize.fsolve(__manufacturer_derivation_case_2, x0=0.5, full_output=True)
        wn = opt[0][0]
        pn = (1 + wn) / 2
        roh = 1
        return {'wn' : wn, 'pn' : pn, 'roh' : roh, 'qn' : 1 - pn}
    
def check_args(args):
    tau, a, s, cn = args['tau'], args['a'], args['s'], args['cn']
    assert 0 <= tau <= 1
    assert 0 <= 0.01
    assert s <= cn <= 1
    assert 0 <= s <= cn
    assert s <= cn
    
    
def build_args(tau, a, s, cn):
    args =  {
        'tau' : tau,
        'a'   : a,
        's'   : s,
        'cn'  : cn
    }
    check_args(args)
    return args
    
    
if __name__ == '__main__':
    print('You cannot call this file directly')
    solver = ModelOneNumericalSolver()
    pass