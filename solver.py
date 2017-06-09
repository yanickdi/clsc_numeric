import sys
from math import sqrt
import scipy.optimize


MODEL_1, MODEL_2 = 1, 2
_CASE_ONE, _CASE_TWO = 1, 2
DECIMALS_ALLOW_NN = 15

def is_prof_pos(prof):
    """ checks whether a given profit is positive - it allows also a -.1*10^15 as positive! """
    return round(prof, 15) >= 0
    
class ModelTwoNumericalSolver:
    """
        This class offers methods to solve Model 2 (With Online Store of the Manufacturer) numerically
    """
    def __init__(self):
        pass
        
    
    def optimize(self, const_args):
        """
        This is the core method of this class. It will return all six 
        decision variables to maximize the retailer's profit (with respect to the profit maximization
        condition of the retailer) and the cond 0 <= qr <= (tau/roh)*qn
        
        Args:
            const_args (dict): A dictionary having set all parameters, e.g. 
        
        Returns:
            dict: A dictionary of decision_vars {pn, pr, wn, roh, qn, qr} or None if the solution is not possible
                  The dictionary also contains a key named `_dbg` where some calculation info is stored to test this method
        """
        raise NotImplementedError()
    
    def calc_profits(self, const_args, dec_vars):
        """
            Returns the numeric result of the profit of the manufacturer and the retailer (a tuple containing first manufacturer, second retailer)
            having set all decision variables
            
            This method checks whether `dec_vars` is not None. If its None - It will return a tuple of (None, None)
        """
        if dec_vars == None:
            return (None, None)
        tau, a, s, cr, cn, delta = const_args['tau'], const_args['a'], const_args['s'], const_args['cr'], const_args['cn'], const_args['delta']
        wn, pn, roh, qn, qr, pr = dec_vars['wn'], dec_vars['pn'], dec_vars['roh'], dec_vars['qn'], dec_vars['qr'], dec_vars['pr']
        manu_profit = qn * (wn * (1- tau/roh) - cn) + qr*(pr-cr) + ((tau/roh)*qn-qr)*s
        retailer_profit = None
        return manu_profit, retailer_profit
        
    def _optimize_case_one_a(self, const_args):
        """ helper function that solves the case roh >= 1 and qr = 0 """
        tau, a, s, cr, cn, delta = const_args['tau'], const_args['a'], const_args['s'], const_args['cr'], const_args['cn'], const_args['delta']
        dec = {
            'wn' : (cn+1)/2 - ((2-delta)/2) * sqrt((a*tau) / (1-delta)),
            'pr'  : ( delta*(-2*sqrt(tau*a*(1-delta))+ delta*(sqrt(tau*a*(1-delta))+2*delta-5) - cn*delta+cn+3 ) ) / (2*(delta-2)*(delta-1))
        }
        dec['roh'] = self.__roh_case_one(dec['wn'], delta, dec['pr'], tau, a)
        dec['pn'] = self.__pn_case_one(dec['pr'], delta, dec['pr'])
        dec['qn'] = self.__qn_case_one(dec['wn'], delta, dec['pr'])
        dec['qr'] = self.__qr_case_one(dec['wn'], delta, dec['pr'])
        return dec
        
    def _optimize_case_one_b(self, const_args):
        """ helper function that solves the case roh >= 1 and 0 < qr < tau/roh*qn """
        tau, a, s, cr, cn, delta = const_args['tau'], const_args['a'], const_args['s'], const_args['cr'], const_args['cn'], const_args['delta']
        dec = {
            'wn' : (1+cn)/2 -  ((2-delta)/2) * sqrt((a*tau)/(1-delta)),
            'pr'  : (cr+delta+s)/2 -  (delta/2) * sqrt((a*tau)/(1-delta))
        }
        dec['roh'] = self.__roh_case_one(dec['wn'], delta, dec['pr'], tau, a)
        dec['pn'] = self.__pn_case_one(dec['pr'], delta, dec['pr'])
        dec['qn'] = self.__qn_case_one(dec['wn'], delta, dec['pr'])
        dec['qr'] = self.__qr_case_one(dec['wn'], delta, dec['pr'])
        return dec
        
    def __qn_case_one(self, wn, delta, pr):
        return 1 - (self.__pn_case_one(wn, delta, pr) - pr)/(1-delta)
        
    def __qr(self, pn, pr, delta):
        return (pn-pr)/(1-delta) - pr/delta
        
    def __qr_case_one(self, wn, delta, pr):
        return self.__qr(self.__pn_case_one(wn, delta, pr), pr, delta)
        
    def __roh_case_one(self, wn, delta, pr, tau, a):
        return .5 * (1-delta+pr-wn)*sqrt(tau/(a*(1-delta)))
        
    def __pn_case_one(self, wn, delta, pr):
        return .5 * (1+wn-delta+pr)
    
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
            
            This method checks whether `dec_vars` is not None. If its None - It will return a tuple of (None, None)
        """
        if dec_vars == None:
            return (None, None)
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
        
        Returns:
        A dictionary of decision_vars {wn, pn, roh and qn} or None if the solution is not possible
        """
        ## test two cases:
        #       case 1 - roh is >= 1
        #       case 2 - roh is == 1
        
        dec_vars_case_1 = self._optimize_case_one(const_args)
        prof_man_case_1, prof_ret_case_1 = self.calc_profits(const_args, dec_vars_case_1)
        dec_vars_case_2 = self._optimize_case_two(const_args)
        prof_man_case_2, prof_ret_case_2 = self.calc_profits(const_args, dec_vars_case_2)
        
        case = None
        
        if dec_vars_case_1['roh'] < 1:
            if is_prof_pos(prof_man_case_2) and is_prof_pos(prof_ret_case_2):
                case = _CASE_TWO
            else:
                # case one and two not possible
                return None
        else:
            # roh is greater than 1, we have to check both -
            # the manufacturer decides:
            if is_prof_pos(prof_man_case_1) and is_prof_pos(prof_ret_case_1) and is_prof_pos(prof_man_case_2) and is_prof_pos(prof_ret_case_2):
                # both possible
                case = _CASE_ONE if prof_man_case_1 >= prof_man_case_2 else _CASE_TWO
            else:
                if is_prof_pos(prof_man_case_1) and is_prof_pos(prof_ret_case_1):
                    case = _CASE_ONE
                if is_prof_pos(prof_man_case_2) and is_prof_pos(prof_ret_case_2):
                    case = _CASE_TWO

        if case == None:
            return None
        else:
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
    
def check_args(model, args):
    if model == MODEL_1:
        tau, a, s, cn = args['tau'], args['a'], args['s'], args['cn']
        assert 0 <= tau <= 1
        assert 0 <= 0.01
        assert s <= cn <= 1
        assert 0 <= s <= cn
        assert s <= cn
    else:
        #TODO check model 2
        pass
    
    
def build_args(model, tau=None, a=None, s=None, cr=None, cn=None, delta=None):
    if model == MODEL_1:
        args =  {
            'tau' : tau,
            'a'   : a,
            's'   : s,
            'cn'  : cn
        }
    elif model == MODEL_2:
        args =  {
            'tau' : tau,
            'a'   : a,
            's'   : s,
            'cr'  : cr,
            'cn'  : cn,
            'delta' : delta
        }
    else:
        raise RuntimeError(str(model) + ' not allowed.')
    check_args(model, args)
    return args
    
    
if __name__ == '__main__':
    print('You cannot call this file directly')
    solver = ModelOneNumericalSolver()
    pass