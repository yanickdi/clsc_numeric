import sys
from math import sqrt
import scipy.optimize


MODEL_1, MODEL_2 = 1, 2
_CASE_ONE, _CASE_TWO = 1, 2
_CASE_ONE_A, _CASE_ONE_B, _CASE_ONE_C, _CASE_TWO_A, _CASE_TWO_B, _CASE_TWO_C = '1a','1b','1c','2a','2b','2c'
DECIMALS_ALLOW_NN = 14

def is_prof_pos(prof):
    """ checks whether a given profit is positive - it allows also a -.1*10^14 as positive! """
    return round(prof, 15) >= 0
    
def is_almost_equal(one, two):
    """ compares two floats, allowing a deviation after 14 decimal places """
    return round(one, DECIMALS_ALLOW_NN) == round(two, DECIMALS_ALLOW_NN) \
        or (0 in (round(one, DECIMALS_ALLOW_NN), round(two, DECIMALS_ALLOW_NN)) and \
            abs(round(one, DECIMALS_ALLOW_NN)) == abs(round(two, DECIMALS_ALLOW_NN)))
    
class ModelTwoNumericalSolver:
    """
        This class offers methods to solve Model 2 (With Online Store of the Manufacturer) numerically
    """
    def __init__(self):
        pass
        
    
    def optimize(self, par):
        """
        This is the core method of this class. It will return a Solution object having stored
        the profit of the manufacturer, the retailer, the case that led to the solution and all decision vars
        
        Args:
            par (Parameter): A Parameter object of type MODEL_2
        
        Returns:
            Solution: An object of class Solution
                      or None if there is no solution possible
        """
        cases = (
            (_CASE_ONE_A, self._optimize_case_one_a),
            (_CASE_ONE_B, self._optimize_case_one_b),
            (_CASE_ONE_C, self._optimize_case_one_c),
            (_CASE_TWO_A, self._optimize_case_two_a),
            (_CASE_TWO_B, self._optimize_case_two_b),
            (_CASE_TWO_C, self._optimize_case_two_c))
        
        valid_solutions = []
        for case_id, case_fct in cases:
            try:
                dec = case_fct(par)
                profit_man, profit_ret = self.calc_profits(par, dec)
                sol = Solution(dec, profit_man, profit_ret, case_id)
            except ZeroDivisionError:
                # there occured a zero division, so this solution will be skipped
                sol = None
            # check valid
            if sol is not None and self._is_valid(par, sol):
                valid_solutions.append(sol)
        if len(valid_solutions) > 0:
            # take the best valid solution (manufacturer decides)
            return max(valid_solutions, key=lambda sol: (sol.profit_man, sol.dec.lambda1, sol.dec.lambda2))
        else:
            return None
    
    def _is_valid(self, par, sol):
        """ Tests whether a given solution is feasible regarding to all model subjects """
        #TODO: assert all decision vars are positive in case of valid solution
        # check all variables positive
        for var in (sol.dec.pn, sol.dec.pr, sol.dec.wn, sol.dec.qn, sol.dec.qr):
            if var < -10**-DECIMALS_ALLOW_NN:
                return False
        # check case constraints
        if not (sol.dec.roh >= 1):
            return False
        if sol.case in (_CASE_ONE_A, _CASE_TWO_A):
            if not is_almost_equal(sol.dec.qr, 0):
                raise Exception()
                return False
        elif sol.case in (_CASE_ONE_B, _CASE_TWO_B):
            if not (-10**-DECIMALS_ALLOW_NN <= sol.dec.qr <= ((par.tau/sol.dec.roh) * sol.dec.qn)+10**-DECIMALS_ALLOW_NN):
                return False
        elif sol.case in (_CASE_ONE_C, _CASE_TWO_C):
            if not (is_almost_equal(sol.dec.qr, (par.tau/sol.dec.roh) * sol.dec.qn)):
                return False
                
        # check profits
        if not (sol.profit_man >= -10**-DECIMALS_ALLOW_NN and sol.profit_ret >= -10**-DECIMALS_ALLOW_NN):
            return False
        return True
    
    def calc_profits(self, par, dec):
        """
            Returns the numeric result of the profit of the manufacturer and the retailer (a tuple containing first manufacturer, second retailer)
            having set all decision variables
            
            This method checks whether `dec` is not None. If its None - It will return a tuple of (None, None)
        """
        if dec == None:
            return (None, None)
        manu_profit = dec.qn * (dec.wn * (1- par.tau/dec.roh) - par.cn) + dec.qr*(dec.pr-par.cr) + ((par.tau/dec.roh)*dec.qn-dec.qr)*par.s
        retailer_profit = dec.qn * (dec.pn - dec.wn) * (1- par.tau/dec.roh) - par.a*dec.roh
        return manu_profit, retailer_profit
        
    def _optimize_case_one_a(self, par):
        """ helper function that solves the case roh >= 1 and qr = 0 """
        dec = DecisionVariables(MODEL_2,
            wn = (par.cn+1)/2 - ((2-par.delta)/2) * sqrt((par.a*par.tau) / (1-par.delta)),
            pr = ( par.delta*(-2*sqrt(par.tau*par.a*(1-par.delta))+ par.delta*(sqrt(par.tau*par.a*(1-par.delta))+2*par.delta-5) - par.cn*par.delta+par.cn+3 ) ) / (2*(par.delta-2)*(par.delta-1)),
            lambda1 = 1 + par.cn + par.cr + par.s + (2 * (-1 + par.cn))/(-2 + par.delta) - par.delta,
            lambda2 = 0
            )
        dec.roh = self.__roh_case_one(dec.wn, par.delta, dec.pr, par.tau, par.a)
        dec.pn = self.__pn_case_one(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_one(dec.wn, par.delta, dec.pr)
        #dec.qr = self.__qr_case_one(dec.wn, par.delta, dec.pr)
        dec.qr = 0
        return dec
        
    def _optimize_case_one_b(self, par):
        """ helper function that solves the case roh >= 1 and 0 <= qr <= (tau/roh)*qn """
        dec = DecisionVariables(MODEL_2,
            wn = (1+par.cn)/2 -  ((2-par.delta)/2) * sqrt((par.a*par.tau)/(1-par.delta)),
            pr = (par.cr+par.delta+par.s)/2 -  (par.delta/2) * sqrt((par.a*par.tau)/(1-par.delta)),
            lambda1=0, lambda2=0
        )
        dec.roh = self.__roh_case_one(dec.wn, par.delta, dec.pr, par.tau, par.a)
        dec.pn = self.__pn_case_one(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_one(dec.wn, par.delta, dec.pr)
        dec.qr = self.__qr_case_one(dec.wn, par.delta, dec.pr)
        return dec
        
    def _optimize_case_one_c(self, par):
        """ helper function that solves the case roh >= 1 and qr = (tau/roh)*qn """
        dec = DecisionVariables(MODEL_2,
            wn = (par.cn+1)/2 - ((2-par.delta)/2)*sqrt((par.a*par.tau)/(1-par.delta)),
            pr = (par.delta*(-6*sqrt(par.tau*par.a*(1-par.delta))+par.delta*(5*sqrt(par.tau*par.a*(1-par.delta))+2*par.delta-5)-par.cn*par.delta + par.cn + 3)) / ( 2*(par.delta-2)*(par.delta-1)),
            lambda1 = 0,
            lambda2 = (-par.cr * (-2 + par.delta) - par.s *(-2 + par.delta) + par.delta *(-par.cn - (-1 + par.delta)*(-1 + 4*par.a*sqrt(par.tau/(par.a - par.a *par.delta)))))/(-2 + par.delta)
        )
        dec.roh = self.__roh_case_one(dec.wn, par.delta, dec.pr, par.tau, par.a)
        dec.pn = self.__pn_case_one(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_one(dec.wn, par.delta, dec.pr)
        dec.qr = self.__qr_case_one(dec.wn, par.delta, dec.pr)
        return dec
        
    def _optimize_case_two_a(self, par):
        """ helper function that solves the case roh = 1 and qr = 0 """
        dec = DecisionVariables(MODEL_2,
            wn = (1/(1-par.tau))*((1+par.cn)/2 - (par.tau*(1+par.s))/2),
            pr = (par.delta * (par.cn+ 2*par.delta*(par.tau-1)-(par.s+3)*par.tau + 3) ) / (2*(par.delta-2)*(par.tau-1)),
            qr = 0, roh = 1,
            lambda1 = (2 * par.cr * (-2 + par.delta) * (-1 + par.tau) + par.delta * (-2 + 2 * par.delta + par.cn * (-2 + par.tau) + par.tau - 2 * par.delta * par.tau + par.tau**2) - par.s * (4 * (-1 + par.tau) + par.delta * (2 + (-4 + par.tau) * par.tau)))/(2 * (-2 + par.delta) * (-1 + par.tau)),
            lambda2 = 0
        )
        dec.pn = self.__pn_case_two(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_two(dec.wn, par.delta, dec.pr)
        return dec
        
    def _optimize_case_two_b(self, par):
        """ helper function that solves the case roh = 1 and 0 <= qr <= (tau/roh)*qn """
        dec = DecisionVariables(MODEL_2,
            wn =  (par.tau*(-par.delta*(par.cn+5*par.s+5)-par.cr*(par.delta-2)+par.delta**2+6*par.s+4)+4*(par.cn+1)*(par.delta-1)+par.delta*par.s*par.tau**2) / (par.delta*((par.tau-8)*par.tau+8)+8*(par.tau-1)),
            pr = (-par.tau*(par.delta*(par.cn+5*par.delta+3*par.s-5)+par.cr*(3*par.delta-4)-4*par.s)+4*(par.delta-1)*(par.cr+par.delta+par.s)+par.delta*par.tau**2*(par.delta+par.s-1)) / (par.delta*((par.tau-8)*par.tau+8)+8*(par.tau-1)),
            lambda1 = 0, lambda2 = 0
            )
        dec.roh = self.__roh_case_two()
        dec.pn = self.__pn_case_two(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_two(dec.wn, par.delta, dec.pr)
        dec.qr = self.__qr_case_two(dec.wn, par.delta, dec.pr)
        return dec
    
    def _optimize_case_two_c(self, par):
        """ helper function that solves the case roh = 1 and qr = (tau/roh) * qn """
        dec = DecisionVariables(MODEL_2,
            wn =  (par.cn*(par.delta*(par.tau-1)+2)+par.delta*(par.tau*(par.cr*(par.tau-1)+par.delta*(-par.tau)+par.delta+par.tau+2)-1)+2*(par.cr-1)*par.tau+2) / (par.delta*(6*par.tau-2)-4*par.tau+4),
            pr = (par.delta*(par.tau*(par.cn+par.cr+5*par.delta-4)+par.cn+par.tau**2*(par.cr-par.delta+1)-2*par.delta+3)) / (par.delta*(6*par.tau-2)-4*par.tau+4),
            lambda1 = 0,
            lambda2 = (-4 * (par.cr + par.s) + 2 * (1 + par.cn + par.cr + par.s) * par.delta - 2 * par.delta**2 + (4 * (par.cr + par.s) + (-3 + par.cn - 4 * par.cr - 6 * par.s) * par.delta + 4 * par.delta**2) * par.tau + (1 + par.cr - par.delta) * par.delta * par.tau**2)/(4 - 4 * par.tau + par.delta * (-2 + 6 * par.tau))
            )
        dec.roh = self.__roh_case_two()
        dec.pn = self.__pn_case_two(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_two(dec.wn, par.delta, dec.pr)
        dec.qr = self.__qr_case_two(dec.wn, par.delta, dec.pr)
        return dec
        
    def __qr(self, pn, pr, delta):
        return (pn-pr)/(1-delta) - pr/delta
        
    def __qn(self, pn, pr, delta):
        return 1 - (pn-pr)/(1-delta)
        
    def __qn_case_one(self, wn, delta, pr):
        return self.__qn(self.__pn_case_one(wn, delta, pr), pr, delta)
        
    def __qr_case_one(self, wn, delta, pr):
        return self.__qr(self.__pn_case_one(wn, delta, pr), pr, delta)
        
    def __roh_case_one(self, wn, delta, pr, tau, a):
        return .5 * (1-delta+pr-wn)*sqrt(tau/(a*(1-delta)))
        
    def __pn_case_one(self, wn, delta, pr):
        return .5 * (1+wn-delta+pr)
        
    def __roh_case_two(self):
        return 1
    
    def __pn_case_two(self, wn, delta, pr):
        return .5 * (1+wn-delta+pr)
        
    def __qn_case_two(self, wn, delta, pr):
        return self.__qn(self.__pn_case_two(wn, delta, pr), pr, delta)
        
    def __qr_case_two(self, wn, delta, pr):
        return self.__qr(self.__pn_case_two(wn, delta, pr), pr, delta)
    
class ModelOneNumericalSolver:
    """
        This class offers methods to solve Model 1 (Without Online Store of the Manufacturer) numerically
    """
    def __init__(self):
        pass
    
    def calc_profits(self, par, dec_vars):
        """
            Returns the numeric result of the profit of the manufacturer and the retailer (a tuple containing first manufacturer, second retailer)
            having set all decision variables
            
            This method checks whether `dec_vars` is not None. If its None - It will return a tuple of (None, None)
        """
        if dec_vars == None:
            return (None, None)
        wn, pn, roh, qn = dec_vars.wn, dec_vars.pn, dec_vars.roh, dec_vars.qn
        manu_profit = qn * (wn * (1- par.tau/roh) - par.cn + (par.tau/roh) * par.s)
        retailer_profit = qn * (pn - wn) * (1 - par.tau/roh) - par.a * roh
        return manu_profit, retailer_profit
    
    def optimize(self, par):
        """
        This is the core method of this class. It will return all four
        decision variables to maximize the retailer's profit (with
        respect to the profit maximization condition of the retailer)
        
        Returns:
        A Solution object or None if the solution is not possible
        """
        ## test two cases:
        #       case 1 - roh is >= 1
        #       case 2 - roh is == 1
        
        dec_vars_case_1 = self._optimize_case_one(par)
        prof_man_case_1, prof_ret_case_1 = self.calc_profits(par, dec_vars_case_1)
        dec_vars_case_2 = self._optimize_case_two(par)
        prof_man_case_2, prof_ret_case_2 = self.calc_profits(par, dec_vars_case_2)
        
        case = None
        
        if dec_vars_case_1.roh < 1:
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
                    
        if case == _CASE_ONE:
            return Solution(dec_vars_case_1, prof_man_case_1, prof_ret_case_1, case)
        elif case == _CASE_TWO:
            return Solution(dec_vars_case_2, prof_man_case_2, prof_ret_case_2, case)
        else:
            return None
    
    def _optimize_case_one(self, par):
        """
        Returns a dec_vars dict having wn,pn,roh and qn stored
        """
        
        # defining our helper function:
        def __manufacturer_derivation_case_1(wn):
            return (-1/2) - (par.cn/2) + wn + par.a * (par.tau/par.a)**(1/2)
        
        # let scipy do the job:
        opt = scipy.optimize.fsolve(__manufacturer_derivation_case_1, x0=0.5, full_output=True)
        wn = opt[0][0]
        pn = (1 + wn) / 2
        roh = (1-wn) * (1/2) * sqrt(par.tau/par.a)
        return DecisionVariables(MODEL_1, wn=wn, pn=pn, roh=roh, qn=1 - pn)
        
    def _optimize_case_two(self, par):
        """
        Returns a dec_vars dict having wn,pn,roh and qn stored
        """
        #helper function
        def __manufacturer_derivation_case_2(wn):
            return (1/2) * (-1 -par.cn -2*wn * (-1 + par.tau) + par.tau + par.s*par.tau)
        
        # hello scipy:
        opt = scipy.optimize.fsolve(__manufacturer_derivation_case_2, x0=0.5, full_output=True)
        wn = opt[0][0]
        pn = (1 + wn) / 2
        roh = 1
        return DecisionVariables(MODEL_1, wn=wn, pn=pn, roh=roh, qn=1 - pn)

class Solution:
    def __init__(self, dec, profit_man, profit_ret, case):
        self.dec, self.profit_man, self.profit_ret, self.case = dec, profit_man, profit_ret, case
        
    
class Parameter:
    """
        An object of this class is a struct like wrapper for all Model Input Parameter (constants)
    """
    
    def __init__(self, model, tau=None, a=None, s=None, cr=None, cn=None, delta=None):
        if model == MODEL_1:
            self.tau, self.a, self.s, self.cn = tau, a, s, cn
        elif model == MODEL_2:
            self.tau, self.a, self.s, self.cr, self.cn, self.delta = tau, a, s, cr, cn, delta
        else:
            raise RuntimeError(str(model) + ' not allowed.')
        self.model = model
            
    def is_valid(self):
        #TODO
        raise NotImplementedError()
        
    def __str__(self):
        if self.model == MODEL_1:
            return 'tau={:.2f}, a={:.2f}, s={:.2f}, cn={:.2f}'.format(self.tau, self.a, self.s, self.cn)
        elif self.model == MODEL_2:
            return 'tau={:.2f}, a={:.2f}, s={:.2f}, cr={:.2f}, cn={:.2f}, delta={:.2f}'.format(
                self.tau, self.a, self.s, self.cr, self.cn, self.delta)
        else:
            return '?'
            
    def __repr__(self):
        return self.__str__()
            
class DecisionVariables:
    """
        An object of this class is a struct like wrapper for all Model Decision variables
    """
    
    def __init__(self, model, pn=None, pr=None, wn=None, roh=None, qn=None, qr=None, lambda1=None, lambda2=None):
        if model == MODEL_1:
            self.wn, self.pn, self.roh, self.qn = wn, pn, roh, qn
        elif model == MODEL_2:
            self.pn, self.pr, self.wn, self.roh, self.qn, self.qr = pn, pr, wn, roh, qn, qr
            self.lambda1, self.lambda2 = lambda1, lambda2
        else:
            raise RuntimeError(str(model) + ' not allowed.')
        self.model = model
        
    def __str__(self):
        if self.model == MODEL_1:
            return 'wn={:.5f}, pn={:.5f}, roh={:.5f}, qn={:.5f}'.format(self.wn, self.pn, self.roh, self.qn)
        else:
            return 'pn={:.5f}, pr={:.5f}, wn={:.5f}, roh={:.5f}, qn={:.5f}, qr={:.5f}'.format(self.pn, self.pr, self.wn, self.roh, self.qn, self.qr)
    
    def __repr__(self):
        return self.__str__()
    
if __name__ == '__main__':
    sys.exit('You cannot call this file directly')