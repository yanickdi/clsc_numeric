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
        
    
    def optimize(self, par):
        """
        This is the core method of this class. It will return all six 
        decision variables to maximize the retailer's profit (with respect to the profit maximization
        condition of the retailer) and the cond 0 <= qr <= (tau/roh)*qn
        
        Args:
            par (Parameter): A Parameter object of type MODEL_2
        
        Returns:
            DecisionVariables: An object of type MODEL_2 that stores all decision_vars (pn, pr, wn, roh, qn, qr) or None if the solution is not possible
        """
        raise NotImplementedError()
    
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
            pr = ( par.delta*(-2*sqrt(par.tau*par.a*(1-par.delta))+ par.delta*(sqrt(par.tau*par.a*(1-par.delta))+2*par.delta-5) - par.cn*par.delta+par.cn+3 ) ) / (2*(par.delta-2)*(par.delta-1))
            )
        dec.roh = self.__roh_case_one(dec.wn, par.delta, dec.pr, par.tau, par.a)
        dec.pn = self.__pn_case_one(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_one(dec.wn, par.delta, dec.pr)
        dec.qr = self.__qr_case_one(dec.wn, par.delta, dec.pr)
        return dec
        
    def _optimize_case_one_b(self, par):
        """ helper function that solves the case roh >= 1 and 0 < qr < tau/roh*qn """
        dec = DecisionVariables(MODEL_2,
            wn = (1+par.cn)/2 -  ((2-par.delta)/2) * sqrt((par.a*par.tau)/(1-par.delta)),
            pr = (par.cr+par.delta+par.s)/2 -  (par.delta/2) * sqrt((par.a*par.tau)/(1-par.delta))
        )
        dec.roh = self.__roh_case_one(dec.wn, par.delta, dec.pr, par.tau, par.a)
        dec.pn = self.__pn_case_one(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_one(dec.wn, par.delta, dec.pr)
        dec.qr = self.__qr_case_one(dec.wn, par.delta, dec.pr)
        return dec
        
    def _optimize_case_one_c(self, par):
        """ helper function that solves the case roh >= 1 and qr = tau/roh*qn """
        dec = DecisionVariables(MODEL_2,
            wn = (par.cn+1)/2 - ((2-par.delta)/2)*sqrt((par.a*par.tau)/(1-par.delta)),
            pr = (par.delta*(-6*sqrt(par.tau*par.a*(1-par.delta))+par.delta*(5*sqrt(par.tau*par.a*(1-par.delta))+2*par.delta-5)-par.cn*par.delta + par.cn + 3)) / ( 2*(par.delta-2)*(par.delta-1))
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
            pr = (par.delta * (par.cn+ 2*par.delta*(par.tau-1)-(par.s+3)*par.tau + 3) ) / (2*(par.delta-2)*(par.tau-1))
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
        A DecisionVariables object of type MODEL_1 or None if the solution is not possible
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

        if case == None:
            return None
        else:
            return dec_vars_case_1 if case == _CASE_ONE else dec_vars_case_2
    
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
        else:
            return 'Parameter obj of type model 2'
            
class DecisionVariables:
    """
        An object of this class is a struct like wrapper for all Model Decision variables
    """
    
    def __init__(self, model, pn=None, pr=None, wn=None, roh=None, qn=None, qr=None):
        if model == MODEL_1:
            self.wn, self.pn, self.roh, self.qn = wn, pn, roh, qn
        elif model == MODEL_2:
            self.pn, self.pr, self.wn, self.roh, self.qn, self.qr = pn, pr, wn, roh, qn, qr
        else:
            raise RuntimeError(str(model) + ' not allowed.')
        self.model = model
        
    def __str__(self):
        if self.model == MODEL_1:
            return 'wn={:.5f}, pn={:.5f}, roh={:.5f}, qn={:.5f}'.format(self.wn, self.pn, self.roh, self.qn)
        else:
            return 'pn={:.5f}, pr={:.5f}, wn={:.5f}, roh={:.5f}, qn={:.5f}, qr={:.5f}'.format(self.pn, self.pr, self.wn, self.roh, self.qn, self.qr)
    
if __name__ == '__main__':
    print('You cannot call this file directly')
    solver = ModelOneNumericalSolver()
    pass