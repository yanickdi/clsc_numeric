import sys
from math import sqrt
import sqlite3

import numpy as np
from scipy.optimize import fsolve, root, minimize

import utils


MODEL_1, MODEL_2, MODEL_2_QUAD, MODEL_1_QUAD, MODEL_NB = 1, 2, 3, 4, 5
_CASE_ONE, _CASE_TWO = '1', '2'
_UNDEFINED_CASE = '0'
_CASE_ONE_A, _CASE_ONE_B, _CASE_ONE_C, _CASE_TWO_A, _CASE_TWO_B, _CASE_TWO_C = '1a','1b','1c','2a','2b','2c'
ALL_CASES = [_CASE_ONE_A, _CASE_ONE_B, _CASE_ONE_C, _CASE_TWO_A, _CASE_TWO_B, _CASE_TWO_C]
DECIMALS_ALLOW_NN = 10

def is_prof_pos(prof):
    """ checks whether a given profit is positive - it allows also a -.1*10^14 as positive! """
    return round(prof, 15) >= 0
    
def is_almost_equal(one, two):
    """ compares two floats, allowing a deviation after 14 decimal places """
    return round(one, DECIMALS_ALLOW_NN) == round(two, DECIMALS_ALLOW_NN) \
        or (0 in (round(one, DECIMALS_ALLOW_NN), round(two, DECIMALS_ALLOW_NN)) and \
            abs(round(one, DECIMALS_ALLOW_NN)) == abs(round(two, DECIMALS_ALLOW_NN)))
            
class ModelTwoFracQuad:
    """
        This class offers methods to solve the Model 2 where rho is quadratic in the retailer function
    """
    
    @staticmethod
    def create_dec_vars(wn, pr, case, par):
        try:
            if case in (_CASE_ONE_A, _CASE_ONE_B, _CASE_ONE_C):  #case 1: rho >= 1
                pn = ModelTwoFracQuad._pn_case_one(wn, pr, par)
                rho = ModelTwoFracQuad._rho_case_one(wn, pr, par)
            else: #case 2: rho = 1
                pn = ModelTwoFracQuad._pn_case_two(wn, pr, par)
                rho = 1
            qn = ModelTwoFracQuad._qn(pn, pr, par)
            qr = ModelTwoFracQuad._qr(pn, pr, par)
        except ZeroDivisionError as e:
            raise e
        dec = DecisionVariables(MODEL_2, pn=pn, pr=pr, wn=wn, rho=rho, qn=qn, qr=qr)
        return dec
    
    @staticmethod
    def calc_profits(par, dec):
        if dec == None:
            return (None, None)
        
        manu_profit = dec.qn * (dec.wn * (1- par.tau/dec.rho - par.cn)) + dec.qr * (dec.pr - par.cr) + ((par.tau/dec.rho)*dec.qn - dec.qr)*par.s
        retailer_profit = dec.qn * (dec.pn - dec.wn) * (1 - par.tau/dec.rho) - (par.a*dec.rho - 1)**2
        return manu_profit, retailer_profit
        
    @staticmethod
    def _qn(pn, pr, par):
        return 1 - (pn - pr) / (1 - par.delta)
        
    def _qr(pn, pr, par):
        return (pn - pr) / (1 - par.delta) - pr /par.delta
        
    @staticmethod
    def _pn_case_one(wn, pr, par):
        return .5 * (1-par.delta+pr+wn)
    
    @staticmethod
    def _rho_case_one(wn, pr, par):
        return (par.tau**(1/3) * (-1+par.delta-pr+wn)**(2/3))/(2 * (par.a *(1-par.delta))**(1/3))
        
    @staticmethod
    def _pn_case_two(wn, pr, par):
        return .5 * (1-par.delta+pr+wn)
        
    @staticmethod
    def is_valid(par, sol):
        """ Tests whether a given solution is feasible regarding to all model subjects """
        # check all variables positive
        for var in (sol.dec.pn, sol.dec.pr, sol.dec.wn, sol.dec.qn, sol.dec.qr):
            if var < -10**-DECIMALS_ALLOW_NN:
                return False
                
        # check profits
        if not (sol.profit_man >= -10**-DECIMALS_ALLOW_NN and sol.profit_ret >= -10**-DECIMALS_ALLOW_NN):
            return False
            
        # check rho
        if not (sol.dec.rho >= 1):
            return False
        # check qr
        if not (-10**-DECIMALS_ALLOW_NN <= sol.dec.qr <= ((par.tau/sol.dec.rho) * sol.dec.qn)+10**-DECIMALS_ALLOW_NN):
            return False
            
        return True
        
    @staticmethod
    def is_valid_and_case_exact(par, sol):
        """ Tests whether a given solution is feasible regarding to all model subjects and right case """
        if not is_valid(par, sol):
            return False
            
        # check case constraints
        if not (sol.dec.rho >= 1):
            return False
        if sol.case in (_CASE_ONE_A, _CASE_TWO_A):
            if not is_almost_equal(sol.dec.qr, 0):
                raise Exception()
                return False
        elif sol.case in (_CASE_ONE_B, _CASE_TWO_B):
            if not (-10**-DECIMALS_ALLOW_NN <= sol.dec.qr <= ((par.tau/sol.dec.rho) * sol.dec.qn)+10**-DECIMALS_ALLOW_NN):
                return False
        elif sol.case in (_CASE_ONE_C, _CASE_TWO_C):
            if not (is_almost_equal(sol.dec.qr, (par.tau/sol.dec.rho) * sol.dec.qn)):
                return False
        return True
    
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
        if par.cn < par.cr: return None
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
            #valid_solutions.sort(key=lambda sol: )
            # take the best valid solution (manufacturer decides)
            return max(valid_solutions, key=lambda sol: (sol.profit_man, sol.profit_ret))
        else:
            return None
    
    def _is_valid(self, par, sol):
        """ Tests whether a given solution is feasible regarding to all model subjects """
        # check all variables positive
        for var in (sol.dec.pn, sol.dec.pr, sol.dec.wn, sol.dec.qn, sol.dec.qr):
            if var < -10**-DECIMALS_ALLOW_NN:
                return False
        # check case constraints
        if not (sol.dec.rho >= 1):
            return False
        if sol.case in (_CASE_ONE_A, _CASE_TWO_A):
            if not is_almost_equal(sol.dec.qr, 0):
                raise Exception()
                return False
        elif sol.case in (_CASE_ONE_B, _CASE_TWO_B):
            if not (-10**-DECIMALS_ALLOW_NN <= sol.dec.qr <= ((par.tau/sol.dec.rho) * sol.dec.qn)+10**-DECIMALS_ALLOW_NN):
                return False
        elif sol.case in (_CASE_ONE_C, _CASE_TWO_C):
            if not (is_almost_equal(sol.dec.qr, ((par.tau/sol.dec.rho) * sol.dec.qn))):
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
        manu_profit = dec.qn * (dec.wn * (1- par.tau/dec.rho) - par.cn) + dec.qr*(dec.pr-par.cr) + ((par.tau/dec.rho)*dec.qn-dec.qr)*par.s
        retailer_profit = dec.qn * (dec.pn - dec.wn) * (1- par.tau/dec.rho) - par.a*(dec.rho-1)
        return manu_profit, retailer_profit
        
    def _optimize_case_one_a(self, par):
        """ helper function that solves the case rho >= 1 and qr = 0 """
        dec = DecisionVariables(MODEL_2,
            wn = (par.cn+1)/2 - ((2-par.delta)/2) * sqrt((par.a*par.tau) / (1-par.delta)),
            pr = ( par.delta*(-2*sqrt(par.tau*par.a*(1-par.delta))+ par.delta*(sqrt(par.tau*par.a*(1-par.delta))+2*par.delta-5) - par.cn*par.delta+par.cn+3 ) ) / (2*(par.delta-2)*(par.delta-1)),
            lambda1 = 1 + par.cn + par.cr + par.s + (2 * (-1 + par.cn))/(-2 + par.delta) - par.delta,
            lambda2 = 0
            )
        dec.rho = self.__rho_case_one(dec.wn, par.delta, dec.pr, par.tau, par.a)
        dec.pn = self.__pn_case_one(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_one(dec.wn, par.delta, dec.pr)
        #dec.qr = self.__qr_case_one(dec.wn, par.delta, dec.pr)
        dec.qr = 0
        return dec
        
    def _optimize_case_one_b(self, par):
        """ helper function that solves the case rho >= 1 and 0 <= qr <= (tau/rho)*qn """
        dec = DecisionVariables(MODEL_2,
            wn = (1+par.cn)/2 -  ((2-par.delta)/2) * sqrt((par.a*par.tau)/(1-par.delta)),
            pr = (par.cr+par.delta+par.s)/2 -  (par.delta/2) * sqrt((par.a*par.tau)/(1-par.delta)),
            lambda1=0, lambda2=0
        )
        dec.rho = self.__rho_case_one(dec.wn, par.delta, dec.pr, par.tau, par.a)
        dec.pn = self.__pn_case_one(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_one(dec.wn, par.delta, dec.pr)
        dec.qr = self.__qr_case_one(dec.wn, par.delta, dec.pr)
        return dec
        
    def _optimize_case_one_c(self, par):
        """ helper function that solves the case rho >= 1 and qr = (tau/rho)*qn """
        dec = DecisionVariables(MODEL_2,
            wn = (par.cn+1)/2 - ((2-par.delta)/2)*sqrt((par.a*par.tau)/(1-par.delta)),
            pr = (par.delta*(-6*sqrt(par.tau*par.a*(1-par.delta))+par.delta*(5*sqrt(par.tau*par.a*(1-par.delta))+2*par.delta-5)-par.cn*par.delta + par.cn + 3)) / ( 2*(par.delta-2)*(par.delta-1)),
            lambda1 = 0,
            lambda2 = (-par.cr * (-2 + par.delta) - par.s *(-2 + par.delta) + par.delta *(-par.cn - (-1 + par.delta)*(-1 + 4*par.a*sqrt(par.tau/(par.a - par.a *par.delta)))))/(-2 + par.delta)
        )
        dec.rho = self.__rho_case_one(dec.wn, par.delta, dec.pr, par.tau, par.a)
        dec.pn = self.__pn_case_one(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_one(dec.wn, par.delta, dec.pr)
        dec.qr = self.__qr_case_one(dec.wn, par.delta, dec.pr)
        return dec
        
    def _optimize_case_two_a(self, par):
        """ helper function that solves the case rho = 1 and qr = 0 """
        dec = DecisionVariables(MODEL_2,
            wn = (1/(1-par.tau))*((1+par.cn)/2 - (par.tau*(1+par.s))/2),
            pr = (par.delta * (par.cn+ 2*par.delta*(par.tau-1)-(par.s+3)*par.tau + 3) ) / (2*(par.delta-2)*(par.tau-1)),
            qr = 0, rho = 1,
            lambda1 = (2 * par.cr * (-2 + par.delta) * (-1 + par.tau) + par.delta * (-2 + 2 * par.delta + par.cn * (-2 + par.tau) + par.tau - 2 * par.delta * par.tau + par.tau**2) - par.s * (4 * (-1 + par.tau) + par.delta * (2 + (-4 + par.tau) * par.tau)))/(2 * (-2 + par.delta) * (-1 + par.tau)),
            lambda2 = 0
        )
        dec.pn = self.__pn_case_two(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_two(dec.wn, par.delta, dec.pr)
        return dec
        
    def _optimize_case_two_b(self, par):
        """ helper function that solves the case rho = 1 and 0 <= qr <= (tau/rho)*qn """
        dec = DecisionVariables(MODEL_2,
            wn =  (par.tau*(-par.delta*(par.cn+5*par.s+5)-par.cr*(par.delta-2)+par.delta**2+6*par.s+4)+4*(par.cn+1)*(par.delta-1)+par.delta*par.s*par.tau**2) / (par.delta*((par.tau-8)*par.tau+8)+8*(par.tau-1)),
            pr = (-par.tau*(par.delta*(par.cn+5*par.delta+3*par.s-5)+par.cr*(3*par.delta-4)-4*par.s)+4*(par.delta-1)*(par.cr+par.delta+par.s)+par.delta*par.tau**2*(par.delta+par.s-1)) / (par.delta*((par.tau-8)*par.tau+8)+8*(par.tau-1)),
            lambda1 = 0, lambda2 = 0
            )
        dec.rho = self.__rho_case_two()
        dec.pn = self.__pn_case_two(dec.wn, par.delta, dec.pr)
        dec.qn = self.__qn_case_two(dec.wn, par.delta, dec.pr)
        dec.qr = self.__qr_case_two(dec.wn, par.delta, dec.pr)
        return dec
    
    def _optimize_case_two_c(self, par):
        """ helper function that solves the case rho = 1 and qr = (tau/rho) * qn """
        dec = DecisionVariables(MODEL_2,
            wn =  (par.cn*(par.delta*(par.tau-1)+2)+par.delta*(par.tau*(par.cr*(par.tau-1)+par.delta*(-par.tau)+par.delta+par.tau+2)-1)+2*(par.cr-1)*par.tau+2) / (par.delta*(6*par.tau-2)-4*par.tau+4),
            pr = (par.delta*(par.tau*(par.cn+par.cr+5*par.delta-4)+par.cn+par.tau**2*(par.cr-par.delta+1)-2*par.delta+3)) / (par.delta*(6*par.tau-2)-4*par.tau+4),
            lambda1 = 0,
            lambda2 = (-4 * (par.cr + par.s) + 2 * (1 + par.cn + par.cr + par.s) * par.delta - 2 * par.delta**2 + (4 * (par.cr + par.s) + (-3 + par.cn - 4 * par.cr - 6 * par.s) * par.delta + 4 * par.delta**2) * par.tau + (1 + par.cr - par.delta) * par.delta * par.tau**2)/(4 - 4 * par.tau + par.delta * (-2 + 6 * par.tau))
            )
        dec.rho = self.__rho_case_two()
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
        
    def __rho_case_one(self, wn, delta, pr, tau, a):
        return .5 * (1-delta+pr-wn)*sqrt(tau/(a*(1-delta)))
        
    def __pn_case_one(self, wn, delta, pr):
        return .5 * (1+wn-delta+pr)
        
    def __rho_case_two(self):
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
        wn, pn, rho, qn = dec_vars.wn, dec_vars.pn, dec_vars.rho, dec_vars.qn
        manu_profit = qn * (wn * (1- par.tau/rho) - par.cn + (par.tau/rho) * par.s)
        retailer_profit = qn * (pn - wn) * (1 - par.tau/rho) - par.a * (rho-1)
        return manu_profit, retailer_profit
        
    def _is_valid(self, par, sol):
        """ Tests whether a given solution is feasible regarding to all model subjects """
        # check all variables positive
        for var in (sol.dec.pn, sol.dec.qn):
            if var < -10**-DECIMALS_ALLOW_NN:
                return False
        # check case constraints
        if not (sol.dec.rho >= 1):
            return False
                
        # check profits
        if not (sol.profit_man >= 10**-DECIMALS_ALLOW_NN and sol.profit_ret >= 10**-DECIMALS_ALLOW_NN):
            return False
        return True
    
    def optimize(self, par):
        """
        This is the core method of this class. It will return all four
        decision variables to maximize the retailer's profit (with
        respect to the profit maximization condition of the retailer)
        
        Returns:
        A Solution object or None if the solution is not possible
        """
        ## test two cases:
        #       case 1 - rho is >= 1
        #       case 2 - rho is == 1
        if par.a == 0: return None
        
        dec_vars_case_1 = self._optimize_case_one(par)
        prof_man_case_1, prof_ret_case_1 = self.calc_profits(par, dec_vars_case_1)
        dec_vars_case_2 = self._optimize_case_two(par)
        prof_man_case_2, prof_ret_case_2 = self.calc_profits(par, dec_vars_case_2)
        
        case = None
        sol1 = Solution(dec_vars_case_1, prof_man_case_1, prof_ret_case_1, _CASE_ONE)
        sol2 = Solution(dec_vars_case_2, prof_man_case_2, prof_ret_case_2, _CASE_TWO)
        valid_solutions = []
        for sol in (sol1, sol2):
            if self._is_valid(par, sol):
                valid_solutions.append(sol)
        if len(valid_solutions) > 0:
            # take the best valid solution (manufacturer decides)
            return max(valid_solutions, key=lambda sol: (sol.profit_man, sol.profit_ret))
        else:
            return None
        
        if dec_vars_case_1.rho < 1:
            if is_prof_pos(prof_man_case_2) and is_prof_pos(prof_ret_case_2):
                case = _CASE_TWO
            else:
                # case one and two not possible
                return None
        else:
            # rho is greater than 1, we have to check both -
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
        Returns a dec_vars dict having wn,pn,rho and qn stored
        """
        
        # defining our helper function:
        def __manufacturer_derivation_case_1(wn):
            return (-1/2) - (par.cn/2) + wn + par.a * (par.tau/par.a)**(1/2)
        
        # let scipy do the job:
        opt = fsolve(__manufacturer_derivation_case_1, x0=0.5, full_output=True)
        wn = opt[0][0]
        pn = (1 + wn) / 2
        rho = (1-wn) * (1/2) * sqrt(par.tau/par.a)
        return DecisionVariables(MODEL_1, wn=wn, pn=pn, rho=rho, qn=1 - pn)
        
    def _optimize_case_two(self, par):
        """
        Returns a dec_vars dict having wn,pn,rho and qn stored
        """
        #helper function
        def __manufacturer_derivation_case_2(wn):
            return (1/2) * (-1 -par.cn -2*wn * (-1 + par.tau) + par.tau + par.s*par.tau)
        
        # hello scipy:
        opt = fsolve(__manufacturer_derivation_case_2, x0=0.5, full_output=True)
        wn = opt[0][0]
        pn = (1 + wn) / 2
        rho = 1
        return DecisionVariables(MODEL_1, wn=wn, pn=pn, rho=rho, qn=1 - pn)
        
    @staticmethod
    def get_retailer_profit(par, wn, pn, rho):
        """ returns the the retailer's profit having set wn, pn and rho (i.e. all decision variables are set) """
        qn = 1 - pn
        return qn * (pn - wn) * (1-par.tau/rho) - par.a * (rho - 1)
        
    @staticmethod
    def get_manufacturer_profit(par, wn, case=_UNDEFINED_CASE):
        pn, rho = ModelOneNumericalSolver.get_retailer_decision(par, wn, case=case)
        if pn is None or rho is None: return None
        qn = 1 - pn
        return qn * (wn*(1-par.tau/rho) - par.cn + (par.tau/rho) * par.s)
        
    @staticmethod
    def get_retailer_decision(par, wn, case=_UNDEFINED_CASE):
        """ Returns a tuple of (pn, rho). Having wn set to a fixed value, retailer optimizes profits """
        pn_1 = (1/2) * (1+wn)
        rho_1 = (par.tau**(1/2) * (1-wn)) / (2*par.a**(1/2))
        pn_2 = (1/2) * (1+wn)
        rho_2 = 1
        all = []
        pn_rhos = ((pn_1, rho_1), (pn_2, rho_2))
        if case != _UNDEFINED_CASE:
            if case == _CASE_ONE:   pn_rhos = ((pn_1, rho_1),)
            elif case == _CASE_TWO: pn_rhos = ((pn_2, rho_2),)
        for pn, rho in pn_rhos:
            profit = ModelOneNumericalSolver.get_retailer_profit(par, wn, pn, rho)
            if 0 <= pn <= 1 and rho >= 1 and profit > 0:
                all.append([profit, (pn, rho)])
        if len(all) == 0: return None, None
        return max(all, key=lambda k: k[0])[1]    
    
    singleton = None
    @staticmethod
    def solve(par):
        if ModelOneNumericalSolver.singleton is None: ModelOneNumericalSolver.singleton = ModelOneNumericalSolver()
        return ModelOneNumericalSolver.singleton.optimize(par)
    

class Solution:
    def __init__(self, dec, profit_man, profit_ret, case):
        self.dec, self.profit_man, self.profit_ret, self.case = dec, profit_man, profit_ret, case
        
    def __str__(self):
        if self.dec.b is not None:
            return 'Solution object (wn={}, b={})'.format(self.dec.wn, self.dec.b)
        return 'Solution object (wn={}, pr={})'.format(self.dec.wn, self.dec.pr)
        

        
class Parameter:
    """
        An object of this class is a struct like wrapper for all Model Input Parameter (constants)
    """
    
    def __init__(self, model, tau=None, a=None, s=None, cr=None, cn=None, delta=None):
        #if model == MODEL_1 or model == MODEL_1_QUAD:
        #    self.tau, self.a, self.s, self.cn = tau, a, s, cn
        #    self.cr, self.delta = None, None
        #else:
        #    self.tau, self.a, self.s, self.cr, self.cn, self.delta = tau, a, s, cr, cn, delta
        self.tau, self.a, self.s, self.cr, self.cn, self.delta = tau, a, s, cr, cn, delta
        self.model = model
        
    def __str__(self):
        if self.model == MODEL_1 or self.model == MODEL_1_QUAD or self.model == MODEL_NB:
            return 'tau={:.2f}, a={:.4f}, s={:.2f}, cn={:.4f}'.format(self.tau, self.a, self.s, self.cn)
        else:
            return 'tau={:.4f}, a={:.5f}, s={:.4f}, cr={:.4f}, cn={:.4f}, delta={:.4f}'.format(
                self.tau, self.a, self.s, self.cr, self.cn, self.delta)
            
    def __repr__(self):
        return self.__str__()
            
class DecisionVariables:
    """
        An object of this class is a struct like wrapper for all Model Decision variables
    """
    
    def __init__(self, model, pn=None, pr=None, wn=None, rho=None, qn=None, qr=None, b=None, lambda1=None, lambda2=None):
        self.pn, self.pr, self.wn, self.rho, self.qn, self.qr, self.b = pn, pr, wn, rho, qn, qr, b
        self.lambda1, self.lambda2 = lambda1, lambda2
        #if model == MODEL_1 or model == MODEL_:
        #    self.wn, self.pn, self.rho, self.qn = wn, pn, rho, qn
        #    self.pr, self.qr = None, None
        #else:
        #    self.pn, self.pr, self.wn, self.rho, self.qn, self.qr = pn, pr, wn, rho, qn, qr
        #    self.lambda1, self.lambda2 = lambda1, lambda2
        self.model = model
        
    def __str__(self):
        if self.model == MODEL_1:
            return 'wn={:.5f}, pn={:.5f}, rho={:.5f}, qn={:.5f}'.format(self.wn, self.pn, self.rho, self.qn)
        elif self.model == MODEL_NB:
            return 'wn={:.5f}, b={:.3f}, pn={:.5f}, rho={:.5f}, qn={:.5f}'.format(self.wn, self.b, self.pn, self.rho, self.qn)
        else:
            return 'pn={:.5f}, pr={:.5f}, wn={:.5f}, rho={:.5f}, qn={:.5f}, qr={:.5f}'.format(self.pn, self.pr, self.wn, self.rho, self.qn, self.qr)
    
    def __repr__(self):
        return self.__str__()
        
def drange(start, end, step_size):
    """ A floating point range from [start, end] with step size step_size"""
    r = start
    while r <= end:
        yield r
        r += step_size
        
def almost_equal(first, sec, tol=0.001):
    return abs(first-sec) < tol
    
class ModelTwoGridTester:
    def __init__(self, par):
        self.par = par
        self.points = []
        
    def plot_profit_surface(self):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        sol = self.search()
        np_arr = np.zeros([len(self.points), 3])
        # create numpy arrays
        for i, (wn, pr, prof) in enumerate(self.points):
            prof_val = prof if prof is not None else -100
            np_arr[i, :] = [wn, pr, prof_val]
        best_points = np.zeros((len(self.caller_info), 2))
        for i, iter_inf in enumerate(self.caller_info):
            best_points[i:] = (iter_inf['best_x'], iter_inf['best_y'])
            
        from matplotlib import pyplot as plt
        valids = np_arr[np_arr[:, 2] > -100]
        non_valids =np_arr[np_arr[:, 2] <= -100]
        #plt.plot(valids[:, 0], valids[:, 1], linestyle='', marker='o', color='green')
        #ax.plot(valids[:, 0], valids[:, 1], valids[:, 2], linestyle='', marker='o', color='green')
        ax.plot_trisurf(valids[:, 0], valids[:, 1], valids[:, 2])
        #plt.plot(non_valids[:, 0], non_valids[:, 1], linestyle='', marker='x', color='red')
        #plt.plot(best_points[:, 0], best_points[:, 1], linestyle='', marker='*', color='blue')
        #ax.show3d()
        ax.plot([0.5385829601790], [0.4856783722067], [0.1607711822739], linestyle='', marker='o', color='blue')
        ax.plot([sol.dec.wn], [sol.dec.pr], [sol.profit_man], linestyle='', marker='o', color='red')
        #ax.set_xlim(0.5, 0.6)
        #ax.set_ylim(0.4, 0.6)
        #ax.set_zlim(0, 0.3)
        plt.show()
        
        
    def plot(self):
        sol = self.search()
        np_arr = np.zeros([len(self.points), 3])
        # create numpy arrays
        for i, (wn, pr, prof) in enumerate(self.points):
            prof_val = prof if prof is not None else -100
            np_arr[i, :] = [wn, pr, prof_val]
        best_points = np.zeros((len(self.caller_info), 2))
        for i, iter_inf in enumerate(self.caller_info):
            best_points[i:] = (iter_inf['best_x'], iter_inf['best_y'])
            
        from matplotlib import pyplot as plt
        valids = np_arr[np_arr[:, 2] > -100]
        non_valids =np_arr[np_arr[:, 2] <= -100]
        plt.plot(valids[:, 0], valids[:, 1], linestyle='', marker='o', color='green', label='looked at')
        plt.plot(non_valids[:, 0], non_valids[:, 1], linestyle='', marker='x', color='red', label='no sol')
        plt.plot(best_points[:, 0], best_points[:, 1], linestyle='', marker='*', color='blue', label='best fval in it')
        plt.plot(0.5385829601790, 0.4856783722067, label='global opt', marker='8', color='magenta', linestyle='')
        
        plt.legend()
        plt.show()
        
    def search(self):
        if self.par.a == 0: return None
        wn_range = [0, 1]
        pr_range = [0, 1]
        raster_size = 1
        raster_start_size = 201
        iter = 1
        self.caller_info = []
        sol_tuple = ImprovedGridSearch2D.maximize(self._grid_search_func,
            self.par, wn_range, pr_range, raster_size, iter, raster_start_size=raster_start_size, caller_info=self.caller_info)
        if sol_tuple is None: return None
        dec = DecisionVariables(MODEL_2_QUAD, pn=sol_tuple[2], pr=sol_tuple[1],
            wn=sol_tuple[0], rho=sol_tuple[3], qn=sol_tuple[4], qr=sol_tuple[5])
        # get case:
        qr_rel = dec.qr / (par.tau / dec.rho)*dec.qn
        if dec.rho > 1:
            case = _CASE_ONE
        else:
            case = _CASE_TWO
        sol = Solution(dec, sol_tuple[6], sol_tuple[7], case)
        return sol
        
    def _grid_search_func(self, par, wn, pr):
        ret_dec = ModelTwoGridSearch._retailer_decision(par, wn, pr)
        if ret_dec is not None:
            pn, rho, qn, qr, ret_prof = ret_dec
            man_profit = qn*(wn*(1-self.par.tau/rho)-par.cn) + qr*(pr-self.par.cr)+((self.par.tau/rho)*qn - qr)*self.par.s
            self.points.append([wn, pr, man_profit])
            return man_profit, (wn, pr, pn, rho, qn, qr, man_profit, ret_prof)
        else:
            self.points.append([wn, pr, None])
            return None, None
        
        
        
        
        
class ModelTwoQuadSolver:
    """ This class solves Model Two Quad by first using the GridSearch - then a local search algorithm (nelder-mead) """
    
    @staticmethod
    def solve(par, resolution='high'):
        if par.a == 0 or par.cn < par.cr: return None
        case_solutions = []
        if resolution == 'very-high':
            raster_size = 2000
            wn_range = [0.52, 0.55]
            pr_range = [0.48, 0.49]
            #search_cases = (_CASE_ONE, _CASE_TWO)
            # only case one interesting
            search_cases = (_CASE_ONE, )
            do_local_search = True
        elif resolution == 'low':
            raster_size = 51
            search_cases = (_UNDEFINED_CASE,) # all cases once
            do_local_search = False
            #wn_range = [0.52, 0.55]
            #pr_range = [0.48, 0.52]
            wn_range = [0, 1]
            pr_range = [0, 1]
        elif resolution == 'high' or (resolution == 'super high' and par.cn < 0.66):
            raster_size = 51
            search_cases = (_CASE_ONE, _CASE_TWO)  # all cases sepearte
            do_local_search = True
            wn_range, pr_range = [0, 1], [0, 1]
        elif resolution == 'super high' and par.cn >= 0.66:
            raster_size = 51
            search_cases = (_CASE_ONE, _CASE_TWO)  # all cases sepearte
            do_local_search = True
            wn_range, pr_range = [0.94, 1], [0, 0.3]
        # we have to check both cases:
        for case in search_cases:
            startP = ModelTwoQuadGridSearch.search(par, iter=1, raster_size=raster_size, case=case, wn_range=wn_range, pr_range=pr_range)
            if startP is None:
                #print('warning at point', par, ' in case', case, ',startP is None')
                continue
            start_vec = [startP.dec.wn, startP.dec.pr]
            # improve wn, pr
            result = minimize(ModelTwoQuadSolver._minimize_func, args={'par': par, 'case': case}, x0=start_vec,
                              method='Nelder-Mead', options={'xatol': 0.00000000000000001, 'fatol': 0.00000000000000001})
            # build a new solution using wn, pr
            sol = ModelTwoQuadSolver.getSolution(par, float(result.x[0]), float(result.x[1]))
            # assert that wn/pr combination lies really in case `case`
            if sol is not None:
                case_solutions.append(sol)

        if len(case_solutions) >= 1:
            return max(case_solutions, key=lambda k: k.profit_man)
        else: return None
        
    @staticmethod
    def _minimize_func(x, args):
        wn, pr = float(x[0]), float(x[1])
        par = args['par']
        case = args['case']
        man_profit, data = ModelTwoQuadSolver.profit(par, wn, pr, case)
        if man_profit is None:
            return sys.maxsize
        return -man_profit
        
    @staticmethod
    def getSolution(par, wn, pr, case=_UNDEFINED_CASE):
        man_profit, sol_tuple = ModelTwoQuadSolver.profit(par, wn, pr, case)
        if sol_tuple is None or man_profit <= 0 or sol_tuple[7] <= 0: return None
        dec = DecisionVariables(MODEL_2_QUAD, pn=sol_tuple[2], pr=sol_tuple[1],
            wn=sol_tuple[0], rho=sol_tuple[3], qn=sol_tuple[4], qr=sol_tuple[5])
        # get case:
        qr_rel = dec.qr / ((par.tau / dec.rho)*dec.qn)
        if dec.rho > 1:
            if   qr_rel < .01: case = _CASE_ONE_A
            elif qr_rel > .99: case = _CASE_ONE_C
            else: case = _CASE_ONE_B
        else:
            if   qr_rel < .01: case = _CASE_TWO_A
            elif qr_rel > .99: case = _CASE_TWO_C
            else: case = _CASE_TWO_B
        #case = _CASE_ONE if dec.rho > 1 else _CASE_TWO
        sol = Solution(dec, sol_tuple[6], sol_tuple[7], case)
        return sol
        
    @staticmethod
    def _retailer_decision(par, wn, pr, case=_UNDEFINED_CASE):
        assert type(wn) == float and type(pr) == float
        tau, a, s, cr, cn, delta = par.tau, par.a, par.s, par.cr, par.cn, par.delta
        rho_1 = -(((-1)**(2/3) * par.tau**(1/3) *(-1 + par.delta - pr + wn)**(2/3))/( 2 *(par.a *(-1 + par.delta))**(1/3)))
        pn_1  = .5 * (1+wn-delta+pr)
        rho_2 = 1
        pn_2 = .5 * (1+wn-delta+pr)
        valids = []

        if type(rho_1) == complex or type(rho_1) == np.complex128:
            if abs(rho_1.imag) < .000001 and rho_1.real > 1:
                rho_1 = rho_1.real
            else:
                rho_1 = -1
          
        for pn, rho in ((pn_1, rho_1), (pn_2, rho_2)):
            qn = 1 - (pn - pr)/(1-delta)
            qr = (pn-pr)/(1-delta) - pr/delta
            if rho == 0: continue
            ret_profit = qn*(pn-wn)*(1-par.tau/rho)-par.a*(rho**2-1)
            if (0 <= qr <= (par.tau/rho)*qn) and rho >= 1 and ret_profit >= 0:
                valids.append([pn, rho, qn, qr, ret_profit])
        if len(valids) > 0:
             best_valid = max(valids, key=lambda k: k[4])
             # only return solution if retailer decision is in case `case`
             best_rho = best_valid[1]
             if case == _CASE_ONE:
                return best_valid if best_rho > 1 else None
             if case == _CASE_TWO:
                return best_valid if best_rho == 1 else None
             if case == _UNDEFINED_CASE:
                return best_valid
        else:
            return None
        
    @staticmethod
    def profit(par, wn, pr, case=_UNDEFINED_CASE):
        ret_dec = ModelTwoQuadSolver._retailer_decision(par, wn, pr, case=case)
        if ret_dec is None: return None, None
        pn, rho, qn, qr, ret_prof = ret_dec
        man_profit = qn*(wn*(1-par.tau/rho)-par.cn) + qr*(pr-par.cr)+((par.tau/rho)*qn - qr)*par.s
        return man_profit, (wn, pr, pn, rho, qn, qr, man_profit, ret_prof)
        
class ModelTwoQuadGridSearch: 
    @staticmethod
    def search(par, iter=1, raster_size=51, case=_UNDEFINED_CASE, wn_range=[0, 1], pr_range=[0, 1]):
        if par.a == 0: return None
        raster_start_size = raster_size
        sol_tuple = ImprovedGridSearch2D.maximize(ModelTwoQuadGridSearch._grid_search_func,
            {'par': par, 'case': case}, wn_range, pr_range, raster_size, iter, raster_start_size=raster_start_size)
        if sol_tuple is None: return None
        dec = DecisionVariables(MODEL_2_QUAD, pn=sol_tuple[2], pr=sol_tuple[1],
            wn=sol_tuple[0], rho=sol_tuple[3], qn=sol_tuple[4], qr=sol_tuple[5])
        case = _CASE_ONE if sol_tuple[3] > 1 else _CASE_TWO
        sol = Solution(dec, sol_tuple[6], sol_tuple[7], case)
        return sol
    
    @staticmethod
    def _grid_search_func(options, wn, pr):
        par, case = options['par'], options['case']
        return ModelTwoQuadSolver.profit(par, wn, pr, case)        
        
        
        
class ModelOneQuadGridSearch:
    @staticmethod
    def _retailer_decision(par, wn):
        tau, a, s, cn = par.tau, par.a, par.s, par.cn
        #rho_1 = (par.tau**(1/3) * (-1 + wn)**(2/3))/(2 *par.a**(1/3))
        rho_1 = -(((-1)**(1/3) * par.tau**(1/3) * (-1 + wn)**(2/3))/(2 *par.a**(1/3)))
        #rho_1 = ((-1)**(2/3) * par.tau**(1/3) *(-1 + wn)**(2/3))/(2 * a**(1/3))
        pn_1  = (1+wn)/2
        rho_2 = 1
        pn_2 = (1+wn)/2
        valids = []
        for pn, rho in ((pn_1, rho_1), (pn_2, rho_2)):
            qn = 1 - pn
            if type(rho) == complex and round(rho.imag, 10) == 0: rho = rho.real
            if rho == 0: continue
            ret_profit = qn*(pn-wn)*(1-tau/rho)-a*(rho**2-1)
            if (rho >= 1 and ret_profit >= 0):
                valids.append([pn, rho, qn, ret_profit])
        if len(valids) > 0:
            return max(valids, key=lambda k: k[3])
        else:
            return None
        
    @staticmethod
    def _grid_search_func(par, wn, _):
        ret_dec = ModelOneQuadGridSearch._retailer_decision(par, wn)
        if ret_dec is None: return None, None
        pn, rho, qn, ret_prof = ret_dec
        man_profit = qn*(wn*(1-par.tau/rho)- par.cn + (par.tau/rho) * par.s)
        return man_profit, (wn, pn, rho, qn, man_profit, ret_prof)
        
    def search(self, par, resolution='high'):
        if par.a == 0: return None
        wn_range = [0, 1]
        y_range = [0, 0] #not there
        raster_size = 30
        iter = 20
        sol_tuple = GridSearch2D.maximize(ModelOneQuadGridSearch._grid_search_func,
            par, wn_range, y_range, raster_size, iter)
        if sol_tuple is None or sol_tuple[4] <= 0: return None
        dec = DecisionVariables(MODEL_1_QUAD, pn=sol_tuple[1], wn=sol_tuple[0],
                rho=sol_tuple[2], qn=sol_tuple[3])
        case = _CASE_ONE if sol_tuple[2] > 1 else _CASE_TWO
        sol = Solution(dec, sol_tuple[4], sol_tuple[5], case)
        if sol is not None and sol.profit_ret > 0.0001:
            return sol
        else: return None
    
    def search_old(self, par):
        """ Returns a solution object or None if no solution found """
        max_man_profit = -1
        #            0:wn  1:pn   2:rho  3:qn  4:ret_profit
        max_param = None
        for wn in np.arange(0, 1+.01, .01):
            wn = float(wn)
            ret_dec = self._retailer_decision(par, wn)
            if ret_dec is None: continue
            pn, rho, qn, ret_prof = ret_dec
            man_profit = qn*(wn*(1-par.tau/rho))- par.cn + (par.tau/rho) * par.s
            if man_profit > max_man_profit:
                max_man_profit = man_profit
                max_param = [wn, pn, rho, qn, ret_prof]
        
        # if we found something, build a solution object and return
        if max_man_profit >= 0:
            dec = DecisionVariables(MODEL_1_QUAD, pn=max_param[1], wn=max_param[0],
                rho=max_param[2], qn=max_param[3])
            case = _CASE_ONE if max_param[2] > 1 else _CASE_TWO
            sol = Solution(dec, max_man_profit, max_param[4], case)
            return sol
        else:
            return None

class ModelNBZeroSolver:
    @staticmethod
    def solve(par, resolution='high'):
        startP = ModelNBGridSearch.search(par, resolution, b_range=[0, 0])
        #return startP
        #try to improve
        if startP is None: return None
        start_vec = [startP.dec.wn]
        result = minimize(ModelNBSolver._minimize_func_but_set_b_to_zero,
            args={'par': par},
            x0=start_vec, method='Nelder-Mead')
        return ModelNBSolver.getSolution(par, float(result.x[0]), 0)

class ModelNBSolver:
    """ This class solves Model NB Quad by first using the GridSearch - then a local search algorithm (nelder-mead) """
    
    @staticmethod
    def solveOld(par, resolution='high'):
        raise RuntimeError()
        solPositiveB = ModelNBSolver.solvePositiveB(par, resolution)
        solZeroB = ModelNBZeroSolver.solve(par, resolution)
        
        if solPositiveB is not None and solPositiveB.profit_man > solZeroB.profit_man:
            diff = solPositiveB.profit_man - solZeroB.profit_man
            if diff < 0.000001:
                return solZeroB
        else:
            return solZeroB
            
    @staticmethod
    def solve(par, resolution='high'):
        if par.a == 0: return None
        if resolution == 'high':
            wn_lower_lim = par.cn
            wn_upper_lim = 1
            wn_count = 100
            b_count = 100
        elif resolution == 'case study':
            wn_lower_lim = 0.463
            wn_upper_lim = 0.555
            wn_count = 1000
            b_count = 3
        best_fval = -sys.maxsize
        best_fobj = None
        # case 1
        for wn in np.linspace(wn_lower_lim, wn_upper_lim, wn_count):
            for b in np.linspace(par.s, wn, b_count):
                wn, b = float(wn), float(b)
                # calc manufacturer profit
                man_prof_tuple = ModelNBSolver._manufacturer_profit(par, wn, b)
                if man_prof_tuple[0] is not None and man_prof_tuple[0] > best_fval:
                    # update
                    best_fval = man_prof_tuple[0]
                    best_fobj = man_prof_tuple[1]
        # case 2
        wn_2 = (-1 - par.cn + par.tau + par.s* par.tau)/(2 *(-1 + par.tau))
        b_2 = wn_2
        #b_2 = par.s
        #wn_2 = (1+par.cn)/2 - (par.tau * (1+par.s) )/2 + b_2 * par.tau
        pn_2 = (1 - par.tau - b_2 * par.tau + wn_2)/(2 * (1 - par.tau))
        qn_2 = 1 - pn_2
        profitr_case_2 = qn_2 * (pn_2*(1 - par.tau/1) - wn_2 + par.tau/ 1 * b_2)
        man_profit_case_2 = qn_2*(wn_2 - par.cn - (par.tau/1)*b_2 + (par.tau/1)*par.s)
        # use best solution (case one or two?)
        if man_profit_case_2 > best_fval and profitr_case_2 >= 0:
            best_fval = man_profit_case_2
            best_fobj = (wn_2, b_2, pn_2, 1, qn_2, man_profit_case_2, profitr_case_2)
        
        if best_fobj is None: return None
        wn, b, pn, rho, qn, man_profit, ret_prof = best_fobj
        dec = DecisionVariables(MODEL_NB, pn=pn, b=b, wn=wn, rho=rho, qn=qn)
        case = _CASE_ONE if rho > 1 else _CASE_TWO
        sol = Solution(dec, man_profit, ret_prof, case)
        
        #return sol
        # do not try localsearch
        #return ModelNBSolver.local_search(par, sol)
        
    @staticmethod
    def local_search(par, sol):
        start_vec = [sol.dec.wn, sol.dec.b]
        result = minimize(ModelNBSolver.__local_search_minimize_func,
            args={'par': par},
            x0=start_vec, method='Nelder-Mead')
        
        wn, b = result.x; wn, b = float(wn), float(b)
        fval, fobj = ModelNBSolver._manufacturer_profit(par, wn, b)
        wn, b, pn, rho, qn, man_profit, ret_prof = fobj
        dec = DecisionVariables(MODEL_NB, pn=pn, b=b, wn=wn, rho=rho, qn=qn)
        case = _CASE_ONE if rho > 1 else _CASE_TWO
        sol = Solution(dec, man_profit, ret_prof, case)
        return sol
            
    @staticmethod
    def __local_search_minimize_func(x, args):
        wn, b, par = float(x[0]), float(x[1]), args['par']
        fval, f_obj = ModelNBSolver._manufacturer_profit(par, wn, b)
        if fval is None:
            return sys.maxsize
            
        return -fval
                    
    def _manufacturer_profit(par, wn, b):
        if b < par.s: return None, None
        ret_prof_tuple = ModelNBSolver._retailer_profit(par, wn, b)
        if ret_prof_tuple is None: return None, None
        ret_prof, pn, rho = ret_prof_tuple
        qn = 1 - pn
        man_profit = qn*(wn - par.cn - (par.tau/rho)*b + (par.tau/rho)*par.s)
        return man_profit, (wn, b, pn, rho, qn, man_profit, ret_prof)
                
    @staticmethod
    def _retailer_profit(par, wn, b):
        valids = []
        for case in (1,):
            if case == 1:
                pn, rho = utils.model_nb_pn_rho_case_one(par, wn, b)
            qn = 1 - pn
            if not(0 <= qn <= 1 and rho >= 1):
                continue
            profitr = qn * (pn*(1 - par.tau/rho) - wn + par.tau/ rho * b) - par.a*(rho - 1)
            if profitr >= 0:
                valids.append( (profitr, pn, rho) )
        if len(valids) > 0:
            return max(valids, key=lambda k: k[0])
        return None
    
    @staticmethod
    def _pn_rho_case_one(par, wn, b):
        rho = (1 - b)* (1/((wn - b)**2 + 4*par.a/par.tau))**(1/2)
        pn = (rho - par.tau - b * par.tau + rho * wn)/(2 * (rho - par.tau))
        
        if wn == b:
            rho_model_n = (par.tau**(1/2) *(1 - wn)) /(2 * par.a**(1/2))
            pn_model_n = (1 + wn) / 2
            assert (rho_model_n - rho) < 0.00001 and (pn_model_n - pn) < 0.000001
        
        return pn, rho
        
    @staticmethod
    def _pn_rho_case_two(par, wn):
        pn = (1 - par.tau - wn * par.tau + wn)/(2 * (1 - par.tau))
        rho = 1
        pn_model_n = (1+wn)/2
        
        assert(pn - pn_model_n < 0.00001)
        return pn_model_n, rho
    
    @staticmethod
    def solvePositiveB(par, resolution='high'):
        startP = ModelNBGridSearch.search(par, resolution)
        if resolution == 'middle':
            # also try to improve:
            start_vec = [startP.dec.wn, startP.dec.b]
            result = minimize(ModelNBSolver._minimize_func,
                args={'par': par},
                x0=start_vec, method='Nelder-Mead')
            sol = ModelNBSolver.getSolution(par, float(result.x[0]), float(result.x[1]))
            return sol
        return startP

    @staticmethod
    def _minimize_func(x, args):
        wn, b  = float(x[0]), float(x[1])
        par = args['par']
        man_profit, data = ModelNBSolver.profit(par, wn, b)
        if man_profit is None:
            return sys.maxsize
        return -man_profit
        
    @staticmethod
    def _minimize_func_but_set_b_to_zero(x, args):
        wn, b  = float(x[0]), 0
        par = args['par']
        man_profit, data = ModelNBSolver.profit(par, wn, b)
        if man_profit is None:
            return sys.maxsize
        return -man_profit

    @staticmethod
    def getSolution(par, wn, b):
        man_profit, sol_tuple = ModelNBSolver.profit(par, wn, b)
        if sol_tuple is None: return None
        if sol_tuple is None: return None
        dec = DecisionVariables(MODEL_NB, pn=sol_tuple[2], b=sol_tuple[1],
            wn=sol_tuple[0], rho=sol_tuple[3], qn=sol_tuple[4])
        case = _CASE_ONE if sol_tuple[3] > 1 else _CASE_TWO
        sol = Solution(dec, sol_tuple[5], sol_tuple[6], case)
        return sol


    @staticmethod
    def _retailer_decision(par, wn, b):
        a = par.a; s = par.s; tau = par.tau
        valid_pn_rhos = utils.model_nb_pn_rho(wn, b, a, s, tau)
        valids = []
        for pn, rho, case_ in valid_pn_rhos:
            qn = 1 - pn
            ret_profit = (1-pn) *(pn *(1 - tau/rho) - wn + (tau/rho)*b) - a *(rho - 1)
            if (rho >= 1 and ret_profit >= 0):
                valids.append([pn, rho, qn, ret_profit])
        if len(valids) > 0:
            return max(valids, key=lambda k: k[3])
        else:
            return None

    @staticmethod
    def profit(par, wn, b, case=_UNDEFINED_CASE):
        if b > wn or b < 0: return None, None
        ret_dec = ModelNBSolver._retailer_decision(par, wn, b)
        if ret_dec is None: return None, None
        pn, rho, qn, ret_prof = ret_dec
        man_profit = qn*(wn - par.cn - (par.tau/rho)*b + (par.tau/rho)*par.s)
        return man_profit, (wn, b, pn, rho, qn, man_profit, ret_prof)

class ModelNBGridSearch:
    @staticmethod
    def search(par, resolution='high', b_range=[0,1]):
        if par.a == 0: return None
        wn_range = [0, 1]
        raster_size = 20
        iter = 5 #high auf 10?
        if resolution == 'low':
            raster_size = 10
            iter = 1
        elif resolution == 'very high':
            raster_size = 30
            iter = 10
        elif resolution == 'middle':
            raster_size = 50
            iter = 3
        
        grid_search_func = ModelNBSolver.profit
        sol_tuple = ImprovedGridSearch2D.maximize(grid_search_func,
            par, wn_range, b_range, raster_size, iter)
        if sol_tuple is None: return None
        dec = DecisionVariables(MODEL_NB, pn=sol_tuple[2], b=sol_tuple[1],
            wn=sol_tuple[0], rho=sol_tuple[3], qn=sol_tuple[4])
        case = _CASE_ONE if sol_tuple[3] > 1 else _CASE_TWO
        sol = Solution(dec, sol_tuple[5], sol_tuple[6], case)
        return sol
        

  
class GridSearch2D:
    """
    This class offers a grid search maximizing a
    for 2 dimensional concave function
    
    The function to optimize `func` must accept 3 arguments:
        (func_arg, x, y) - where func_arg is passed through
    The function must return 2 values:
        the first is a numerical object that supports the > operator
        the second is an object that will be finally be returned in case of a maximium
        
    The function must return None, None if there is no solution for a given x,y
    """
    
    @staticmethod
    def maximize(func, func_arg, x_start_range, y_start_range, raster_size, iterations, raster_start_size=None):
        raster_start_size = raster_size if raster_start_size is None else raster_start_size
        x_lim, y_lim = x_start_range, y_start_range
        x_range, y_range = x_start_range[:], y_start_range[:]
        best_f_val, best_f_obj = None, None
        best_x, best_y = None, None
        for iter in range(iterations):
            raster_it_size = raster_size if iter >= 1 else raster_start_size
            f_val, f_obj, f_x, f_y = GridSearch2D._search_raster(
                        func, func_arg, x_range, y_range, raster_it_size)
            if f_val is None: return best_f_obj
            
            # update best values if found
            if best_f_val is None or f_val > best_f_val:
                best_f_val, best_f_obj, best_x, best_y = f_val, f_obj, f_x, f_y
                
            # update new search ranges
            x_range = GridSearch2D._range(best_x, x_range, raster_it_size, x_lim[0], x_lim[1])
            y_range = GridSearch2D._range(best_y, y_range, raster_it_size, y_lim[0], y_lim[1])
        return best_f_obj
        
    @staticmethod
    def _range(point, old_range, raster_size, limit_low, limit_up):
        """ Returns a new 1-dimensional search range for deeper search
            
        Args:
            point (float): The (old) best solution was found at this point
            old_range ((float, float)): A tuple that describes the old range
            raster_size (int): The amount of search points at the (old) range
            limit_low (float): This is the lower limit of our global search range
            limit_up (float): Upper limit of global search range
            
        Returns:
            [new_low, new_up]
        """
        # new distance must at least contain the old neighbour of `point`
        old_distance_of_one_raster = (old_range[1] - old_range[0]) / raster_size
        new_distance = old_distance_of_one_raster * 2
        # we must not exceed global limits
        new_low = max(point - new_distance/2, limit_low)
        new_up = min(point + new_distance/2, limit_up)
        new_range = [new_low, new_up]
        
        return new_range
        
    @staticmethod
    def _search_raster(func, func_arg, x_range, y_range, raster_size):
        """
            this helper method creates the raster, calls the func and returns
            both returned best values and corresponding x and y val
        """
        best_f_val, best_f_obj, best_x, best_y = None, None, None, None
        y_num = raster_size if y_range[0] != y_range[1] else 1
        for x in np.linspace(x_range[0], x_range[1], num=raster_size):
            for y in np.linspace(y_range[0], y_range[1], num=y_num):
                x, y = float(x), float(y)
                # call our function
                f_val, f_obj = func(func_arg, x, y)
                if f_val == None: continue
                if best_f_val is None or f_val > best_f_val:
                    best_f_val, best_f_obj = f_val, f_obj
                    best_x, best_y = x, y
        return best_f_val, best_f_obj, best_x, best_y
        
class ImprovedGridSearch2D:
    """
    This class offers a grid search maximizing a
    for 2 dimensional concave function
    
    The function to optimize `func` must accept 3 arguments:
        (func_arg, x, y) - where func_arg is passed through
    The function must return 2 values:
        the first is a numerical object that supports the > operator
        the second is an object that will be finally be returned in case of a maximium
        
    The function must return None, None if there is no solution for a given x,y
    """
    
    @staticmethod
    def maximize(func, func_arg, x_start_range, y_start_range, raster_size, iterations, raster_start_size=None, caller_info=[]):
        raster_start_size = raster_size if raster_start_size is None else raster_start_size
        x_lim, y_lim = x_start_range, y_start_range
        x_range, y_range = x_start_range[:], y_start_range[:]
        best_f_val, best_f_obj = None, None
        best_x, best_y = None, None
        for iter in range(iterations):
            raster_it_size = raster_size if iter >= 1 else raster_start_size
            f_val, f_obj, f_x, f_y = ImprovedGridSearch2D._search_raster(
                        func, func_arg, x_range, y_range, raster_it_size)
            if f_val is None: return best_f_obj
            # update best values if found
            if best_f_val is None or f_val > best_f_val:
                best_f_val, best_f_obj, best_x, best_y = f_val, f_obj, f_x, f_y
                
            # update new search ranges
            x_range, y_range = ImprovedGridSearch2D._range((best_x, best_y), (x_range, y_range), raster_it_size, (x_lim[0], y_lim[0]), (x_lim[1], y_lim[1]), func_arg)
            if x_range[0] == x_range[1] and y_range[0] == y_range[1]:
                break
            
            # log some info for debugging purpose
            caller_info.append({'best_x': best_x, 'best_y' : best_y})
        return best_f_obj
        
    @staticmethod
    def _range(point, old_range, raster_size, limit_low, limit_up, func_arg):
        """ Returns a new 1-dimensional search range for deeper search
            
        Args:
            point (float, float): The (old) best solution was found at this point (x,y)
            old_range ((float, float), (float, float)): A tuple that describes the old range ((x_low, x_up), (y_low, y_up))
            raster_size (int): The amount of search points at the (old) range
            limit_low (float, float): This is the lower limit of our global search range
            limit_up (float, float): Upper limit of global search range (x,y)
            
        Returns:
            [[new_x_low, new_x_up], [new_y_low, new_y_up]]
        """
        assert old_range[0][0] <= point[0] <= old_range[0][1]
        # new distance must at least contain the old neighbour of `point`
        old_distance_of_one_raster = (old_range[0][1] - old_range[0][0]) / (raster_size-1)
        # are we at the edge? if yes - do not zoom in!
        if point[0] == old_range[0][0] or point[0] == old_range[0][1] or point[1] == old_range[1][0] or point[1] == old_range[1][1]:
            new_low_x = max(point[0] - old_distance_of_one_raster*(raster_size-1), limit_low[0])
            new_up_x = min(point[0] + old_distance_of_one_raster*(raster_size-1), limit_up[0])
            new_low_y = max(point[1] - old_distance_of_one_raster*(raster_size-1), limit_low[1])
            new_up_y = min(point[1] + old_distance_of_one_raster*(raster_size-1), limit_up[1])
        else:
            # zoom in
            new_distance = old_distance_of_one_raster * 4
            # we must not exceed global limits
            new_low_x = max(point[0] - new_distance/2, limit_low[0])
            new_up_x = min(point[0] + new_distance/2, limit_up[0])
            new_low_y = max(point[1] - new_distance/2, limit_low[1])
            new_up_y = min(point[1] + new_distance/2, limit_up[1])
        new_range = [(new_low_x, new_up_x), (new_low_y, new_up_y)]
        return new_range
        
    @staticmethod
    def _search_raster(func, func_arg, x_range, y_range, raster_size):
        """
            this helper method creates the raster, calls the func and returns
            four values: (i) best fval (ii) corresp. f_obj (iii) best_x, (iv) best_y
        """
        best_f_val, best_f_obj, best_x, best_y = None, None, None, None
        y_num = raster_size if y_range[0] != y_range[1] else 1
        for x in np.linspace(x_range[0], x_range[1], num=raster_size):
            for y in np.linspace(y_range[0], y_range[1], num=y_num):
                x, y = float(x), float(y)
                # call our function
                f_val, f_obj = func(func_arg, x, y)
                if f_val == None: continue
                if best_f_val is None or f_val > best_f_val:
                    best_f_val, best_f_obj = f_val, f_obj
                    best_x, best_y = x, y
        return best_f_val, best_f_obj, best_x, best_y
        
class ModelTwoSolver:
    solver = ModelTwoNumericalSolver()
    
    @staticmethod
    def solve(par):
        return ModelTwoSolver.solver.optimize(par)

class SolverProxy:
    def __init__(self):
        self.db = Database()
        self.model_1_solver = ModelOneNumericalSolver()
        self.model_2_solver = ModelTwoNumericalSolver()
        self.model_1_quad_search = ModelOneQuadGridSearch()
        self.model_2_quad_search = ModelTwoQuadGridSearch()
        self.model_nb_search = ModelNBGridSearch()
        #self.model_2_test = ModelTwoGridSearch()
        
    def read_calculation(self, par):
        """ Reads a calculation from database. Raises a CalculationNotFoundError if not found. """
        state, sol = db.read_calculation(par)
        if state == Database.NOT_IN_DB:
            raise CalculationNotFoundError()
        return sol
        
    def read_or_calc_and_write(self, par, comment=None, resolution='high'):
        """ doesn't commit after write! """
        state, sol = self.db.read_calculation(par)
        if state == Database.NOT_IN_DB:
            sol = self.calculate(par, resolution)
            self.db.write_calculation(par, sol, comment)
        return sol
        
    def commit(self):
        self.db.commit()
    
    def calculate(self, par, resolution='high'):
        """ Tries to find a solution for problem of type par"""
        if par.model == MODEL_1:
            sol = self.model_1_solver.optimize(par)
        elif par.model == MODEL_2:
            sol = self.model_2_solver.optimize(par)
        elif par.model == MODEL_2_QUAD:
            sol = ModelTwoQuadSolver.solve(par, resolution)
        elif par.model == MODEL_1_QUAD:
            sol = self.model_1_quad_search.search(par, resolution)
        elif par.model == MODEL_NB:
            sol = ModelNBSolver.solve(par, resolution)
        return sol
        
    def beginWrite(self):
        self.db.beginWrite()
        
    def endWrite(self):
        self.db.endWrite()
        
class CalculationNotFoundError (RuntimeError):
    pass
    
class Database:
    def __init__(self):
        self.conn = sqlite3.connect('database.db')
        # turn on foreign key support
        #self.conn.execute('''PRAGMA foreign_keys = ON;''')
        # create tables
        with self.conn as c:
            c.execute('''CREATE TABLE IF NOT EXISTS calculation
                ( calc_id integer primary key,
                  tau real, a real, s real, cr real, cn real, delta real,
                  wn real, pr real, pn real, rho real, qn real, qr real,
                  profit_man real, profit_ret real,
                  sol_case text, model integer, comment text, lastmodified text)''')
    
    instance = None
    FOUND = 0
    NO_SOLUTION = 1
    NOT_IN_DB = 2
        
    @staticmethod
    def getInstance():
        if Database.instance is None:
            Database.instance = Database()
        return Database.instance
        
    def beginWrite(self):
        self.conn.execute('BEGIN')
        
    def endWrite(self):
        self.conn.execute('END')
        
    def read_calculation(self, par):
        """ Returns two values: (state and solution)"""
        tau, a, s, cr, cn, delta = par.tau, par.a, par.s, par.cr, par.cn, par.delta
        model = par.model
        query_str = '''
            SELECT tau, a, s, cr, cn, delta, 
                    wn, pr, pn, rho, qn, qr,
                    profit_man, profit_ret,
                    sol_case, model, comment
            FROM calculation
            WHERE tau is ? and a is ? and s is ? and cr is ? and cn is ? and delta is ? and model is ?'''
        cur = self.conn.execute(query_str,
            (tau, a, s, cr, cn, delta, model))
        row = cur.fetchone()
        if row is None:
            # nothing found - never calculated
            return Database.NOT_IN_DB, None
        wn, pr, pn, rho, qn, qr, profit_man, profit_ret, sol_case, model, comment = (row[i] for i in range(6, 17))
        if profit_man is None:
            # found but the solution is None
            return Database.NO_SOLUTION, None
        if model == MODEL_NB:
            dec = DecisionVariables(model=model, pn=pn, pr=pr, wn=wn, rho=rho, qn=qn, b=qr)
        else:
            dec = DecisionVariables(model=model, pn=pn, pr=pr, wn=wn, rho=rho, qn=qn, qr=qr)
        sol = Solution(dec, profit_man, profit_ret, sol_case)
        return Database.FOUND, sol
        
    def write_calculation(self, par, sol, comment=None):
        tau = par.tau
        a = par.a
        s = par.s
        cr = par.cr
        cn = par.cn
        delta = par.delta
        model = par.model
        if sol == None:
            wn, pr, pn, rho, qn, qr, profit_man, profit_ret, sol_case = None, None, None, None, None, None, None, None, None
        else:
            wn = sol.dec.wn
            pr = sol.dec.pr
            pn = sol.dec.pn
            rho = sol.dec.rho
            qn = sol.dec.qn
            if model == MODEL_NB:
                qr = sol.dec.b
            else:
                qr = sol.dec.qr
            profit_man = sol.profit_man
            profit_ret = sol.profit_ret
            sol_case = sol.case
        cur = self.conn.execute(
            '''INSERT INTO calculation (tau, a, s, cr, cn, delta, 
                    wn, pr, pn, rho, qn, qr,
                    profit_man, profit_ret,
                    sol_case, model, comment, lastmodified)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, datetime('now'))''',
               (tau, a, s, cr, cn, delta, wn, pr, pn,
                rho, qn, qr, profit_man, profit_ret, sol_case, model, comment))
    
    def commit(self):
        self.conn.commit()
        
    def delete_where_comment(self, comment):
        self.conn.execute('''
            DELETE FROM calculation
            WHERE comment = ?''', (comment, ))
        self.commit()

        
if __name__ == '__main__':
    par = Parameter(MODEL_2_QUAD, tau=0.09, a=.0008163265306122449, s=0.04000000000000001, cn=0.1, cr=0.04000000000000001, delta=0.7956)
    blub = ModelTwoGridTester(par)
    blub.plot()
    #blub.plot_profit_surface()
      