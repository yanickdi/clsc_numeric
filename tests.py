import unittest
from math import sqrt

from solver import ModelOneNumericalSolver, is_prof_pos, Parameter, MODEL_1, MODEL_2
from generator import Generator, MemoryOutputFile

class TestModelOneNumericalSolver(unittest.TestCase):
    def test_case_1a(self):
        solver = ModelOneNumericalSolver()
        # this parameter should lead to case one (roh is gte 1)
        par = Parameter(MODEL_1, tau=0.1, a=0.005, s=0.0005, cn=0.01)
        # self checking the my test input variables..
        self.assertTrue(self.__input_is_in_case_1(par))
        
        analyitcal_profit_manufacturer = ((1-par.cn)**2 / 8) - (1/2)*(1+par.cn-2*par.s) * (par.tau*par.a)**(1/2) + (1/2)*par.a*par.tau
        profit_solver_manufacturer,_ = solver.calc_profits(par, solver.optimize(par))
        self.assertAlmostEqual(analyitcal_profit_manufacturer, profit_solver_manufacturer)
        
    def test_case_1a_dec_vars(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case one (roh is gte 1)
        par = Parameter(MODEL_1, tau=0.1, a=0.005, s=0.0005, cn=0.01)
        # self checking the test input variables..
        # if the following condition is true, it must lead to a case a optimization
        self.assertTrue(self.__input_is_in_case_1(par) and not self.__input_is_in_case_2(par))
        
        solver_dec_vars = solver.optimize(par)
        self.assertAlmostEqual(solver_dec_vars['pn'], (3+par.cn)/4 - (1/2)*(par.a*par.tau)**(1/2), msg='pn not the same')
        self.assertAlmostEqual(solver_dec_vars['wn'], (1+par.cn)/2 - (par.a*par.tau)**(1/2), msg='wn not the same')
        self.assertAlmostEqual(solver_dec_vars['roh'], par.tau/2 + (1-par.cn)/4 * (par.tau/par.a)**(1/2), msg='roh not the same')

    def test_case_2a(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case two (roh is equal to 1)
        par = Parameter(MODEL_1, tau=0.1, a=0.006, s=0.005, cn=0.3)
        # self checking the test input variables..
        # if the following condition is true, it must lead to a case b optimization
        self.assertTrue(self.__input_is_in_case_2(par) and not self.__input_is_in_case_1(par))
        
        analyitcal_profit_manufacturer = ((1-par.cn-par.tau+par.s*par.tau)**2)/(8*(1-par.tau))
        profit_solver_manufacturer,_ = solver.calc_profits(par, solver.optimize(par))
        self.assertAlmostEqual(analyitcal_profit_manufacturer, profit_solver_manufacturer)
        
        
    def test_case_1_or_2(self):
        solver = ModelOneNumericalSolver()
        par = Parameter(MODEL_1, tau=0.3, a=0.01, s=0, cn=0.3)
        # self checking if input vars not in case 1 and not in case 2:
        self.assertTrue(self.__input_is_in_case_1(par) and self.__input_is_in_case_2(par))
        prof_solver_man, prof_solver_ret = solver.calc_profits(par, solver.optimize(par))
        self.assertAlmostEqual(prof_solver_ret, 0.00428571) # would be case 2 solution, because is higher than case 1 solution
        self.assertAlmostEqual(prof_solver_man, 0.02857143) 
        
    def test_case_2_dec_vars(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case two (roh is equal to 1)
        par = Parameter(MODEL_1, tau=0.1, a=0.006, s=0.005, cn=0.3)
        # self checking the test input variables
        self.assertTrue(self.__input_is_in_case_2(par))
        
        solver_dec_vars = solver.optimize(par)
        self.assertAlmostEqual(solver_dec_vars['pn'], 0.83319444, msg='pn not the same')
        self.assertAlmostEqual(solver_dec_vars['wn'], (1/(1-par.tau)) * ((1+par.cn)/2 - (par.tau*(1+par.s))/2), msg='wn not the same')
        self.assertAlmostEqual(solver_dec_vars['roh'], 1.0, msg='roh not the same')
    
    def test_qn(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case one (roh is gte 1)
        par = Parameter(MODEL_1, tau=0.1, a=0.005, s=0.0005, cn=0.01)
        dec_vars = solver.optimize(par)
        self.assertAlmostEqual(dec_vars['qn'], 1 - dec_vars['pn'])
        # TODO: also test a case leading to roh == 1
        
    def test_sol_not_possible(self):
        solver = ModelOneNumericalSolver()
        par = Parameter(MODEL_1, tau=1, a=0.01, s=0, cn=0.8)
        self.assertIsNone(solver.optimize(par))
        
    def __input_is_in_case_1(self, par):
        return par.cn <= 1 - 4*(1-par.tau/2)*sqrt(par.a/par.tau)
        
    def __input_is_in_case_2(self, par):
        return par.cn >= 1 - par.tau*(1-par.s) - 4*(1-par.tau) * sqrt(par.a/par.tau)
        
class TestGenerator(unittest.TestCase):
    def test_model_1_compare_analytical(self):
        mof = MemoryOutputFile()
        generator = Generator(MODEL_1, mof)
        solver = ModelOneNumericalSolver()
        ana_solver = AnalyticalSolver()
        
        generator.generate()        
        for solution in mof.getSolutions():
            par, solver_dec_vars = solution['par'], solution['dec_vars']
            solver_prof_man, solver_prof_ret = solution['profit_man'], solution['profit_ret']
            assert par != None
            dec_vars, prof_man, prof_ret = ana_solver.calcModelOne(par)
            if dec_vars == None:
                self.assertIsNone(prof_man)
                self.assertIsNone(prof_ret)
                self.assertIsNone(solver_dec_vars)
                self.assertIsNone(solver_prof_man)
                self.assertIsNone(solver_prof_ret)
            elif solver_dec_vars == None:
                self.assertIsNone(solver_prof_man)
                self.assertIsNone(solver_prof_ret)
                if dec_vars != None:
                    print(const_args)
                self.assertIsNone(dec_vars)
                self.assertIsNone(prof_man)
                self.assertIsNone(prof_ret)
            else:
                if round(solver_dec_vars['pn'], 7) != round(dec_vars['pn'], 7):
                    print(solver_dec_vars)
                    print(dec_vars)
                    print(const_args)
                    
                self.assertAlmostEqual(solver_dec_vars['pn'], dec_vars['pn'])
                self.assertAlmostEqual(solver_dec_vars['wn'], dec_vars['wn'])
                self.assertAlmostEqual(solver_dec_vars['roh'], dec_vars['roh'])
                self.assertAlmostEqual(solver_dec_vars['qn'], dec_vars['qn'])
                self.assertAlmostEqual(solver_prof_man, prof_man)
                self.assertAlmostEqual(solver_prof_ret, prof_ret)
        

class AnalyticalSolver:
    """
        This class is used to try to give analytical solutions of model input data.
        It is used to test the Solver's Solution and should only be used by UnitTests
    """
    
    def calcModelOne(self, par):
        """
        Returns a tuple (dec_vars, profit_manufacturer, profit_retailer)
        
        If the solution isnt possible, this method will return a tuple of (None, None, None)
        """
        
        case_1_pn  = (3+par.cn)/4 - (1/2)*(par.a*par.tau)**(1/2)
        case_1_wn  = (1+par.cn)/2 - (par.a*par.tau)**(1/2)
        case_1_roh = par.tau/2 + (1-par.cn)/4 * (par.tau/par.a)**(1/2)
        case_1_qn  = (1-par.cn)/4 + (1/2)*(par.a*par.tau)**(1/2)
        if case_1_roh != 0:
            # i can skip this solution..
            case_1_prof_man = case_1_qn * ( case_1_wn * (1- par.tau/case_1_roh) - par.cn + (par.tau/case_1_roh)*par.s)
            case_1_prof_ret = case_1_qn * (case_1_pn - case_1_wn)*(1- par.tau/case_1_roh) - par.a*case_1_roh
        
        if par.tau != 1:
            case_2_pn  = (1/(1-par.tau)) * ( (3+par.cn)/(4) - (par.tau*(3+par.s)/(4)))
            case_2_wn  = (1/(1-par.tau)) * ((1+par.cn)/(2) - (par.tau*(1+par.s))/(2) )
            case_2_roh = 1
            case_2_qn  = (1/(1-par.tau)) * ( (1-par.cn)/(4) - (par.tau*(1-par.s))/(4) )
            case_2_prof_man = case_2_qn * ( case_2_wn * (1- par.tau/case_2_roh) - par.cn + (par.tau/case_2_roh)*par.s)
            case_2_prof_ret = case_2_qn * (case_2_pn - case_2_wn)*(1- par.tau/case_2_roh) - par.a*case_2_roh
        else:
            # par.tau == 1 leads to division by zero
            pass
        
        if round(case_1_roh, 7) >= 1:
            # i can take both solutions
            if par.tau != 1 and case_2_prof_man > case_1_prof_man and case_2_prof_man >= 0 and case_2_prof_ret >= 0:
                sol = 'CASE_2'
            else:
                sol = 'CASE_1'
        else:
            # only case two would be possible
            if par.tau == 1:
                # no analytical solution possible:
                #raise RuntimeError('no analytical solution possible')
                return (None, None, None)
            # have to fall back on case 2
            else:
                # only case 2 analytical is possible:
                ret_val = ({'pn' : case_2_pn, 'wn' : case_2_wn, 'roh' : case_2_roh, 'qn' : case_2_qn}, case_2_prof_man, case_2_prof_ret)
                if is_prof_pos(ret_val[1]) and is_prof_pos(ret_val[2]): return ret_val
                else: return (None, None, None)
        
        if sol == 'CASE_1':
            ret_val = ({'pn' : case_1_pn, 'wn' : case_1_wn, 'roh' : case_1_roh, 'qn' : case_1_qn}, case_1_prof_man, case_1_prof_ret)
        else:
            ret_val = ({'pn' : case_2_pn, 'wn' : case_2_wn, 'roh' : case_2_roh, 'qn' : case_2_qn}, case_2_prof_man, case_2_prof_ret)
            
        if not is_prof_pos(ret_val[2]) or not is_prof_pos(ret_val[1]):
            if ret_val[2] >= 0 or ret_val[1] >= 0:
                pass
            return (None, None, None)
            
        return ret_val
        
if __name__ == '__main__':
    unittest.main()