import unittest
from math import sqrt

from solver import build_args, check_args, ModelOneNumericalSolver
from generator import Generator, MemoryOutputFile, MODEL_1, MODEL_2

class TestModelOneNumericalSolver(unittest.TestCase):
    def test_case_a(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case one (roh is gte 1)
        const_args = build_args(tau=0.1, a=0.005, s=0.0005, cn=0.01)
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        # self checking the my test input variables..
        self.assertTrue(cn <= 1 - 4*(1-tau/2)*(a/tau)**(1/2))
        
        analyitcal_profit_manufacturer = ((1-cn)**2 / 8) - (1/2)*(1+cn-2*s) * (tau*a)**(1/2) + (1/2)*a*tau
        profit_solver_manufacturer,_ = solver.calc_profits(const_args, solver.optimize(const_args))
        self.assertAlmostEqual(analyitcal_profit_manufacturer, profit_solver_manufacturer)
        
    def test_case_a_dec_vars(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case one (roh is gte 1)
        const_args = build_args(tau=0.1, a=0.005, s=0.0005, cn=0.01)
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        # self checking the test input variables..
        # if the following condition is true, it must lead to a case a optimization
        self.assertTrue(cn <= 1 - 4*(1-tau/2)*(a/tau)**(1/2))
        
        solver_dec_vars = solver.optimize(const_args)
        self.assertAlmostEqual(solver_dec_vars['pn'], (3+cn)/4 - (1/2)*(a*tau)**(1/2), msg='pn not the same')
        self.assertAlmostEqual(solver_dec_vars['wn'], (1+cn)/2 - (a*tau)**(1/2), msg='wn not the same')
        self.assertAlmostEqual(solver_dec_vars['roh'], tau/2 + (1-cn)/4 * (tau/a)**(1/2), msg='roh not the same')

    def test_case_b(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case two (roh is equal to 1)
        const_args = build_args(tau=0.1, a=0.006, s=0.005, cn=0.3)
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        # self checking the test input variables..
        # if the following condition is true, it must lead to a case b optimization
        self.assertTrue(cn > 1 - 4*(1-tau/2)*(a/tau)**(1/2))
        
        analyitcal_profit_manufacturer = ((1-cn-tau+s*tau)**2)/(8*(1-tau))
        profit_solver_manufacturer,_ = solver.calc_profits(const_args, solver.optimize(const_args))
        self.assertAlmostEqual(analyitcal_profit_manufacturer, profit_solver_manufacturer)
        
    def test_case_b_dec_vars(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case two (roh is equal to 1)
        const_args = build_args(tau=0.1, a=0.006, s=0.005, cn=0.3)
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        # self checking the test input variables..
        # if the following condition is true, it must lead to a case b optimization
        self.assertTrue(cn > 1 - 4*(1-tau/2)*(a/tau)**(1/2))
        
        solver_dec_vars = solver.optimize(const_args)
        self.assertAlmostEqual(solver_dec_vars['pn'], (1/(1-tau)) * ((3+cn)/4 - (tau*(3+s))/4), msg='pn not the same')
        self.assertAlmostEqual(solver_dec_vars['wn'], (1/(1-tau)) * ((1+cn)/2 - (tau*(1+s))/2), msg='wn not the same')
        self.assertAlmostEqual(solver_dec_vars['roh'], 1.0, msg='roh not the same')
    
    def test_qn(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case one (roh is gte 1)
        const_args = build_args(tau=0.1, a=0.005, s=0.0005, cn=0.01)
        dec_vars = solver.optimize(const_args)
        self.assertAlmostEqual(dec_vars['qn'], 1 - dec_vars['pn'])
        # TODO: also test a case leading to roh == 1
        
class TestGenerator(unittest.TestCase):
    def test_model_1_compare_analytical(self):
        mof = MemoryOutputFile()
        generator = Generator(MODEL_1, mof)
        solver = ModelOneNumericalSolver()
        ana_solver = AnalyticalSolver()
        
        generator.generate()        
        for solution in mof.getSolutions():
            const_args, solver_dec_vars = solution['const_args'], solution['dec_vars']
            solver_prof_man, solver_prof_ret = solution['profit_man'], solution['profit_ret']
        
            dec_vars, prof_man, prof_ret = ana_solver.calcModelOne(const_args)
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
    
    def calcModelOne(self, const_args):
        """ Returns a tuple (dec_vars, profit_manufacturer, profit_retailer) """
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        
        case_1_pn  = (3+cn)/4 - (1/2)*(a*tau)**(1/2)
        case_1_wn  = (1+cn)/2 - (a*tau)**(1/2)
        case_1_roh = tau/2 + (1-cn)/4 * (tau/a)**(1/2)
        case_1_qn  = (1-cn)/4 + (1/2)*(a*tau)**(1/2)
        if case_1_roh != 0:
            # i can skip this solution..
            case_1_prof_man = case_1_qn * ( case_1_wn * (1- tau/case_1_roh) - cn + (tau/case_1_roh)*s)
            case_1_prof_ret = case_1_qn * (case_1_pn - case_1_wn)*(1- tau/case_1_roh) - a*case_1_roh
        
        case_2_pn  = (1/(1-tau)) * ( (3+cn)/(4) - (tau*(3+s)/(4)))
        case_2_wn  = (1/(1-tau)) * ((1+cn)/(2) - (tau*(1+s))/(2) )
        case_2_roh = 1
        case_2_qn  = (1/(1-tau)) * ( (1-cn)/(4) - (tau*(1-s))/(4) )
        case_2_prof_man = case_2_qn * ( case_2_wn * (1- tau/case_2_roh) - cn + (tau/case_2_roh)*s)
        case_2_prof_ret = case_2_qn * (case_2_pn - case_2_wn)*(1- tau/case_2_roh) - a*case_2_roh
        
        if case_1_roh >= 1:
            # i can take both solutions
            if case_2_prof_man > case_1_prof_man:
                sol = 'CASE_2'
            else:
                sol = 'CASE_1'
        else:
            # have to fall back on case 2
            sol = 'CASE_2'
        
        if sol == 'CASE_1':
            return ({'pn' : case_1_pn, 'wn' : case_1_wn, 'roh' : case_1_roh, 'qn' : case_1_qn}, case_1_prof_man, case_1_prof_ret)
        else:
            return ({'pn' : case_2_pn, 'wn' : case_2_wn, 'roh' : case_2_roh, 'qn' : case_2_qn}, case_2_prof_man, case_2_prof_ret)
        
if __name__ == '__main__':
    unittest.main()