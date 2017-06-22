import unittest
from math import sqrt

from solver import ModelOneNumericalSolver, ModelTwoNumericalSolver, is_prof_pos, Parameter, DecisionVariables, MODEL_1, MODEL_2
from generator import Generator, MemoryOutputFile
from plotter import _rho_diff_func

class TestModelOneNumericalSolver(unittest.TestCase):
    def test_case_1a(self):
        solver = ModelOneNumericalSolver()
        # this parameter should lead to case one (rho is gte 1)
        par = Parameter(MODEL_1, tau=0.1, a=0.005, s=0.0005, cn=0.01)
        # self checking the my test input variables..
        self.assertTrue(self.__input_is_in_case_1(par))
        
        analyitcal_profit_manufacturer = ((1-par.cn)**2 / 8) - (1/2)*(1+par.cn-2*par.s) * (par.tau*par.a)**(1/2) + (1/2)*par.a*par.tau
        sol = solver.optimize(par)
        self.assertAlmostEqual(analyitcal_profit_manufacturer, sol.profit_man)
        
    def test_case_1a_dec_vars(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case one (rho is gte 1)
        par = Parameter(MODEL_1, tau=0.1, a=0.005, s=0.0005, cn=0.01)
        # self checking the test input variables..
        # if the following condition is true, it must lead to a case a optimization
        self.assertTrue(self.__input_is_in_case_1(par) and not self.__input_is_in_case_2(par))
        
        sol = solver.optimize(par)
        self.assertAlmostEqual(sol.dec.pn, (3+par.cn)/4 - (1/2)*(par.a*par.tau)**(1/2), msg='pn not the same')
        self.assertAlmostEqual(sol.dec.wn, (1+par.cn)/2 - (par.a*par.tau)**(1/2), msg='wn not the same')
        self.assertAlmostEqual(sol.dec.rho, par.tau/2 + (1-par.cn)/4 * (par.tau/par.a)**(1/2), msg='rho not the same')

    def test_case_2a(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case two (rho is equal to 1)
        par = Parameter(MODEL_1, tau=0.1, a=0.006, s=0.005, cn=0.3)
        # self checking the test input variables..
        # if the following condition is true, it must lead to a case b optimization
        self.assertTrue(self.__input_is_in_case_2(par) and not self.__input_is_in_case_1(par))
        
        analyitcal_profit_manufacturer = ((1-par.cn-par.tau+par.s*par.tau)**2)/(8*(1-par.tau))
        sol = solver.optimize(par)
        self.assertAlmostEqual(analyitcal_profit_manufacturer, sol.profit_man)
        
        
    def test_case_1_or_2(self):
        solver = ModelOneNumericalSolver()
        par = Parameter(MODEL_1, tau=0.3, a=0.01, s=0, cn=0.3)
        # self checking if input vars not in case 1 and not in case 2:
        self.assertTrue(self.__input_is_in_case_1(par) and self.__input_is_in_case_2(par))
        sol = solver.optimize(par)
        self.assertAlmostEqual(sol.profit_ret, 0.00428571) # would be case 2 solution, because is higher than case 1 solution
        self.assertAlmostEqual(sol.profit_man, 0.02857143) 
        
    def test_case_2_dec_vars(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case two (rho is equal to 1)
        par = Parameter(MODEL_1, tau=0.1, a=0.006, s=0.005, cn=0.3)
        # self checking the test input variables
        self.assertTrue(self.__input_is_in_case_2(par))
        
        sol = solver.optimize(par)
        self.assertAlmostEqual(sol.dec.pn, 0.83319444, msg='pn not the same')
        self.assertAlmostEqual(sol.dec.wn, (1/(1-par.tau)) * ((1+par.cn)/2 - (par.tau*(1+par.s))/2), msg='wn not the same')
        self.assertAlmostEqual(sol.dec.rho, 1.0, msg='rho not the same')
    
    def test_qn(self):
        solver = ModelOneNumericalSolver()
        # this args should lead to case one (rho is gte 1)
        par = Parameter(MODEL_1, tau=0.1, a=0.005, s=0.0005, cn=0.01)
        sol = solver.optimize(par)
        self.assertAlmostEqual(sol.dec.qn, 1 - sol.dec.pn)
        # TODO: also test a case leading to rho == 1
        
    def test_sol_not_possible(self):
        solver = ModelOneNumericalSolver()
        par = Parameter(MODEL_1, tau=1, a=0.01, s=0, cn=0.8)
        self.assertIsNone(solver.optimize(par))
        
    def __input_is_in_case_1(self, par):
        return par.cn <= 1 - 4*(1-par.tau/2)*sqrt(par.a/par.tau)
        
    def __input_is_in_case_2(self, par):
        return par.cn >= 1 - par.tau*(1-par.s) - 4*(1-par.tau) * sqrt(par.a/par.tau)

class TestModelTwoNumericalSolver(unittest.TestCase):
    def test_case_one_a(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to case one (a), but i will check the output anyway
        par = Parameter(MODEL_2, tau=.5, a=.01, s=.1, cr=.2, cn=.3, delta=.4)
        # TODO, check the case
        #self.assertTrue(self.__input_is_in_case_1(par))
        dec = solver._optimize_case_one_a(par)
        profit_man,profit_ret = solver.calc_profits(par, dec)
        
        self.assertAlmostEqual(dec.pn, .735606453541)
        self.assertAlmostEqual(dec.pr, .294242581416)
        self.assertAlmostEqual(dec.wn, .576970325666)
        self.assertAlmostEqual(dec.rho, 1.448143094543)
        self.assertAlmostEqual(dec.qr, 0)
        self.assertAlmostEqual(dec.qn, 0.264393546459)
        self.assertAlmostEqual(profit_man, 0.0296879322287)
        self.assertAlmostEqual(profit_ret, 0.0129795065546)
        self.assertAlmostEqual(dec.lambda1, 2.07500000000)
        self.assertAlmostEqual(dec.lambda2, 0)
        
    def test_case_one_b(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to case one (b), but i will check the output anyway
        par = Parameter(MODEL_2, tau=.1, a=.05, s=.1, cr=.2, cn=.3, delta=.8)
        # TODO, check the case
        dec = solver._optimize_case_one_b(par)
        profit_man, profit_ret = solver.calc_profits(par, dec)
        
        self.assertAlmostEqual(dec.pn, .62094305850)
        self.assertAlmostEqual(dec.pr, 0.48675444680)
        self.assertAlmostEqual(dec.wn, .55513167019)
        self.assertAlmostEqual(dec.rho, .20811388301)
        self.assertAlmostEqual(dec.qn, .32905694150)
        self.assertAlmostEqual(dec.qr, 0.06250000000)
        self.assertAlmostEqual(profit_man, .0236623643454)
        self.assertAlmostEqual(profit_ret, .0008443058496)
        self.assertAlmostEqual(dec.lambda1, 0)
        self.assertAlmostEqual(dec.lambda2, 0)
        
    def test_case_one_c(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to case one (b), but i will check the output anyway
        par = Parameter(MODEL_2, tau=.1, a=.05, s=.1, cr=.2, cn=.3, delta=.8)
        # TODO, check the case
        dec = solver._optimize_case_one_c(par)
        profit_man, profit_ret = solver.calc_profits(par, dec)
        
        self.assertAlmostEqual(dec.pn, .6081945407613)
        self.assertAlmostEqual(dec.pr, .4612574113277)
        self.assertAlmostEqual(dec.wn, .5551316701949)
        self.assertAlmostEqual(dec.qn, .2653143528319)
        self.assertAlmostEqual(dec.qr, .1581138830084)
        self.assertAlmostEqual(profit_man,  .0212244937790)
        self.assertAlmostEqual(profit_ret, -.0027016118539)
        self.assertAlmostEqual(dec.lambda1, 0)
        self.assertAlmostEqual(dec.lambda2, -0.050994071)
      
    def test_case_two_a(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to two a, but i will check the output anyway
        par = Parameter(MODEL_2, tau=.1, a=.05, s=.1, cr=.2, cn=.3, delta=.8)
        # TODO, check the case
        dec = solver._optimize_case_two_a(par)
        profit_man, profit_ret = solver.calc_profits(par, dec)
        
        self.assertAlmostEqual(dec.wn, .661111111)
        self.assertAlmostEqual(dec.pr, .574074074)
        self.assertAlmostEqual(dec.pn, .717592593)
        self.assertAlmostEqual(dec.rho, 1)
        self.assertAlmostEqual(dec.qn, .282407407)
        self.assertAlmostEqual(dec.qr, .0)
        self.assertAlmostEqual(profit_man,  .0861342592593)
        self.assertAlmostEqual(profit_ret, -.0356442901235)
        self.assertAlmostEqual(dec.lambda1, -.070740741)
        self.assertAlmostEqual(dec.lambda2, 0)
        
    def test_case_two_b(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to two a, but i will check the output anyway
        par = Parameter(MODEL_2, tau=.1, a=.05, s=.1, cr=.2, cn=.3, delta=.8)
        # TODO, check the case
        dec = solver._optimize_case_two_b(par)
        profit_man, profit_ret = solver.calc_profits(par, dec)
        
        self.assertAlmostEqual(dec.pn, .704748603)
        self.assertAlmostEqual(dec.pr, .542458101)
        self.assertAlmostEqual(dec.wn, .667039106)
        self.assertAlmostEqual(dec.rho, 1)
        self.assertAlmostEqual(dec.qn, .188547486)
        self.assertAlmostEqual(dec.qr, .133379888)
        self.assertAlmostEqual(profit_man,  .0908519553073)
        self.assertAlmostEqual(profit_ret, -.0436009721919)
        self.assertAlmostEqual(dec.lambda1, 0)
        self.assertAlmostEqual(dec.lambda2, 0)
        
    def test_case_two_c(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to two a, but i will check the output anyway
        par = Parameter(MODEL_2, tau=.1, a=.05, s=.1, cr=.2, cn=.3, delta=.8)
        # TODO, check the case
        dec = solver._optimize_case_two_c(par)
        profit_man, profit_ret = solver.calc_profits(par, dec)
        
        self.assertAlmostEqual(dec.pn, .71258064516)
        self.assertAlmostEqual(dec.pr, .56580645161)
        self.assertAlmostEqual(dec.wn, .65935483871)
        self.assertAlmostEqual(dec.rho, 1)
        self.assertAlmostEqual(dec.qn, .26612903226)
        self.assertAlmostEqual(dec.qr, .02661290323)
        self.assertAlmostEqual(profit_man,  .0878225806452)
        self.assertAlmostEqual(profit_ret, -.0372515608741)
        self.assertAlmostEqual(dec.lambda1, 0)
        self.assertAlmostEqual(dec.lambda2, 0.052903225806451626)

    def test_optimize_instance_a(self):
        solver = ModelTwoNumericalSolver()
        par = Parameter(MODEL_2, tau=.1, a=.05, s=.1, cr=.2, cn=.3, delta=.8)
        sol = solver.optimize(par)
        self.assertIsNone(sol)
        
    def test_optimize_instance_b(self):
        solver = ModelTwoNumericalSolver()
        par = Parameter(MODEL_2, tau=.1, a=.01, s=.1, cr=.2, cn=.3, delta=.8)
        sol = solver.optimize(par)
        self.assertIsNotNone(sol)
        self.assertAlmostEqual(sol.profit_man, 0.0878225806452)
        
    def test_optimize_instance_c(self):
        solver = ModelTwoNumericalSolver()
        par = Parameter(MODEL_2, tau=.0, a=.01, s=.0, cr=.0, cn=.1, delta=.8)
        sol = solver.optimize(par)
        self.assertIsNotNone(sol)
        self.assertAlmostEqual(sol.profit_man, 0.1687500000000)
        self.assertAlmostEqual(sol.profit_ret, 0.0181250000000)
        
    def test_optimize_instance_d(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to two a, but i will check the output anyway
        par = Parameter(MODEL_2, tau=.4, a=.01, s=.0, cr=.0, cn=.0, delta=.9)
        sol = solver.optimize(par)
        self.assertIsNotNone(sol)
        self.assertAlmostEqual(sol.profit_man, 0.1382727272727)
        
    def test_optimize_instance_e(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to two a, but i will check the output anyway
        par = Parameter(MODEL_2, tau=.7, a=.01, s=.4, cr=.4, cn=.6, delta=.9)
        sol = solver.optimize(par)
        self.assertIsNone(sol)
        
    def test_optimize_instance_f(self):
        solver = ModelTwoNumericalSolver()
        par = Parameter(MODEL_2, tau=.9, a=.1, s=.2, cr=.4, cn=.8, delta=.4)
        sol = solver.optimize(par)
        self.assertIsNone(sol)
             
class TestGenerator(unittest.TestCase):
    def test_model_1_compare_analytical(self):
        mof = MemoryOutputFile()
        generator = Generator(MODEL_1, mof)
        solver = ModelOneNumericalSolver()
        ana_solver = AnalyticalSolver()
        
        generator.generate()        
        for solution in mof.getSolutions():
            par, sol = solution['par'], solution['sol']
            assert par != None
            dec_vars, prof_man, prof_ret = ana_solver.calcModelOne(par)
            if dec_vars == None:
                self.assertIsNone(prof_man)
                self.assertIsNone(prof_ret)
                self.assertIsNone(sol)
            elif sol == None:
                if dec_vars != None:
                    print(const_args)
                self.assertIsNone(dec_vars)
                self.assertIsNone(prof_man)
                self.assertIsNone(prof_ret)
            else:
                self.assertAlmostEqual(sol.dec.pn, dec_vars.pn)
                self.assertAlmostEqual(sol.dec.wn, dec_vars.wn)
                self.assertAlmostEqual(sol.dec.rho, dec_vars.rho)
                self.assertAlmostEqual(sol.dec.qn, dec_vars.qn)
                self.assertAlmostEqual(sol.profit_man, prof_man)
                self.assertAlmostEqual(sol.profit_ret, prof_ret)
        
class TestContour(unittest.TestCase):
    par1 = Parameter(MODEL_1, tau=0.3, a=0.013, s=0.1, cn=0.3)
    par2 = Parameter(MODEL_2, tau=0.3, a=0.013, s=0.1, cn=0.3, cr=0.1, delta=0.4)
    solver_1, solver_2 = ModelOneNumericalSolver(), ModelTwoNumericalSolver()
    sol_1, sol_2 = solver_1.optimize(par1), solver_2.optimize(par2)
    diff = _rho_diff_func(sol_1, sol_2)
    print(diff)
        
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
        case_1_rho = par.tau/2 + (1-par.cn)/4 * (par.tau/par.a)**(1/2)
        case_1_qn  = (1-par.cn)/4 + (1/2)*(par.a*par.tau)**(1/2)
        if case_1_rho != 0:
            # i can skip this solution..
            case_1_prof_man = case_1_qn * ( case_1_wn * (1- par.tau/case_1_rho) - par.cn + (par.tau/case_1_rho)*par.s)
            case_1_prof_ret = case_1_qn * (case_1_pn - case_1_wn)*(1- par.tau/case_1_rho) - par.a*case_1_rho
        
        if par.tau != 1:
            case_2_pn  = (1/(1-par.tau)) * ( (3+par.cn)/(4) - (par.tau*(3+par.s)/(4)))
            case_2_wn  = (1/(1-par.tau)) * ((1+par.cn)/(2) - (par.tau*(1+par.s))/(2) )
            case_2_rho = 1
            case_2_qn  = (1/(1-par.tau)) * ( (1-par.cn)/(4) - (par.tau*(1-par.s))/(4) )
            case_2_prof_man = case_2_qn * ( case_2_wn * (1- par.tau/case_2_rho) - par.cn + (par.tau/case_2_rho)*par.s)
            case_2_prof_ret = case_2_qn * (case_2_pn - case_2_wn)*(1- par.tau/case_2_rho) - par.a*case_2_rho
        else:
            # par.tau == 1 leads to division by zero
            pass
        
        if round(case_1_rho, 7) >= 1:
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
                ret_val = (DecisionVariables(MODEL_1, pn=case_2_pn, wn=case_2_wn, rho=case_2_rho, qn=case_2_qn), case_2_prof_man, case_2_prof_ret)
                if is_prof_pos(ret_val[1]) and is_prof_pos(ret_val[2]): return ret_val
                else: return (None, None, None)
        
        if sol == 'CASE_1':
            ret_val = (DecisionVariables(MODEL_1, pn=case_1_pn, wn=case_1_wn, rho=case_1_rho, qn=case_1_qn), case_1_prof_man, case_1_prof_ret)
        else:
            ret_val = (DecisionVariables(MODEL_1, pn=case_2_pn, wn=case_2_wn, rho=case_2_rho, qn=case_2_qn), case_2_prof_man, case_2_prof_ret)
            
        if not is_prof_pos(ret_val[2]) or not is_prof_pos(ret_val[1]):
            if ret_val[2] >= 0 or ret_val[1] >= 0:
                pass
            return (None, None, None)
            
        return ret_val
        
if __name__ == '__main__':
    unittest.main()