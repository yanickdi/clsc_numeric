import unittest

from numeric import build_args, check_args, ModelNormalNumericalSolver

class TestModelNormalNumericalSolver(unittest.TestCase):

    def test_case_a(self):
        solver = ModelNormalNumericalSolver()
        # this args should lead to case one (roh is gte 1)
        const_args = build_args(tau=0.1, a=0.005, s=0.0005, cn=0.01)
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        # self checking the my test input variables..
        self.assertTrue(cn <= 1 - 4*(1-tau/2)*(a/tau)**(1/2))
        
        analyitcal_profit_manufacturer = ((1-cn)**2 / 8) - (1/2)*(1+cn-2*s) * (tau*a)**(1/2) + (1/2)*a*tau
        profit_solver_manufacturer = solver.profit_manufacturer(solver.optimize(const_args), const_args)
        #print(analyitcal_profit_manufacturer)
        #print(profit_solver_manufacturer)
        self.assertAlmostEqual(analyitcal_profit_manufacturer, profit_solver_manufacturer)
        
        
if __name__ == '__main__':
    unittest.main()