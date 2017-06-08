"""
    This file stores some tests for classes that are in development right now
"""

import unittest
from math import sqrt

from solver import build_args, check_args, ModelTwoNumericalSolver, is_prof_pos
from generator import Generator, MemoryOutputFile, MODEL_1, MODEL_2

class TestModelTwoNumericalSolver(unittest.TestCase):
    def test_case_one_a(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to case one (a), but i will check the output anyway
        const_args = build_args(tau=0.5, a=0.01, s=0.1, cr=0.2, cn=0.3, delta=0.4)
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        # self checking the my test input variables..
        self.assertTrue(self.__input_is_in_case_1(const_args))
        
        analyitcal_profit_manufacturer = ((1-cn)**2 / 8) - (1/2)*(1+cn-2*s) * (tau*a)**(1/2) + (1/2)*a*tau
        profit_solver_manufacturer,_ = solver.calc_profits(const_args, solver.optimize(const_args))
        self.assertAlmostEqual(analyitcal_profit_manufacturer, profit_solver_manufacturer)
        
if __name__ == '__main__':
    unittest.main()