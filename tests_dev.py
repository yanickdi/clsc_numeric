"""
    This file stores some tests for classes that are in development right now
"""

import unittest
from math import sqrt

from solver import ModelTwoNumericalSolver, Parameter, is_prof_pos
from generator import Generator, MemoryOutputFile, MODEL_1, MODEL_2

class TestGeneratorModelTwo(unittest.TestCase):
    def test_model_two(self):
        #mof = MemoryOutputFile()
        generator = Generator(MODEL_2, 'stdout')
        generator.generate()
        
    def test_optimize_instance_b(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to two a, but i will check the output anyway
        par = Parameter(MODEL_2, tau=.1, a=.01, s=.1, cr=.2, cn=.3, delta=.8)
        sol = solver.optimize(par)
        self.assertIsNotNone(sol)
        self.assertAlmostEqual(sol.profit_man, 0.0878225806452)
   
if __name__ == '__main__':
    unittest.main()