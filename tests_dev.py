"""
    This file stores some tests for classes that are in development right now
"""

import unittest
from math import sqrt
from tempfile import TemporaryDirectory
from os.path import join, isfile

from solver import ModelTwoFracQuadNumericalSolver, Parameter, is_prof_pos
from generator import Generator, MemoryOutputFile, MODEL_1, MODEL_2

class TestModelTwoFracQuadNumericalSolver(unittest.TestCase):
    def test_optimize(self):
        par = Parameter(MODEL_2, tau=0.1, a=0.05, s=0.1, cr=0.2, cn=0.3, delta=0.5)
        model = ModelTwoFracQuadNumericalSolver()
        sol = model.optimize(par)
        self.assertAlmostEqual(sol.dec.wn,  0.575)
        self.assertAlmostEqual(sol.dec.pr,  0.3583333333333332)
        self.assertAlmostEqual(sol.dec.mu1, 0.033333333333333326)
  
if __name__ == '__main__':
    unittest.main()