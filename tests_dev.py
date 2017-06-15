"""
    This file stores some tests for classes that are in development right now
"""

import unittest
from math import sqrt

from solver import ModelTwoNumericalSolver, Parameter, is_prof_pos
from generator import Generator, MemoryOutputFile, MODEL_1, MODEL_2

class TestGeneratorModelTwo(unittest.TestCase):
    def test_model_two(self):
        def __test_callback(par, dec_vars, profit_man, profit_ret):
            if dec_vars is not None:
                self.assertTrue(dec_vars.lambda1 >= 0, par)
            else:
                print('blub')
        mof = MemoryOutputFile(callback=__test_callback)
        generator = Generator(MODEL_2, mof)
        generator.generate()
        
   
if __name__ == '__main__':
    unittest.main()