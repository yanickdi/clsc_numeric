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
   
if __name__ == '__main__':
    unittest.main()