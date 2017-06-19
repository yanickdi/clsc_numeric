"""
    This file stores some tests for classes that are in development right now
"""

import unittest
from math import sqrt
from tempfile import TemporaryDirectory
from os.path import join, isfile

from solver import ModelTwoNumericalSolver, Parameter, is_prof_pos
from generator import Generator, MemoryOutputFile, MODEL_1, MODEL_2, HtmlOutputFile

class TestGeneratorModelTwo(unittest.TestCase):
    def test_model_two(self):
        generator = Generator(MODEL_1, 'test_out.html')
        generator.generate()

class TestHtmlOutputFile(unittest.TestCase):
    def test_template_generation(self):
        with TemporaryDirectory() as tmpdir:
            test_filen = join(tmpdir, 'test.html')
            file = HtmlOutputFile(test_filen, MODEL_1)
            file.open()
            file.close()
            self.assertTrue(isfile(test_filen))
  
if __name__ == '__main__':
    unittest.main()