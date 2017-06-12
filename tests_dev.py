"""
    This file stores some tests for classes that are in development right now
"""

import unittest
from math import sqrt

from solver import ModelTwoNumericalSolver, Parameter, is_prof_pos
from generator import Generator, MemoryOutputFile, MODEL_1, MODEL_2

class TestModelTwoNumericalSolver(unittest.TestCase):
    def test_case_one_a(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to case one (a), but i will check the output anyway
        par = Parameter(MODEL_2, tau=.5, a=.01, s=.1, cr=.2, cn=.3, delta=.4)
        # TODO, check the case
        #self.assertTrue(self.__input_is_in_case_1(par))
        dec_vars = solver._optimize_case_one_a(par)
        profit_man,_ = solver.calc_profits(par, dec_vars)
        
        self.assertAlmostEqual(dec_vars.wn, .5769703256659778)
        self.assertAlmostEqual(dec_vars.pr, .29424258141649456)
        self.assertAlmostEqual(dec_vars.qr, .0)
        self.assertAlmostEqual(profit_man, .029687932228693096)
        #TODO: test retailer profit
        
    def test_case_one_b(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to case one (b), but i will check the output anyway
        par = Parameter(MODEL_2, tau=.1, a=.05, s=.1, cr=.2, cn=.3, delta=.8)
        # TODO, check the case
        dec_vars = solver._optimize_case_one_b(par)
        profit_man,_ = solver.calc_profits(par, dec_vars)
        
        self.assertAlmostEqual(dec_vars.wn, .5551316701949486)
        self.assertAlmostEqual(dec_vars.pr, .48675444679663243)
        self.assertAlmostEqual(dec_vars.qr, .06250000000000011)
        self.assertAlmostEqual(profit_man, .023662364345369585)
        #TODO: test retailer profit
    
    def test_case_one_c(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to case one (b), but i will check the output anyway
        par = Parameter(MODEL_2, tau=.1, a=.05, s=.1, cr=.2, cn=.3, delta=.8)
        # TODO, check the case
        dec_vars = solver._optimize_case_one_c(par)
        profit_man, profit_ret = solver.calc_profits(par, dec_vars)
        
        self.assertAlmostEqual(dec_vars.wn, .5551316701949)
        self.assertAlmostEqual(dec_vars.pr, .6976524961498)
        self.assertAlmostEqual(dec_vars.qr, .1581138830084)
        self.assertAlmostEqual(profit_man,  .0586018385643)
        self.assertAlmostEqual(profit_ret, -.0027016118539)

if __name__ == '__main__':
    unittest.main()