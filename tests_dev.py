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
        dec = solver._optimize_case_one_a(par)
        profit_man,_ = solver.calc_profits(par, dec)
        
        self.assertAlmostEqual(dec.wn, .5769703256659778)
        self.assertAlmostEqual(dec.pr, .29424258141649456)
        self.assertAlmostEqual(dec.qr, .0)
        self.assertAlmostEqual(profit_man, .029687932228693096)
        #TODO: test retailer profit
        
    def test_case_one_b(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to case one (b), but i will check the output anyway
        par = Parameter(MODEL_2, tau=.1, a=.05, s=.1, cr=.2, cn=.3, delta=.8)
        # TODO, check the case
        dec = solver._optimize_case_one_b(par)
        profit_man,_ = solver.calc_profits(par, dec)
        
        self.assertAlmostEqual(dec.wn, .5551316701949486)
        self.assertAlmostEqual(dec.pr, .48675444679663243)
        self.assertAlmostEqual(dec.qr, .06250000000000011)
        self.assertAlmostEqual(profit_man, .023662364345369585)
        #TODO: test retailer profit
        
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
        self.assertAlmostEqual(dec.roh, 1)
        self.assertAlmostEqual(dec.qn, .282407407)
        self.assertAlmostEqual(dec.qr, .0)
        self.assertAlmostEqual(profit_man,  .0861342592593)
        self.assertAlmostEqual(profit_ret, -.0356442901235)
        
    def test_case_two_b(self):
        solver = ModelTwoNumericalSolver()
        # i dont know whether this parms lead to two a, but i will check the output anyway
        par = Parameter(MODEL_2, tau=.1, a=.05, s=.1, cr=.2, cn=.3, delta=.8)
        # TODO, check the case
        dec = solver._optimize_case_two_b(par)
        profit_man, profit_ret = solver.calc_profits(par, dec)
        
        self.assertAlmostEqual(dec.wn, .667039106)
        self.assertAlmostEqual(dec.pr, .542458101)
        self.assertAlmostEqual(dec.pn, .704748603)
        self.assertAlmostEqual(dec.roh, 1)
        self.assertAlmostEqual(dec.qn, .188547486)
        self.assertAlmostEqual(dec.qr, 0.133379888)
        self.assertAlmostEqual(profit_man,  .0908519553073)
        self.assertAlmostEqual(profit_ret, -.0436009721919)
        
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
        self.assertAlmostEqual(dec.roh, 1)
        self.assertAlmostEqual(dec.qn, .26612903226)
        self.assertAlmostEqual(dec.qr, .02661290323)
        self.assertAlmostEqual(profit_man,  .0878225806452)
        self.assertAlmostEqual(profit_ret, -.0372515608741)

if __name__ == '__main__':
    unittest.main()