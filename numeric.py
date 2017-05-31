import scipy.optimize
import sys

class ModelNormalNumericalSolver:
    """
        This class offers methods to solve Model 1 (Without Online Store of the Manufacturer) numerically
    """
    def __init__(self):
        pass
    
    def profit_manufacturer(self, dec_vars, const_args):
        """ Returns the numeric result of the profit of the manufacturer having set all decision variables"""
        print(dec_vars)
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        wn, pn, roh, qn = dec_vars['wn'], dec_vars['pn'], dec_vars['roh'], dec_vars['qn']
        return qn * (wn * (1- tau/roh) - cn + (tau/roh) * s)
    
    def optimize(self, const_args):
        """
        This is the core method of this class. It will return all four
        decision variables to maximize the retailer's profit (with
        respect to the profit maximization condition of the retailer)
        """
        ## test two cases:
        #       case 1 - roh is >= 1
        #       case 2 - roh is == 1
        
        #case 1:
        dec_vars_case_1 = self._optimize_case_one(const_args)
        return dec_vars_case_1
    
    def _optimize_case_one(self, const_args):
        """
        Returns a dec_vars dict having wn,pn and roh stored
        """
        tau, a, s, cn = const_args['tau'], const_args['a'], const_args['s'], const_args['cn']
        
        # defining our helper function:
        def __manufacturer_derivation_case_1(wn):
            return (-1/2) - (cn/2) + wn + a * (tau/a)**(1/2)
        
        # let scipy do the job:
        opt = scipy.optimize.fsolve(__manufacturer_derivation_case_1, x0=0.5, full_output=True)
        wn = opt[0][0]
        pn = (1 + wn) / 2
        roh = (1-wn) * (1/2) * (tau/a)**2
        return {'wn' : wn, 'pn' : pn, 'roh' : roh, 'qn' : 1 - pn}
    
def check_args(args):
    tau, a, s, cn = args['tau'], args['a'], args['s'], args['cn']
    assert 0 <= tau <= 1
    assert 0 <= 0.01
    assert s <= cn <= 1
    assert 0 <= s <= cn
    assert s <= cn
    
def build_args(tau, a, s, cn):
    args =  {
        'tau' : tau,
        'a'   : a,
        's'   : s,
        'cn'  : cn
    }
    check_args(args)
    return args
    
    
if __name__ == '__main__':
    solver = ModelNormalNumericalSolver()
    const_args = build_args(tau=0.1, a=0.005, s=0.0005, cn=0.01)
    dec_vars = solver.optimize(const_args)
    print(solver.profit_manufacturer(dec_vars, const_args))
    sys.exit()
    
    tau, a, s, cn = args['tau'], args['a'], args['s'], args['cn']
    if cn <= 1 - 4*(1-tau/2)*(a/tau)**(1/2):
        print('case 1(a): p >= 1')
        wn = (1+cn)/2 - (a*tau)**(1/2)
        print('optimales wn muesste sein:', wn)
        print('anaL ', (1/2) * (1+cn - 2*a * (tau/a)**(1/2)))
    else:
        print('case 2(a): p == 1')
    # this one could be done analyitcal
    opt = scipy.optimize.fsolve(manufacturer_derivation_case_1, x0=0.5, args=args, full_output=True)
    print(opt[0])
    #print(opt)