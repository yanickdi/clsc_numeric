import scipy.optimize
import math

def manufacturer_derivation_case_1(wn, args):
    tau, a, s, cn = args['tau'], args['a'], args['s'], args['cn']
    #assert wn >= 0
    return (-1/2) - (cn/2) + wn + a * (tau/a)**(1/2)
    
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
    args = build_args(tau=0.1, a=0.005, s=0.0005, cn=0.01)
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