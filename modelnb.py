from cmath import sqrt as Sqrt
import numpy as np, sys
from solver import ModelNBGridSearch, Parameter, MODEL_NB


        
if __name__ == '__main__':
    par = Parameter(MODEL_NB, tau=.9, a=.02, s=0.4*0.4, cn=.4)
    search = ModelNBGridSearch()
    sol = search.search(par)
    print(sol)