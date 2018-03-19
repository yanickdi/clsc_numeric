import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from solver import ModelOneNumericalSolver

from solver import Parameter, MODEL_1, _CASE_ONE, _CASE_TWO

def plot_model_n_rhos():
    # calculate 50 instances:
    NUM_POINTS = 50
    a_values   = np.linspace(0.0075, 0.0085, NUM_POINTS)
    rho_values = np.zeros(NUM_POINTS)
    for i, a in enumerate(a_values):
        par = Parameter(MODEL_1, tau=0.15, cn=0.1, s = 0.5*0.1, a=a)
        solution = ModelOneNumericalSolver.solve(par)
        rho_values[i] = solution.dec.rho
    # plot a curve, where rho is a function of a:
    plt.plot(a_values, rho_values, marker='o', markersize=4, label=r'$\rho$')
    plt.legend(); plt.xlabel('a (costs of effort)'); plt.ylabel(r'$\rho$ (effort)'); plt.show()
    
    
def plot_model_n_manufacturer(showLine=True, points=[], pointText=True):
    par = Parameter(MODEL_1, tau=0.15, cn=0.1, s = 0.5*0.1, a=0.005)
    # print whole profit function of the manufacturer
    NUM_POINTS = 100
    all_wn = np.linspace(0, 1, NUM_POINTS)
    all_profit = np.zeros(NUM_POINTS)
    for i, wn in enumerate(all_wn):
        profit = ModelOneNumericalSolver.get_manufacturer_profit(par, wn)
        all_profit[i] = profit
    if showLine:
        plt.plot(all_wn, all_profit)
        
    # plot points
    if len(points) > 0:
        for wn in points:
            profit = ModelOneNumericalSolver.get_manufacturer_profit(par, wn)
            plt.plot(wn, profit, marker='o', color='red')
            pn, rho = ModelOneNumericalSolver.get_retailer_decision(par, wn)
            latextext = r'${p_n}^*=$'+'{}'.format(pn) + r', ${\rho}^*=$' + '{:.2f}'.format(rho)
            if pointText: plt.text(wn, profit-0.01, latextext, size=8)
    # knee?
    #plt.annotate('knee?', xy=(.63, .081), xytext=(.43, .05),
    #        arrowprops=dict(facecolor='black', shrink=0.05),
    #        )
    
    #plt.text(0.25, 0.05, '$w_n$=0.2')
    #plt.text(0.25, 0.04, r'${p_n}^*(w_n=0.2)= ?$')
    #plt.text(0.25, 0.03, r'${\rho}^*(w_n=0.2)= ?$')
    #plt.vlines(0.2, -1, 1, linestyle='--', color='grey', alpha=.5)
    #plt.yticks([])
    #plt.ylabel(r'profit($w_n, {p_n}^*, {\rho}^*$)')
    
    plt.title(r"Manufacturer's Profit $\tau=0.15, c_n =0.1, s=\frac{c_n}{2}, a=0.005$")
    plt.xlim(0, 1)
    plt.ylim(np.nanmin(all_profit)*1.1, np.nanmax(all_profit)*1.1)
    plt.xlabel('$w_n$')
    plt.show()
    
    
def plot_model_n_manufacturer_both_cases():
    par = Parameter(MODEL_1, tau=0.15, cn=0.1, s = 0.5*0.1, a=0.005)
    # print whole profit function of the manufacturer
    NUM_POINTS = 100
    all_wn = np.linspace(0, 1, NUM_POINTS)
    all_profit_one = np.zeros(NUM_POINTS)
    all_profit_two = np.zeros(NUM_POINTS)
    all_profit = np.zeros(NUM_POINTS)
    knee_point = -1
    for i, wn in enumerate(all_wn):
        all_profit_one[i] = ModelOneNumericalSolver.get_manufacturer_profit(par, wn, case=_CASE_ONE)
        all_profit_two[i] = ModelOneNumericalSolver.get_manufacturer_profit(par, wn, case=_CASE_TWO)
        all_profit[i] = ModelOneNumericalSolver.get_manufacturer_profit(par, wn)
        if all_profit_two[i] == all_profit[i] and knee_point < 0:
            knee_point = i
    plt.plot(all_wn, all_profit)
    plt.text(0.35, 0.085, r'$\rho > 1$', size=8, color='#1f77b4', rotation=38)
    plt.text(0.85, 0.048, r'$\rho = 1$', size=8, color='#1f77b4', rotation=-54)
    plt.plot(all_wn[10:knee_point], all_profit_two[10:knee_point], linestyle='--', color='grey')
    plt.text(0.30, 0.059, r'$\rho = 1$', size=8, color='grey', rotation=38)
    
    plt.title(r"Manufacturer's Profit $\tau=0.15, c_n =0.1, s=\frac{c_n}{2}, a=0.005$")
    plt.xlim(0, 1)
    plt.ylim(np.nanmin(all_profit)*1.1, np.nanmax(all_profit)*1.1)
    plt.xlabel('$w_n$')
    plt.show()
    
def plot_model_n_retailer(wn):
    par = Parameter(MODEL_1, tau=0.15, cn=0.1, s = 0.5*0.1, a=0.005)
    NUM_PN = 10; NUM_RHO = 100
    values = np.zeros([int(NUM_PN*NUM_RHO), 3])
    i = 0
    for pn in np.linspace(0, 1, NUM_PN):
        for rho in np.linspace(1, 20, NUM_RHO):
            profit = ModelOneNumericalSolver.get_retailer_profit(par, wn, pn, rho)
            values[i, :] = [pn, rho, profit]
            i += 1
    # plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_trisurf(values[:,0], values[:,1], values[:,2])
    ax.set_xlabel('pn'); ax.set_ylabel('rho')
    ax.set_title('Retailer\'s Profit,  $w_n$={}'.format(wn), loc='left')    
    # plot optimum
    best_pn, best_rho = ModelOneNumericalSolver.get_retailer_decision(par, wn)
    best_profit = ModelOneNumericalSolver.get_retailer_profit(par, wn, best_pn, best_rho)
    latexlegend = r'${p_n}^*=$'+'{}'.format(best_pn) + r', ${\rho}^*=$' + '{:.2f}'.format(best_rho)
    ax.scatter(best_pn, best_rho, best_profit, values[:,2], color='red', label=latexlegend)
    plt.legend()
    plt.show()
    
def main():
    #plot_model_n_rhos()
    #plot_model_n_manufacturer(showLine=False, points=[], pointText=True)
    #plot_model_n_manufacturer(showLine=True, points=[.2, .4, .6, .8, .9], pointText=False)
    #plot_model_n_retailer(0.2)
    plot_model_n_manufacturer_both_cases()
    
if __name__ == '__main__':
    main()