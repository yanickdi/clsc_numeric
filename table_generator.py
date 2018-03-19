import sys, os

from jinja2 import FileSystemLoader, Environment
from solver import MODEL_1, MODEL_2, MODEL_NB, MODEL_1_QUAD, MODEL_2_QUAD, Parameter, SolverProxy
import utils

ALLOWED_MODELS = ('modelo', 'modelnb', 'modeloq')

class Factorial:
    tau_lh   = lambda : (.15, .5)
    s_lh     = lambda cn: (0, cn/2)
    cr_lh    = lambda cn: (.1*cn, .4*cn)
    delta_lh = lambda : (.5, .85)
    cn_lh    = lambda : (.1, .5)
    a_lh     = lambda : (.001, .01) 

    LOW, HIGH = 'LOW', 'HIGH'
    
    def __init__(self, model):
        self._model = model
        self._solverProxy = SolverProxy()
        self._calcTable()
        
    def _calcTable(self):
        self.__ff_table = [[None for j in range(8)] for i in range(8)]
        lohi = (0, 1)
        for i, (s, cr, a) in enumerate([(s, cr, a) for s in lohi for cr in lohi for a in lohi]):
            for j, (cn, delta, tau) in enumerate([(cn, delta, tau) for cn in lohi for delta in lohi for tau in lohi]):
                tau_val   = Factorial.tau_lh()[tau]
                delta_val = Factorial.delta_lh()[delta]
                cn_val    = Factorial.cn_lh()[cn]
                cr_val    = Factorial.cr_lh(cn_val)[cr]
                s_val     = Factorial.s_lh(cn_val)[s]
                a_val     = Factorial.a_lh()[a]
                
                par_o = Parameter(MODEL_2, tau=tau_val, a=a_val, s=s_val, cr=cr_val, cn=cn_val, delta=delta_val)
                par_n = Parameter(MODEL_1, tau=tau_val, a=a_val, s=s_val, cn=cn_val)
                par_nb = Parameter(MODEL_NB, tau=tau_val, a=a_val, s=s_val, cn=cn_val)
                par_oq = Parameter(MODEL_2_QUAD, tau=tau_val, a=a_val, s=s_val, cr=cr_val, cn=cn_val, delta=delta_val)
                par_nq = Parameter(MODEL_1_QUAD, tau=tau_val, a=a_val, s=s_val, cn=cn_val)
                
                
                if self._model == 'modelo':
                    par_1, par_2 = par_o, par_n
                elif self._model == 'modelnb':
                    par_1, par_2 = par_o, par_nb
                elif self._model == 'modeloq':
                    par_1, par_2 = par_oq, par_nq
                
                # solve
                self.__ff_table[i][j] = {
                    'par_1' : par_1, 'par_2': par_2,
                    'sol_1' : self._solverProxy.read_or_calc_and_write(par_1),
                    'sol_2' : self._solverProxy.read_or_calc_and_write(par_2)
                }
                if self.__ff_table[i][j]['sol_2'] is None or self.__ff_table[i][j]['sol_2'].profit_ret < 0.00005:
                    self.__ff_table[i][j]['sol_2'] = None
                print(par_2)
                self._solverProxy.commit()
                
    def _latex_percentage(self, value):
        if self._isLatex:
            return '{:.2f}\\%'.format(value * 100)
        else:
            return '{:.4f}'.format(value)
                
    def template_table_val(self, table, i, j):
        # sol_1 is online store!
        par_1, par_2 = self.__ff_table[i][j]['par_1'], self.__ff_table[i][j]['par_2']
        sol_1, sol_2 = self.__ff_table[i][j]['sol_1'], self.__ff_table[i][j]['sol_2']
        percentage = self._latex_percentage
        
        if table == 'case':
            case_1 = '-' if sol_1 is None else sol_1.case
            case_2 = '-' if sol_2 is None else sol_2.case
            return '{}/{}'.format(case_1, case_2)
            
        if None in (sol_1, sol_2): return '-'
        if table == 'profits':
            return percentage(sol_1.profit_man / sol_2.profit_man)
        elif table == 'retailerprof':
            return percentage(sol_1.profit_ret / sol_2.profit_ret)
        elif table == 'rho':
            return percentage(sol_1.dec.rho / sol_2.dec.rho)
        elif table == 'price_new':
            return percentage(sol_1.dec.pn / sol_2.dec.pn)
        elif table == 'wholesale_price':
            return percentage(sol_1.dec.wn / sol_2.dec.wn)
        elif table == 'restocking_price':
            return percentage((sol_2.dec.b - par_2.s) / (sol_2.dec.wn - par_2.s))
        elif table == 'retailerprof2':
            return sol_2.profit_ret
     
    def getTemplateVariables(self):
        return {
            'esc' : utils.escape_tex,
            'ft' : self.template_table_val}
            
    def render_template(self, tpl_filename, out_filename, isLatex):
        self._isLatex = isLatex
        env = Environment(loader = FileSystemLoader('templates'))
        template = env.get_template(tpl_filename)
        with open(out_filename, 'w', newline='\n') as f:
            renderedString = template.render(self.getTemplateVariables())
            f.write(renderedString)
            
def get_file_list(model):
    if model == 'modelo':
        return {
            'tex_template' : 'template_table_model_o_vs_n.tex',
            'tex_output'   : 'output/table_model_o_vs_n.tex',
            'csv_template' : 'template_table_model_o_vs_n.csv',
            'csv_output'   : 'output/table_model_o_vs_n.csv'
            }
    elif model == 'modelnb':
        return {
            'tex_template' : 'template_table_model_o_vs_nb.tex',
            'tex_output'   : 'output/table_model_o_vs_nb.tex',
            'csv_template' : 'template_table_model_o_vs_nb.csv',
            'csv_output'   : 'output/table_model_o_vs_nb.csv'
            }
    elif model == 'modeloq':
        return {
            'tex_template' : 'template_table_model_oq_vs_nq.tex',
            'tex_output'   : 'output/table_model_oq_vs_nq.tex',
            'csv_template' : 'template_table_model_oq_vs_nq.csv',
            'csv_output'   : 'output/table_model_oq_vs_nq.csv'
            }

def main():
    if len(sys.argv) != 2 or sys.argv[1] not in ALLOWED_MODELS:
        print('usage {} modelo/modelnb/modeloq'.format(os.path.basename(sys.argv[0])))
    model = sys.argv[1]
    fullFactorial = Factorial(model)
    file_list = get_file_list(model)
    fullFactorial.render_template(file_list['tex_template'], file_list['tex_output'],
        isLatex=True)
    fullFactorial.render_template(file_list['csv_template'], file_list['csv_output'],
        isLatex=False)

if __name__ == '__main__':
    main()