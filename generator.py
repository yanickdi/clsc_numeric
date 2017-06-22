import sys
import argparse

try:
    import scipy
except:
    sys.exit('Failed: You have to install the python scipy package.')
    
try:
    from jinja2 import Template
except:
    sys.exit('Failed: You have to install jinja2 template engine. Please see documentation')

from solver import MODEL_1, MODEL_2, Parameter, DecisionVariables, drange
import solver
    
class Generator:
    """
    This class generates a set of model input data and generates an
    output file where the corresponding numerical model results are stored
    """
    def __init__(self, model, output_file, options={}):
        self.model = model
        assert model in (MODEL_1, MODEL_2)
        self.output_file = output_file
        if issubclass(type(output_file), MemoryOutputFile):
            self.file_writer = output_file
        elif output_file == 'stdout':
            self.file_writer = StdoutFile()
        else:
            suffix = output_file.split('.')[-1]
            if suffix == 'csv':
                german_comma = options.get('german_comma', False)
                self.file_writer = TemplateOutputFile(output_file, 'templates/model.tpl.csv', model, options=options)
            elif suffix == 'html':
                self.file_writer = TemplateOutputFile(output_file, 'templates/model.tpl.html', model, options=options)
            else:
                raise RuntimeError('other output files than stdout not implemented yet.')
                
    def generate(self):
        if self.model == MODEL_1:
            self._generate_model_one()
        elif self.model == MODEL_2:
            self._generate_model_two()
        else:
            raise RuntimeError('model not known')
            
    def _generate_model_two(self):
        """ Generator of model two """
        model_solver = solver.ModelTwoNumericalSolver()
        self.file_writer.open()
        for par in self.__model_two_par_generator():
            sol = model_solver.optimize(par)
            if sol == None:
                self.file_writer.writeSolution(par, None)
            else:
                self.file_writer.writeSolution(par, sol)
        self.file_writer.close()
        
    def _generate_model_one(self):
        """ Do se generation """
        model_solver = solver.ModelOneNumericalSolver()
        self.file_writer.open()
        for par in self.__model_one_par_generator():
            sol = model_solver.optimize(par)
            if sol == None:
                self.file_writer.writeSolution(par, None)
            else:
                self.file_writer.writeSolution(par, sol)
        self.file_writer.close()
    
    def __model_two_par_generator(self):
        """ Helper generator for yielding all combinations of input data for model 2"""
        for tau in drange(0, 1, 0.1):
            for a in drange(0.01, 0.1, 0.01):
                for s in drange(0, 1, .1):
                    for cr in drange(s, 1, .1):
                        for cn in drange(cr, 1, .1):
                            for delta in drange(cr, 1, .1):
                                if delta == 0: continue
                                yield Parameter(
                                    MODEL_2, tau=round(tau, 1), a=round(a, 2), s=round(s, 1),
                                    cr=round(cr, 1), cn=round(cn, 1), delta=round(delta, 1))
                        
    def __model_one_par_generator(self):
        """ Helper generator for yielding all combinations of input data for model 2"""
        for tau in drange(0, 1, 0.1):
            for a in drange(0.01, 0.1, 0.01):
                for s in drange(0, 1, .1):
                    for cn in drange(s, 1, .1):
                        yield Parameter(MODEL_1, tau=round(tau, 1), a=round(a, 2), s=round(s, 1), cn=round(cn, 1))
        
        
class AbstractOutputFile:
    """ This is the abstract base class of all output classes """
    
    def open(self):
        raise NotImplementedError()
        
    def close(self):
        raise NotImplementedError()
    
    def writeSolution(self, par, sol):
        raise NotImplementedError()
        
class TemplateOutputFile(AbstractOutputFile):
    """ This class uses jinja2 template engine to store all output as an text based file """

    def __init__(self, outfile, tplfile, model, options={}):
        self.outfile = outfile
        self.tplfile = tplfile
        self.model = model
        self.options = options
        self.solutions = []
    
    def open(self):
        with open(self.tplfile, 'r') as f:
            self.template = Template(f.read())
            
    def writeSolution(self, par, sol):
        if not (self.options.get('only_valid_solutions', False) and sol == None):
            self.solutions.append( {'par' : par, 'sol' : sol})
        
    def close(self):
        with open(self.outfile, 'w') as f:
            model = 'MODEL_1' if self.model == MODEL_1 else 'MODEL_2'
            f.write(self.template.render({
                'model' : model,
                'solutions' : self.solutions,
                'options' : self.options
                }))
        
class MemoryOutputFile:
    """ This class is used to test and store the outcome of the Generator Object """
    
    def __init__(self, callback=None):
        self.callback = callback
    
    def open(self):
        self._list = []
    
    def writeSolution(self, par, sol):
        if self.callback is not None:
            self.callback(par, sol)
        else:
            self._list.append({'par' : par, 'sol' : sol})
    
    def close(self):
        pass
        
    def getSolutions(self):
        """
        This method will return all stored solutions so far
        Do not call this method *before* the open() call
        
        The list contains dictionaries with the following keywords: const_args, dec_vars, profit_man, profit_ret
        """
        assert self.callback is None
        return self._list
        
        
class StdoutFile(AbstractOutputFile):
    """ This class writes the output to stdout """
    def __init__(self):
        self.sol_nr = 0
        
    def open(self):
        pass
    def close(self):
        pass
    
    def writeSolution(self, par, sol):
        self.sol_nr += 1
        line = '{}: {}, {}'.format(self.sol_nr, par, sol)
        print(line)
    
def __parser_model_one_or_two(string):
    if string.lower() == 'one' or string == '1':
        return MODEL_1
    elif string.lower() == 'two' or string == '2':
        return MODEL_2
    raise argparse.ArgumentTypeError('model has to be either \'one\' or \'two\' (or 1/2)')
    
def __parser_file(string):
    if string == 'stdout':
        return string
    if len(string.split('.')) <= 1:
        raise argparse.ArgumentTypeError('The output file has to have a suffix.')
    suffix = string.split('.')[-1]
    if suffix not in ('html', 'tex', 'csv'):
        raise argparse.ArgumentTypeError('Supported output types are: .html/.tex/.csv')
    return string

if __name__ == '__main__':
    # parse the command line
    parser = argparse.ArgumentParser(description='Numeric solver for Andrea\'s Model')
    parser.add_argument('-model', type=__parser_model_one_or_two, nargs=1, required=True)
    parser.add_argument('-output', type=__parser_file, nargs=1, required=True)
    parser.add_argument('--german-comma', action='store_true')
    parser.add_argument('--only-valid-solutions', action='store_true')
    args = parser.parse_args()
    
    model = args.model[0]
    output_file = args.output[0]
    generator = Generator(model, output_file, options={
        'german_comma' : args.german_comma,
        'only_valid_solutions' : args.only_valid_solutions})
    generator.generate()