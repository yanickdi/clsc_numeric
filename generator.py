import sys
import argparse

try:
    import scipy
except:
    sys.exit('Failed: You have to install the python scipy package.')

from solver import MODEL_1, MODEL_2
import solver
    
class Generator:
    """
    This class generates a set of model input data and generates an
    output file where the corresponding numerical model results are stored
    """
    def __init__(self, model_nr, output_file):
        self.model_nr = model_nr
        self.output_file = output_file
        if issubclass(type(output_file), MemoryOutputFile):
            self.file_writer = output_file
        elif output_file == 'stdout':
            self.file_writer = StdoutFile()
        else:
            suffix = output_file.split('.')[-1]
            if suffix == 'csv':
                self.file_writer = CsvOutputFile(output_file)
            else:
                raise RuntimeError('other output files than stdout not implemented yet.')
        assert model_nr == MODEL_1 #TODO: Implement model 2
        
    def generate(self):
        """ Do se generation """
        model_solver = solver.ModelOneNumericalSolver() #TODO: Switch here between different models
        self.file_writer.open()
        for const_args in self.__model_one_const_generator():
            dec_vars = model_solver.optimize(const_args)
            manu_profit, ret_profit = model_solver.calc_profits(const_args, dec_vars)
            self.file_writer.writeSolution(const_args, dec_vars, manu_profit, ret_profit)
        self.file_writer.close()
    
    def __model_one_const_generator(self):
        """ Helper generator for yielding all combinations of input data """
        for tau in drange(0, 1, 0.1):
            for a in drange(0.01, 0.1, 0.01):
                for s in drange(0, 1, .1):
                    for cn in drange(s, 1, .1):
                        yield {'tau' : round(tau, 1), 'a' : round(a, 2), 's' : round(s, 1), 'cn' : round(cn, 1)}
    
def drange(start, end, step_size):
    """ A floating point range from [start, end] with step size step_size"""
    r = start
    while r <= end:
        yield r
        r += step_size
        
        
class AbstractOutputFile:
    """ This is the abstract base class of all output classes """
    
    def open(self):
        raise NotImplementedError()
        
    def close(self):
        raise NotImplementedError()
    
    def writeSolution(self, const_args, dec_vars, profit_man, profit_ret):
        raise NotImplementedError()
        
class MemoryOutputFile:
    """ This class is used to test and store the outcome of the Generator Object """
    
    def open(self):
        self._list = []
    
    def writeSolution(self, const_args, dec_vars, profit_man, profit_ret):
        self._list.append({'const_args' : const_args, 'dec_vars' : dec_vars, 'profit_man' : profit_man, 'profit_ret' : profit_ret})
    
    def close(self):
        pass
        
    def getSolutions(self):
        """
        This method will return all stored solutions so far
        Do not call this method *before* the open() call
        
        The list contains dictionaries with the following keywords: const_args, dec_vars, profit_man, profit_ret
        """
        return self._list
        
        
class StdoutFile(AbstractOutputFile):
    """ This class writes the output to stdout """
    def __init__(self):
        self.sol_nr = 0
        
    def open(self):
        pass
    def close(self):
        pass
    
    def writeSolution(self, const_args, dec_vars, profit_man, profit_ret):
        self.sol_nr += 1
        line = '{}: {}, {}, profit of manufacturer: {}, profit or retailer: {}'.format(self.sol_nr, const_args, dec_vars, profit_man, profit_ret)
        print(line)

        
class CsvOutputFile(AbstractOutputFile):
    """ This class writes an comma separated, text based file """
    
    def __init__(self, file_name):
        self.file_name = file_name
    
    def open(self):
        self.file = open(self.file_name, 'w') # overrides file if exists
        self.file.write('tau;a;s;cn;pn;wn;roh;qn;manufacturer profit;retailer profit\n')
        
    def writeSolution(self, const_args, dec_vars, profit_man, profit_ret):
        const_str = '{};{};{};{}'.format(const_args['tau'], const_args['a'], const_args['s'], const_args['cn']) #TODO: This code doesnt support model 2
        if dec_vars == None:
            dec_str = '{};{};{};{}'.format(None, None, None, None)
            profit_man = -1.0
            profit_ret = -1.0
        else:
            dec_str =   '{};{};{};{}'.format(dec_vars['pn'], dec_vars['wn'], dec_vars['roh'], dec_vars['qn'])
        line = '{};{};{:.5f};{:.5f}'.format(const_str, dec_str, profit_man, profit_ret)
        self.file.write(line.replace('.',',') + '\n')
        
    def close(self):
        self.file.close()
    
    
def __parser_model_one_or_two(string):
    if string.lower() == 'one' or string == '1':
        return 1
    elif string.lower() == 'two' or string == '2':
        raise argparse.ArgumentTypeError('Model 2 not implemented yet.')
        #return 2
    raise argparse.ArgumentTypeError('model has to be either \'one\' or \'two\' (or 1/2)')
    
def __parser_file(string):
    if string == 'stdout':
        return string
    if len(string.split('.')) <= 1:
        raise argparse.ArgumentTypeError('The output file has to have an suffix.')
    suffix = string.split('.')[-1]
    if suffix not in ('html', 'tex', 'csv'):
        raise argparse.ArgumentTypeError('Supported output types are: .html/.tex/.csv')
    return string

if __name__ == '__main__':
    # parse the command line
    parser = argparse.ArgumentParser(description='Numeric solver for Andrea\'s Model')
    parser.add_argument('--model', type=__parser_model_one_or_two, nargs=1, required=True)
    parser.add_argument('--output', type=__parser_file, nargs=1, required=True)
    args = parser.parse_args()
    
    model_nr = args.model[0]
    output_file = args.output[0]
    generator = Generator(model_nr, output_file)
    generator.generate()