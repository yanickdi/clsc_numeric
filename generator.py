import sys
import argparse

try:
    import scipy
except:
    sys.exit('Failed: You have to install the python scipy package.')
    
import solver
    
class Generator:
    """
    This class generates a set of model input data and generates an
    output file where the corresponding numerical model results are stored
    """
    def __init__(self, model_nr, output_file):
        self.model_nr = model_nr
        self.output_file = output_file
        if output_file == 'stdout':
            self.file_writer = StdoutFile()
        else:
            suffix = output_file.split('.')[-1]
            raise RuntimeError('other output files than stdout not implemented yet.')
        assert model_nr == 1 #TODO: Implement model 2
        
    def generate(self):
        """ Do se generation """
        model_solver = solver.ModelOneNumericalSolver() #TODO: Switch here between different models
        self.file_writer.open()
        for const_args in self.__model_one_const_generator():
            dec_vars = model_solver.optimize(const_args)
            manu_profit, ret_profit = model_solver.calc_profits(const_args, dec_vars) #TODO: HERE
            print(manu_profit)
            self.file_writer.writeSolution(const_args, dec_vars, manu_profit, ret_profit) #TODO: add manufacturer and retailer profit here
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
        line = '{}: {}, {}'.format(self.sol_nr, const_args, profit_man)
        print(line)

        
class CsvOutputFile(AbstractOutputFile):
    """ This class writes an comma separated, text based file """
    
    def __init__(self, output_file):
        self.output_file = output_file
        self.line_nr = 0
    
    def open(self):
        self.file = sys.stdout
        
    def writeSolution(self, const_args, dec_vars, profit_man, profit_ret):
        self.line_nr += 1
        line = '{}: {}, {}'.format(self.line_nr, const_args, profit_man)
        
    def close(self):
        pass
    
    
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