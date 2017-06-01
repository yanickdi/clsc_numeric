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
        assert model_nr == 1 #TODO: Implement model 2
        
    def generate(self):
        """ Do se generation """
        for const_args in self.__model_one_cons_generator():
            print(const_args)
    
    def __model_one_cons_generator(self):
        """ Helper generator for yielding all combinations of input data """
        for tau in drange(0, 1, 0.1):
            for a in drange(0, 0.1, 0.01):
                for s in drange(0, 1, .1):
                    for cn in drange(s, 1, .1):
                        yield {'tau' : tau, 'a' : a, 's' : s, 'cn' : cn}
                        
                        
def build_args(tau, a, s, cn):
    args =  {
        'tau' : tau,
        'a'   : a,
        's'   : s,
        'cn'  : cn
    }
    check_args(args)
    return args
    
def drange(start, end, step_size):
    """ A floating point range from [start, end] with step size step_size"""
    r = start
    while r <= end:
        yield r
        r += step_size

        
class CsvOutputFile:
    """ This class writes an comma separated, text based file """
    
    def __init__(self, output_file):
        self.output_file = output_file
        self.line_nr = 0
    
    def open(self):
        self.file = sys.stdout
        
    def writeSolution(self, const_args, dec_vars, profit_man, profit_ret):
        self.line_nr += 1
        line = '{}: {}, {}'.format(self.line_nr, const_args, profit_man)
        
    def close():
        pass
    
    
def __parser_model_one_or_two(string):
    if string.lower() == 'one' or string == '1':
        return 1
    elif string.lower() == 'two' or string == '2':
        raise argparse.ArgumentTypeError('Model 2 not implemented yet.')
        #return 2
    raise argparse.ArgumentTypeError('model has to be either \'one\' or \'two\' (or 1/2)')
    
def __parser_file(string):
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