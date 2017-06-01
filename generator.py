import sys
import argparse

try:
    import scipy
except:
    sys.exit('Failed: You have to install the python scipy package.')
    
class Generator:
    """
    This class generates a set of model input data and generates an
    output file where the corresponding numerical model results are stored
    """
    def __init__(self, model_nr, output_file):
        self.model_nr = model_nr
        self.output_file = output_file
        assert model_nr == 1
    
    
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