FILENAME = 'calc.txt'

import sys

from solver import Parameter, MODEL_2, Solution, ALL_CASES, ModelTwoFracQuad

def __retrFloat(str):
    return float(str.split(' -> ')[1])

def getParameter(line):
    assert line[0] == '{' and line[-1] == '}'
    strings = line[1:-1].split(', ')
    tau, a, s, cr, cn, delta = (__retrFloat(str) for str in strings)
    return Parameter(MODEL_2, tau, a, s, cr, cn, delta)
    
_CASE_STRINGS = list('Case '+str for str in ('1a', '1b', '1c', '2a', '2b', '2c'))

def __split_solutions_in_line(line):
    solutions = []
    solution_strings = line.split('}, {')
    for case_solution in solution_strings:
        case_solution = case_solution.strip('{}')
        if case_solution != '':
            case_solution_dict = {}
            convertion_error = False
            for pair in case_solution.split(', '):
                splitted = pair.split(' -> ')
                key = splitted[0]
                try:
                    val = float(splitted[1].replace('*^-', 'e-'))
                except:
                    convertion_error = True
                    break
                case_solution_dict[key] = val
            if not convertion_error:
                solutions.append(case_solution_dict)
    return solutions
    
def processSection(par, bufferedLines):
    """ Returns a Solution object """
    case_nr = 0
    solutions = [None for i in range(6)]
    in_case = False
    # scan lines and fill solutions list - this list will be a list of Strings (or None)
    for line in bufferedLines:
        if case_nr <= len(_CASE_STRINGS)-1 and line == _CASE_STRINGS[case_nr]:
            in_case = True
            case_nr += 1
        elif line == '{}' or line.startswith('{{'):
            # now this is a solution line
            solutions[case_nr-1] = line[1:-1]
        else:
            # skip this line
            pass
    
    all_solutions = []
    # now walk trough solutions and convert each element to a Solution object
    for i, line in enumerate(solutions):
        if line is None:
            line = ''
        case_solutions = __split_solutions_in_line(line)
        case = ALL_CASES[i]
        # convert them to solution objects
        for sol in case_solutions:
            wn, pr = sol['wn'], sol['pr']
            try:
                dec = ModelTwoFracQuad.create_dec_vars(wn, pr, case, par)
                prof_m, prof_r = ModelTwoFracQuad.calc_profits(par, dec)
                sol = Solution(dec, prof_m, prof_r, case)
                # validate Solution
                if ModelTwoFracQuad.is_valid(par, sol):
                    all_solutions.append(sol)
                
            except ZeroDivisionError:
                # skip this solution
                pass
    if len(all_solutions) > 0:
        # find best Solutions
        return max(all_solutions, key=lambda sol: sol.profit_man)
    else:
        return None

with open(FILENAME, 'r') as f:
    line = f.readline()
    actParameter = None
    inBetweenLines = []
    nextIsFirstLine = False
    i = 0
    while line != "":
        line = line.strip()
        if line == '--':
            nextIsFirstLine = True
        elif nextIsFirstLine:
            # we are in the first line
            actParameter = getParameter(line)
            inBetweenLines = []
            nextIsFirstLine = False
        else:
            # we are at Solution Section - buffer all lines
            if actParameter is not None:
                inBetweenLines.append(line)
        line = f.readline()
        
        # can we process anything?
        if (nextIsFirstLine or line == "") and actParameter is not None:
            sol = processSection(actParameter, inBetweenLines)
            if sol is not None:
                i += 1
            
    print(i)