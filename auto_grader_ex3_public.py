## auto checker for moodle submissions


import os
import subprocess
from subprocess import PIPE
import sys
import re

required_funcs=["def transcribe(self, dna_seq)", "class MutantCell(StemCell)", "class CancerCell(MutantCell)", "prosite_to_python(self, pattern_dict)", "patterns_to_domains(self, pattern_file)", "def classify(self, seq_list, csv_file)", "class SequenceClassifier", "def __init__(self, pattern_file)", "def __init__(self, genome, num_mutations=0)", "def __init__(self, genome, num_mutations)"]

# Dictionary where keys are input strings (will be transformed to list before running) and values are solutions
# Where an assertion should be raised value is None
example_inputs = {os.path.join(os.path.dirname(sys.argv[0]), "ex3_inputs.json") : re.compile("Original cell: <MutantCell, 2, 3>;+Final number of cells: (200|198);+Protein repertoire size: [0-9]+;+Mutations: 48;", re.MULTILINE), os.path.join(os.path.dirname(sys.argv[0]), "ex3_bad_inputs.json") : None}

# run actual commands and get output
def run_code(python_command, file, curr_input):
    # create initical list
    cmd=[python_command, file]
    # add inputs - must be done separately because it is an in-place function
    cmd.extend(curr_input.split())
    # run
    p = subprocess.run(args=cmd, 
            stdout=PIPE, stderr=PIPE, encoding='utf-8',timeout=15)
    return p.returncode, p.stdout, p.stderr

# wrapper
def run_test(python_command, file, curr_input):
    # initialize
    result = None
    try:
        result = run_code(python_command, file, curr_input)
    # if there is a problem with actual running - notify immediately
    except Exception as e:
        print ("error", e)
    return result

def fix_output(out_str):
    return out_str.replace("'","").replace(" ","").replace("\t", "")

# compare gold standard to tested assignment for given input
def run_tests(python_command, file, curr_input, gold_result):
    # notify and run
    print("Test input:", curr_input)
    return_code, result_stdout, result_stderr = run_test(python_command, file, curr_input)
    
    # make sure return code is same for both files
    # if return code is not 0 and gold strandard is not AssertionError (not None)
    if return_code != 0 and gold_result:
        print("RUN ERROR")
    # else if the expected output is AssertionError
    elif not gold_result:
        # if last line, not including empty lines, is AssertionError, accept
        if [f for f in result_stderr.replace("'","").replace(" ","").replace("\t", "").split("\n") if f] and [f for f in result_stderr.replace("'","").replace(" ","").replace("\t", "").split("\n") if f][-1].startswith("AssertionError:"):
            print("PASS: Assertion")
        elif [f for f in result_stderr.replace("'","").replace(" ","").replace("\t", "").split("\n") if f] and [f for f in result_stderr.replace("'","").replace(" ","").replace("\t", "").split("\n") if f][-1].startswith("AssertionError"):
            print("FAILED: Assertion Must Print Offending Variable")
        else:
            print("RUN ERROR: Assertion Expected")
    else:
        # make sure output is the same
        if gold_result.fullmatch(";".join(result_stdout.split("\n"))):
            print("PASS: Correct Output")
        else:
            # if not - notify
            print("FAILED\n")
            print("\n".join(["Expected Output:", gold_result.pattern.replace(";+", "\n"), "\nGiven Output:", result_stdout]))
    # print extra line for clarity
    print()

if __name__ == '__main__':
        # get input file
        full_submission_file = sys.argv[1]
        
        # get current running python path
        python_command = sys.executable
        
        # iterate over tests and solutions
        for test, solution in example_inputs.items():
            # run test for each input
            run_tests(python_command, full_submission_file, test, solution)

        # check if correct functions/classes were used
        print("Checking function format...")
        # open file
        with open(full_submission_file) as f:
            # read
            f = f.read()
            # print proper name if needed
            [print("Class or function is not properly named, expected:", fun) for fun in required_funcs if not fun + ":" in f]
        
        # final message
        print("\nAll done. Don't forget to check your return values and document properly. Stay calm and keep coding!")
