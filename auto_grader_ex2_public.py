## auto checker for moodle submissions


import os
import subprocess
from subprocess import PIPE
import sys

# path to "gold" solution
gold_solution_path = "/home/alu/fulther/Scripts/scripts_python/HW"
gold_solution_filename1 = os.path.join(gold_solution_path, "exercise2.py")

required_funcs=["def find_srr(self, dna_seq)", "def transcribe(self, dna_seq)", "def translate(self, rna_seq)", "class Polymerase", "def __init__(self, type, error_rate = 0)", "class Ribosome", "def __init__(self, genetic_code, start_codons)", "synthesize(self, rna_seq)", "class Cell", "def __init__(self, name, genome, num_copies, genetic_code, start_codons, division_rate)", "def mitosis(self)", "def meiosis(self)", "def repertoire(self)", "class ProkaryoticCell(Cell)", "class NeuronCell(EukaryoticCell)", "class StemCell(EukaryoticCell)", "class EukaryoticCell(Cell)"]

# Dictionary where keys are input strings (will be transformed to list before running) and values are solutions
# Where an assertion should be raised value is None
example_inputs = {"StemCell 3 100 TTGATCTGCATGTTCATGAT ATCAAATCAAATCAAATCAA TTACATCNTN" : None, 
"StemCell 3 100 TTGATCTGCATGTTCATGAT ATCAAATCAAATCAAATCAA TTACATCAT" : "Original cell: <StemCell, 2, 3>\nFinal number of cells: 27\nRepertoire: [('No simple repeats in DNA sequence', 'AUCAUGAACAUGCAGAUCAA', 'MNMQI'), ('A,3;AAATC,3;AATCA,3;ATCAA,4;CAAAT,3;TCAAA,3', 'UUGAUUUGAUUUGAUUUGAU', 'Non-coding RNA'), ('No simple repeats in DNA sequence', 'AUGAUGUAA', 'MM')]\nUndergoing meiosis...\nFirst cell genome: ['TTGATCTGCATGTTCATGAT', 'ATCAAATCAAATCAAATCAA', 'TTACATCAT']\nFirst cell repertoire: [('No simple repeats in DNA sequence', 'AUCAUGAACAUGCAGAUCAA', 'MNMQI'), ('A,3;AAATC,3;AATCA,3;ATCAA,4;CAAAT,3;TCAAA,3', 'UUGAUUUGAUUUGAUUUGAU', 'Non-coding RNA'), ('No simple repeats in DNA sequence', 'AUGAUGUAA', 'MM')]\nSecond cell genome: ['ATCATGAACATGCAGATCAA', 'TTGATTTGATTTGATTTGAT', 'ATGATGTAA']\nSecond cell repertoire: [('No simple repeats in DNA sequence', 'UUGAUCUGCAUGUUCAUGAU', 'MFM'), ('ATTTG,3;GATTT,3;T,3;TGATT,3;TTGAT,4;TTTGA,3', 'AUCAAAUCAAAUCAAAUCAA', 'Non-coding RNA'), ('No simple repeats in DNA sequence', 'UUACAUCAU', 'Non-coding RNA')]", 
"Eggs 3 100 TTGATCTGCATGTTCATGAT ATCAAATCAAATCAAATCAA TTACATCAT" : None, 
"NeuronCell 3.5 100 TTGATCTGCATGTTCATGAT ATCAAATCAAATCAAATCAA TTACATCAT" : None, 
"NeuronCell 3 10.0 TTGATCTGCATGTTCATGAT ATCAAATCAAATCAAATCAA TTACATCAT" : None, 
"NeuronCell 3 10 TTGATCTGCATGTTCATGAT ATCAAATCAAATCAAATCAA TTACATCAT" : "Original cell: <NeuronCell, 2, 2>\nFinal number of cells: 8\nRepertoire: [('No simple repeats in DNA sequence', 'AUCAUGAACAUGCAGAUCAA', 'MNMQI'), ('A,3;AAATC,3;AATCA,3;ATCAA,4;CAAAT,3;TCAAA,3', 'UUGAUUUGAUUUGAUUUGAU', 'Non-coding RNA'), ('No simple repeats in DNA sequence', 'AUGAUGUAA', 'MM')]\nUndergoing meiosis...\nFirst cell genome: ['TTGATCTGCATGTTCATGAT', 'ATCAAATCAAATCAAATCAA', 'TTACATCAT']\nFirst cell repertoire: [('No simple repeats in DNA sequence', 'AUCAUGAACAUGCAGAUCAA', 'MNMQI'), ('A,3;AAATC,3;AATCA,3;ATCAA,4;CAAAT,3;TCAAA,3', 'UUGAUUUGAUUUGAUUUGAU', 'Non-coding RNA'), ('No simple repeats in DNA sequence', 'AUGAUGUAA', 'MM')]\nSecond cell genome: ['ATCATGAACATGCAGATCAA', 'TTGATTTGATTTGATTTGAT', 'ATGATGTAA']\nSecond cell repertoire: [('No simple repeats in DNA sequence', 'UUGAUCUGCAUGUUCAUGAU', 'MFM'), ('ATTTG,3;GATTT,3;T,3;TGATT,3;TTGAT,4;TTTGA,3', 'AUCAAAUCAAAUCAAAUCAA', 'Non-coding RNA'), ('No simple repeats in DNA sequence', 'UUACAUCAU', 'Non-coding RNA')]", 
"NeuronCell 4 13 TTTTCATTCATATCATGATCA TTACATCAT ATCATCGCTACTATACTATACTATACTA" : "Original cell: <NeuronCell, 2, 2>\nFinal number of cells: 13\nRepertoire: [('T,4', 'UGAUCAUGAUAUGAAUGAAAA', 'MNE'), ('No simple repeats in DNA sequence', 'AUGAUGUAA', 'MM'), ('ACTAT,3;ATACT,3;CTATA,3;TACTA,4;TATAC,3', 'UAGUAUAGUAUAGUAUAGUAGCGAUGAU', 'M')]\nUndergoing meiosis...\nFirst cell genome: ['TTTTCATTCATATCATGATCA', 'TTACATCAT', 'ATCATCGCTACTATACTATACTATACTA']\nFirst cell repertoire: [('T,4', 'UGAUCAUGAUAUGAAUGAAAA', 'MNE'), ('No simple repeats in DNA sequence', 'AUGAUGUAA', 'MM'), ('ACTAT,3;ATACT,3;CTATA,3;TACTA,4;TATAC,3', 'UAGUAUAGUAUAGUAUAGUAGCGAUGAU', 'M')]\nSecond cell genome: ['TGATCATGATATGAATGAAAA', 'ATGATGTAA', 'TAGTATAGTATAGTATAGTAGCGATGAT']\nSecond cell repertoire: [('A,4', 'UUUUCAUUCAUAUCAUGAUCA', 'MI'), ('No simple repeats in DNA sequence', 'UUACAUCAU', 'Non-coding RNA'), ('AGTAT,3;ATAGT,3;GTATA,3;TAGTA,4;TATAG,3', 'AUCAUCGCUACUAUACUAUACUAUACUA', 'Non-coding RNA')]", 
"ProkaryoticCell 100 100 AGCTAGTCAATTGATCTGCATGTTCATGAT" : "Original cell: <ProkaryoticCell, 1, 4>\nFinal number of cells: 100\nRepertoire: [('No simple repeats in DNA sequence', 'AUCAUGAACAUGCAGAUCAAUUGACUAGCU', 'MNMQINULA')]\nCannot undergo meiosis", 
"ProkaryoticCell 5 999 AGCTAGTCAATTGATCTGCATGTTCATGAT ATCGAGTTTGGGTTGTGTGTGAAGTCATC ATCATCGCTACTA" : "Original cell: <ProkaryoticCell, 1, 4> Final number of cells: 997 Repertoire: [('No simple repeats in DNA sequence', 'AUCAUGAACAUGCAGAUCAAUUGACUAGCU', 'MNMQINULA'), ('G,3;GT,3;T,3;TG,4', 'GAUGACUUCACACACAACCCAAACUCGAU', 'MTSHTTQTR'), ('No simple repeats in DNA sequence', 'UAGUAGCGAUGAU', 'M')]\nCannot undergo meiosis"}

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
        if [f for f in result_stderr.replace("'","").replace(" ","").replace("\t", "").split("\n") if f][-1].startswith("AssertionError:"):
            print("PASS: Assertion")
        elif [f for f in result_stderr.replace("'","").replace(" ","").replace("\t", "").split("\n") if f][-1].startswith("AssertionError"):
            print("FAILED: Assertion Must Print Offending Variable")
        else:
            print("RUN ERROR: Assertion Expected")
    else:
        # make sure output is the same
        if "".join(result_stdout.lower().replace("'","").replace(" ","").split()) == "".join(gold_result.lower().replace("'","").replace(" ","").split()):
            print("PASS: Correct Output")
        else:
            # if not - notify
            print("FAILED\n")
            print("\n".join(["Expected Output:", gold_result, "\nGiven Output:", result_stdout]))
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
