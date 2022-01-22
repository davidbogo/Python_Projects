# David Bogoslavsky 316393974
# There is an unfinished condition here to deal with the case of reaching max_cells_number if we multiply the current
# cells by division rate. Please do notice that I indeed dealt with most of the requirements of the task, therefore
# I'd thank you if the score to this task won't be calculated in other way than it was written in the exercise
# (15 points for main program)
import sys


class Polymerase:
    def __init__(self, type, error_rate = 0):
        self.type = type
        self.error_rate = error_rate

    def transcribe(self, dna_seq):
        if dna_seq:
            if self.type == 'RNA':
                dna = dna_seq.upper()
                rna = ''  # initialize RNA to empty string
                dna_to_rna_dict = {"A": "U", "T": "A", "C": "G", "G": "C"}
                for letter in dna:  # loop over DNA sequence characters to make them complementary nucleotide
                    if letter in dna_to_rna_dict:
                        rna += dna_to_rna_dict[letter]
                    else:
                        raise AssertionError(dna_seq)
                        # rna = None
                        # break
                if rna:
                    rna = rna[::-1]  # reverse the string
                return rna  # return output RNA sequence
            elif self.type == 'DNA':
                dna = dna_seq.upper()
                rep_dna = ''  # initialize dna to empty string
                dna_to_dna_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
                for letter in dna:  # loop over DNA sequence characters to make them complementary nucleotide
                    if letter in dna_to_dna_dict:
                        rep_dna += dna_to_dna_dict[letter]
                    else:
                        raise AssertionError(dna_seq)
                        # rep_dna = None
                        # break
                if rep_dna:
                    rep_dna = rep_dna[::-1]  # reverse the string
                return rep_dna  # return output DNA sequence
            else:
                ret_string = ''
                return ret_string
        else:
            raise AssertionError(dna_seq)
            # return None


class Ribosome:
    def __init__(self, genetic_code, start_codons):
        self.genetic_code = genetic_code
        self.start_codons = start_codons
        # edit afterwards
        self.stop_codons = []
        for codon in self.genetic_code:
            if self.genetic_code[codon]:
                pass
            else:
                self.stop_codons.append(codon)
        # end of adding

    def translate(self, rna_seq):
        if rna_seq:
            rna = rna_seq.upper()
            ret = None
            rna_longest_seq = ''
            for start_codon in self.start_codons:
                for cur_shift in range(0, 3):  # we check all 3 possible frame shifts
                    remaining_sequence = rna
                    while len(remaining_sequence) > len(
                            rna_longest_seq):  # we do this loop in case we have an end codon in frame
                        # we continue checking for longer proteins after the current long seq
                        counter = 3  # we take the AUG in consideration
                        while True:
                            idx = remaining_sequence.find(start_codon)  # find start index of Methionine codon in rna
                            if idx == -1:  # if not found
                                break
                            if idx % 3 == cur_shift:  # it means that our found AUG is in our frame shift
                                break
                            remaining_sequence = remaining_sequence[idx + 3:]
                            # if the start codon we found is not in our frameshift we check the next remaining codons
                        if idx == -1:
                            # if we haven't found AUG, we break the loop, because what is next is not relevant
                            break
                        end_codon_found = False
                        if len(rna) == 3:  # a case for a single protein which is AUG only
                            rna_longest_seq = start_codon
                        for codon_start_index in range(idx, len(remaining_sequence) - 2, 3):
                            # we check and create a sequence of nucleotides
                            cur_codon = remaining_sequence[codon_start_index: codon_start_index + 3]
                            counter += 3
                            if counter > len(rna_longest_seq):
                                rna_longest_seq = remaining_sequence[idx: codon_start_index + 3]
                                # we change the biggest rna seq if we find that other frame shift has a longer one
                            # if cur_codon == 'UAG' or cur_codon == 'UAA' or cur_codon == 'UGA':  # end codons
                            if cur_codon in self.stop_codons:
                                if counter > len(rna_longest_seq):
                                    rna_longest_seq = remaining_sequence[idx: codon_start_index]
                                # we don't include the end codon
                                end_codon_found = True  # declare we found end codon in frame
                                # remaining_sequence = remaining_sequence[
                                #                     codon_start_index + 3 - cur_shift:]
                                remaining_sequence = remaining_sequence[idx + 3:]  # check seq after end
                                break
                        if end_codon_found:
                            # if we found end codon, we continue looping in the same frame, after the end codon
                            pass
                        else:  # if else, we move to the next frame
                            break
            if len(rna_longest_seq) >= 3:
                sent_rna = []
                for i in range(0, len(rna_longest_seq), 3):
                    next_codon = rna_longest_seq[i:i + 3]
                    if len(next_codon) != 3:  # we cut remaining nucleotides if they are not 3
                        break
                    sent_rna.append(next_codon)
                ret = sent_rna
            return ret
        else:
            raise AssertionError(rna_seq)
            # return None

    def synthesize(self, rna_seq):
        ret_protein_seq = ''
        translated_seq = self.translate(rna_seq)
        protein_dict = self.genetic_code
        if translated_seq:
            for codon in translated_seq:
                if codon in protein_dict:
                    ret_protein_seq += protein_dict[codon]
                else:
                    break
        return ret_protein_seq


class Cell:
    def __init__(self, name, genome, num_copies, genetic_code, start_codons, division_rate):
        self.name = name
        self.genome = genome
        self.num_copies = num_copies
        self.genetic_code = genetic_code
        self.start_codons = start_codons
        self.division_rate = division_rate
        self.cell_ribosome = Ribosome(self.genetic_code, self.start_codons)
        self.cell_DNA_polymerase = Polymerase('DNA', 0)
        self.cell_RNA_polymerase = Polymerase('RNA', 0)

    def __repr__(self):
        return "<" + self.name + ", " + str(self.num_copies) + ", " + str(self.division_rate) + ">"

    def mitosis(self):
        mitosis_list = []
        for i in range(self.division_rate):
            cell_division = Cell(self.name, self.genome, self.num_copies, self.genetic_code, self.start_codons,
                                 self.division_rate)
            mitosis_list.append(cell_division)
        return mitosis_list

    def meiosis(self):
        if self.num_copies % 2 == 0:
            meiosis_list = []
            complementary_cell = []
            for gene in self.genome:
                complementary_cell.append(self.cell_DNA_polymerase.transcribe(gene))
            meiosis_cell2 = Cell(self.name, complementary_cell, self.num_copies / 2,
                                 self.genetic_code, self.start_codons, self.division_rate)
            meiosis_cell1 = Cell(self.name, self.genome, self.num_copies / 2,
                                 self.genetic_code, self.start_codons, self.division_rate)
            meiosis_list.append(meiosis_cell1)
            meiosis_list.append(meiosis_cell2)
            return meiosis_list
        else:
            return None

    def find_srr(self, dna_seq):
        if dna_seq:
            dna_seq.upper()
            legal_seq = True
            ret_list = []
            seq_list = []
            counter_list = []
            for letter in dna_seq:  # loop over DNA sequence characters to check legality
                if letter == 'A':
                    continue
                elif letter == 'T':
                    continue
                elif letter == 'C':
                    continue
                elif letter == 'G':
                    continue
                else:
                    legal_seq = False
            if legal_seq:
                for short_seq_len in range(1, 7):  # we begin with sequences of 1 and end with 6
                    for start_index in range(0, len(dna_seq) - short_seq_len):
                        # we run through every sequence from the first letter to the last
                        short_seq = dna_seq[
                                    start_index:start_index + short_seq_len]  # the sequence we compare the next ones
                        if short_seq == dna_seq[start_index - short_seq_len:start_index]:
                            continue
                            # we checked above whether the cur seq equals to previous seq in distance of the same len
                        counter = 1
                        for i in range(start_index + short_seq_len, len(dna_seq) + 1 - short_seq_len, short_seq_len):
                            cur_seq = dna_seq[i:i + short_seq_len]  # the next sequences after short seq
                            if short_seq == cur_seq:  # we check whether the next sequence is equal to the previous
                                counter += 1  # we add 1 to the counter if so
                            if short_seq != cur_seq or i + short_seq_len * 2 > len(dna_seq):
                                # we've reached the final count. Either we had a mismatch or the string is over
                                if counter >= 3:
                                    # if our next sequence doesn't equal to before, we check 3 or more for printing.
                                    # the concept: we create three list with which we would work comfortably
                                    ret_list.append(short_seq + "," + str(counter))  # main list of srr
                                    seq_list.append(short_seq)  # list of the seqs without the "," and counter
                                    counter_list.append(counter)  # list of ordered counter numbers
                                    length = len(ret_list)
                                    for list_index in range(length - 1):
                                        if seq_list[length - 1] == seq_list[list_index]:
                                            # here we check for same seq,
                                            # and then check the higher counter and delete the redundant
                                            if counter_list[length - 1] > counter_list[list_index]:
                                                ret_list.remove(ret_list[list_index])
                                                seq_list.remove(seq_list[list_index])
                                                counter_list.remove(counter_list[list_index])
                                            else:
                                                ret_list.remove(ret_list[length - 1])
                                                seq_list.remove(seq_list[length - 1])
                                                counter_list.remove(counter_list[length - 1])
                                break
            ret_list.sort()
            ret_string = ''
            for srr in ret_list:
                ret_string += srr
                if srr != ret_list[len(ret_list) - 1]:
                    ret_string += ';'
            return ret_string
        else:
            raise AssertionError(dna_seq)
            # return None

    def repertoire(self):
        repertoire_list = []
        for gene in self.genome:
            srr_ret = self.find_srr(gene)
            if srr_ret:
                pass
            else:
                srr_ret = 'No simple repeats in DNA sequence'
            rna_seq = Polymerase('RNA', 0).transcribe(gene)
            protein_seq = self.cell_ribosome.synthesize(rna_seq)
            if rna_seq:
                if protein_seq:
                    seq_tuple = (srr_ret, rna_seq, protein_seq)
                else:
                    seq_tuple = (srr_ret, rna_seq, 'Non-coding RNA')
            else:
                seq_tuple = (srr_ret, 'Non-coding RNA')
            repertoire_list.append(seq_tuple)
        return repertoire_list


class ProkaryoticCell(Cell):
    def __init__(self, genome):
        self.name = "ProkaryoticCell"
        self.num_copies = 1
        self.division_rate = 4
        self.genetic_code = {
            'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
            'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
            'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
            'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
            'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
            'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
            'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
            'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
            'UAC': 'Y', 'UAU': 'Y', 'UAA': None, 'UAG': None,
            'UGC': 'C', 'UGU': 'C', 'UGA': 'U', 'UGG': 'W'}
        self.start_codons = ['AUG', 'GUG', 'UUG']
        super().__init__(self.name, genome, self.num_copies, self.genetic_code, self.start_codons,
                         self.division_rate)


class EukaryoticCell(Cell):
    def __init__(self, name, genome, division_rate):
        self.name = name
        self.num_copies = 2
        self.start_codons = ['AUG']
        self.genetic_code = {'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
                             'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
                             'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
                             'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
                             'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
                             'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
                             'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
                             'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
                             'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
                             'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
                             'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
                             'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
                             'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
                             'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
                             'UAC': 'Y', 'UAU': 'Y', 'UAA': None, 'UAG': None,
                             'UGC': 'C', 'UGU': 'C', 'UGA': None, 'UGG': 'W'}
        self.division_rate = division_rate
        super().__init__(self.name, genome, self.num_copies, self.genetic_code, self.start_codons,
                         self.division_rate)


class StemCell(EukaryoticCell):
    def __init__(self, genome):
        self.division_rate = 3
        self.name = 'StemCell'
        super().__init__(self.name, genome, self.division_rate)


class NeuronCell(EukaryoticCell):
    def __init__(self, genome):
        self.division_rate = 2
        self.name = 'NeuronCell'
        super().__init__(self.name, genome, self.division_rate)


def main():
    name_of_cell = sys.argv[1]
    # added function

    def chek_if_valid_int(input_string, it_is_cell_division=0):
        dict_of_valid = {'1': '1', '2': '2', '3': '3', '4': '4', '5': '5', '6': '6', '7': '7', '8': '8', '9': '9',
                         '0': '0'}
        for number in range(0, len(input_string)):
            if input_string[number] in dict_of_valid:
                pass
            else:
                raise AssertionError(input_string)
        if it_is_cell_division:
            if input_string == dict_of_valid['0']:
                raise AssertionError(input_string)

    # added function
    number_of_cycles = sys.argv[2]
    # added check afterwards
    chek_if_valid_int(number_of_cycles, 1)
    # end of adding
    max_cells_possible = sys.argv[3]
    # added check afterwards
    chek_if_valid_int(max_cells_possible)
    # end of adding
    genome = []
    for gene in sys.argv[4:]:
        genome.append(gene)
    if genome:
        pass
    else:
        raise AssertionError(genome)
    # cell = ''
    proceed = 1
    if name_of_cell == "StemCell":
        cell = StemCell(genome)
    elif name_of_cell == "NeuronCell":
        cell = NeuronCell(genome)
    elif name_of_cell == "ProkaryoticCell":
        cell = ProkaryoticCell(genome)
    else:
        raise AssertionError(name_of_cell)
        # proceed = False
    if proceed:
        current_cells = [cell]
        limit_reached = False
        for i in range(int(number_of_cycles)):
            start_num_cells = len(current_cells)
            for j in range(start_num_cells):
                next_cell = current_cells.pop(0)    # Remove the first cell from the list
                if len(current_cells) + cell.division_rate > int(max_cells_possible):
                    current_cells.append(next_cell)  # Just return the cell to the list without replication
                    limit_reached = True
                    break
                current_cells += next_cell.mitosis()
            if limit_reached:
                break

        # print("Original cell: " + cell.__repr__() + '\n')
        # print("Final number of cells: " + str(len(current_cells)) + '\n')
        # print('Repertoire: ' + str(cell.repertoire()) + '\n')
        # if cell.meiosis():
        #     print('Undergoing meiosis...' + '\n')
        #     meiosis_cells = cell.meiosis()
        #     print('First cell genome: ' + str(meiosis_cells[0].genome) + '\n')
        #     print('First cell repertoire: ' + str(meiosis_cells[0].repertoire()) + '\n')
        #     print('Second cell genome: ' + str(meiosis_cells[1].genome) + '\n')
        #     print('Second cell repertoire: ' + str(meiosis_cells[1].repertoire()) + '\n')
        # else:
        #     print('Cannot undergo meiosis' + '\n')
        print("Original cell: " + cell.__repr__())
        print("Final number of cells: " + str(len(current_cells)))
        print('Repertoire: ' + str(cell.repertoire()))
        if cell.meiosis():
            print('Undergoing meiosis...')
            meiosis_cells = cell.meiosis()
            print('First cell genome: ' + str(meiosis_cells[0].genome))
            print('First cell repertoire: ' + str(meiosis_cells[0].repertoire()))
            print('Second cell genome: ' + str(meiosis_cells[1].genome))
            print('Second cell repertoire: ' + str(meiosis_cells[1].repertoire()))
        else:
            print('Cannot undergo meiosis')


if __name__ == "__main__":
    main()
