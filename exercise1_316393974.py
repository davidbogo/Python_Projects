# David Bogoslavsky 316393974
import sys


def find_srr(dna_seq):
    ret_string = None
    for short_seq_len in range(1, 7):  # we begin with sequences of 1 and end with 6
        for start_index in range(0, len(dna_seq) - short_seq_len):
            # we run through every sequence from the first letter to the last
            short_seq = dna_seq[start_index:start_index + short_seq_len]  # the sequence we compare the next ones
            counter = 1
            for i in range(start_index + short_seq_len, len(dna_seq) + 1 - short_seq_len, short_seq_len):
                cur_seq = dna_seq[i:i + short_seq_len]  # the next sequences after short seq
                if short_seq == cur_seq:  # we check whether the next sequence is equal to the previous
                    counter += 1  # we add 1 to the counter if so
                if short_seq != cur_seq or i + short_seq_len * 2 > len(dna_seq):
                    # we've reached the final count. Either we had a mismatch or the string is over
                    if counter >= 3:  # if our next sequence doesn't equal to before, we check 3 or more for printing
                        if ret_string:
                            ret_string += ";"
                        else:
                            ret_string = ""
                        ret_string += short_seq
                        ret_string += ","
                        ret_string += str(counter)
                    break
    return ret_string


def transcribe(dna_seq):
    dna = dna_seq.upper()
    rna = ''  # initialize RNA to empty string
    for letter in dna:  # loop over DNA sequence characters to make them complementary nucleotide
        if letter == 'A':
            rna += 'U'
        elif letter == 'T':
            rna += 'A'
        elif letter == 'C':
            rna += 'G'
        elif letter == 'G':
            rna += 'C'
        else:
            rna += 'N'
    rna = rna[::-1]     # reverse the string
    return rna  # return output RNA sequence


def translate(rna_seq):
    rna = rna_seq.upper()
    ret = None
    rna_longest_seq = ''
    first_printed = False
    for cur_shift in range(0, 3):  # we check all 3 possible frame shifts
        remaining_sequence = rna
        while len(remaining_sequence) > len(rna_longest_seq):  # we do this loop in case we have an end codon in frame
            # we continue checking for longer proteins after the current long seq
            counter = 3   # we take the AUG in consideration
            while True:
                idx = remaining_sequence.find('AUG')  # find start index of Methionine codon in rna
                if idx == -1:  # if not found
                    break
                if idx % 3 == cur_shift:  # it means that our found AUG is in our frame shift
                    break
                remaining_sequence = remaining_sequence[idx + 3:]
                # if the AUG we found is not in our frameshift we check the next remaining codons
            if idx == -1:  # if we haven't found AUG, we break the loop, because what is next is not relevant
                break
            end_codon_found = False
            if len(rna) == 3:  # a case for a single protein which is AUG only
                rna_longest_seq = 'AUG'
            for codon_start_index in range(idx, len(remaining_sequence) - 2, 3):
                # we check and create a sequence of nucleotides
                cur_codon = remaining_sequence[codon_start_index: codon_start_index + 3]
                counter += 3
                if counter > len(rna_longest_seq):
                    rna_longest_seq = remaining_sequence[idx: codon_start_index + 3]
                    # we change the biggest rna seq if we find that other frame shift has a longer one
                if cur_codon == 'UAG' or cur_codon == 'UAA' or cur_codon == 'UGA':  # end codons
                    rna_longest_seq = remaining_sequence[idx: codon_start_index]  # we don't include the end codon
                    end_codon_found = True  # declare we found end codon in frame
                    remaining_sequence = remaining_sequence[codon_start_index + 3 - cur_shift:]  # check seq after end
                    break
            if end_codon_found:  # if we found end codon, we continue looping in the same frame, after the end codon
                pass
            else:  # if else, we move to the next frame
                break
    if len(rna_longest_seq) >= 3:
        sent_rna = ''
        for i in range(0, len(rna_longest_seq), 3):
            next_codon = rna_longest_seq[i:i + 3]
            if len(next_codon) != 3:  # we cut remaining nucleotides if they are not 3
                break
            if first_printed:  # we order the nucleotides by codons as we have been instructed
                sent_rna += ';'
            sent_rna += next_codon
            first_printed = True
        ret = sent_rna
    return ret


def main():
    srr = find_srr(sys.argv[1])
    if srr:
        print(srr)
    else:
        print('No simple repeats in DNA sequence')
    transcribed_rna = transcribe(sys.argv[1])
    print('RNA sequence: ' + transcribed_rna)
    translated = translate(transcribed_rna)
    if translated:
        print('Translation: ' + translated)
    else:
        print('Non-coding RNA')


if __name__ == "__main__":
    main()
