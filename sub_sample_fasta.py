import argparse
from Bio import SeqIO
import random

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--FastaFile', help='Fasta file', type=str, action = 'store', required = True)
parser.add_argument('-o', '--OutputFile', help='Fasta file', type=str, action = 'store', required = True)
parser.add_argument('-n', '--NumberOfSequences', help='Number of sequences to sub-sample from the fasta file.', type=int, action = 'store', required = False)
parser.add_argument('-p', '--ProportionOfSequences', help='Proportion of sequences to sub-sample from the fasta file.', type=int, action = 'store', required = False)

args = parser.parse_args()

fasta_file = args.FastaFile
output_file = args.OutputFile
seq_number = args.NumberOfSequences
seq_proportion = args.ProportionOfSequences

def ParseFasta(file):
    sequences_list = []
    for sequence in SeqIO.parse(file, 'fasta'):
        sequences_list.append(sequence)
    length_fasta = len(sequences_list)
    return(sequences_list, length_fasta)

def ProportionMethod(sequences_list, proportion, length_fasta):
    print('Number of sequences: ', length_fasta)
    number = int(length_fasta * proportion)
    sub_sample = random.sample(sequences_list, number)
    return(sub_sample)

def NumberMethod(sequences_list, number, length_fasta):
    if number > length_fasta:
        print('Cannot sample more sequences than present in the fasta file, exiting...')
        return(None)
    else:
        sub_sample = random.sample(sequences_list, number)
        return(sub_sample)

def SequenceListToFasta(sequence_list, output_file):
    SeqIO.write(sequence_list, output_file, 'fasta')
    print('Fasta file written!')

def main():
    print('--------- rarefy.py ---------', '\n')
    data = ParseFasta(fasta_file)
    name = fasta_file.replace('.fasta','').replace('.fa','')
    sequences = data[0]
    length = data[1]

    if seq_number == None and seq_proportion == None:
        print('Please use an option, -n or -p!')

    elif seq_number == None:
        print('Sub-sampling using the proportion mode:')

        number_kept = int(length_fasta * proportion)
        print('Keeping', seq_proportion + '%','(' + number_kept + ')', 'of the sequences.', '\n')

        sub_sample_proportion = ProportionMethod(fasta_file, seq_proportion, length)
        SequenceListToFasta(sub_sample_proportion)

        print('Done!')

    elif seq_proportion == None:
        print('Sub-sampling using the number mode:')

        proportion_kept = seq_number / length
        print('Keeping', proportion_kept + '%','(' + seq_number + ')', 'of the sequences.', '\n')

        sub_sample_numeric = NumberMethod(fasta_file, seq_number, length)
        SequenceListToFasta(sub_sample_numeric)

        print('Done!')

    else:
        print('You cannot use both proportional and numerical methods!')

if __name__== "__main__":
    main()
