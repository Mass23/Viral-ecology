import os
import gzip
import argparse
import random

from Bio import SeqIO

import pandas as pd
import numpy as np

import seaborn as sns
from matplotlib import pyplot as plt
from Bio.SeqIO.QualityIO import FastqGeneralIterator

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--File', help='File, raw fastq.gz', type=str, action = 'store', required = True)
parser.add_argument('-n', '--ReadsNumber', help='Number of reads to sub-sample', type=int, action = 'store', required = True)

args = parser.parse_args()

fastq_file = args.File
reads_number = args.ReadsNumber

############ FUNCTIONS #########################################################
def CharToPhred(char):
    dict_phred = {'!': 0, '\"': 1, '#': 2, '$': 3, '&': 5, "\'": 6, "(": 7,
                  ")": 8 , "*": 9, "+": 10, ",": 11, "-": 12, ".": 13, "/": 14,
                  "0": 15, "1": 16, "2": 17, "3": 18, "4": 19, "5": 20, "6": 21,
                  "7": 22, "8": 23, "9": 24, ":": 25, ";": 26, "<": 27, "=": 28,
                  ">": 29, "?": 30, "@": 31, "A": 32, "B": 33, "C": 34, "D": 35,
                  "E": 36, "F": 37, "G": 38, "H": 39, "I": 40}
    return(dict_phred[char])

def GCcontent(seq):
    g_count = seq.count('G')
    c_count = seq.count('C')
    a_count = seq.count('A')
    t_count = seq.count('T')
    try:
        gc_content = (g_count + c_count) / (g_count + c_count + a_count + t_count)
    except:
        gc_content = 'NA'
    return(gc_content)

################################################################################
with open(fastq_file.replace('.fastq.gz','').replace('.fq.gz','') + '_qc.tsv', 'w') as out:
    with gzip.open(fastq_file, 'rt') as f:

        print('Parsing fastq file...')
        data = list()

        for title, seq, qual in FastqGeneralIterator(f) :
                data.append([title, seq, qual])

        sub_sample = random.sample(data, reads_number)
        count = 0

        max_read_length = min([len(i[1]) for i in sub_sample])
        out.write('read_number\tlength\tgc_content\t' + '\t'.join([str(i) for i in range(1,max_read_length+1)]) + '\n')

        print('Calculating metrics for the subsample...')
        for i in sub_sample:
            count += 1
            print(' - ', str(count) + ' / ' + str(reads_number) + '\r', end = '')

            length = len(i[1])
            gc_content = GCcontent(i[1])
            quality = i[2]

            out.write(str(count) + '\t' + str(length) + '\t' + str(gc_content) + '\t' + '\t'.join([str(CharToPhred(i)) for i in quality]) + '\n')

        print('Output file done!')

with open(fastq_file.replace('.fastq.gz','').replace('.fq.gz','') + '_qc.tsv', 'r') as f:
    data = pd.read_csv(f, sep = '\t')

    plt.figure(figsize=(7,5))
    sns.barplot(data['length'])
    plt.savefig(fastq_file.replace('.fastq.gz','').replace('.fq.gz','') + '_length.png', dpi=800)
    plt.xlabel('Length (bp)')

    plt.figure(figsize=(7,5))
    sns.distplot(data['gc_content'])
    plt.savefig(fastq_file.replace('.fastq.gz','').replace('.fq.gz','') + '_gc.png', dpi=800)
    plt.xlabel('GC content')
    plt.ylabel('Frequency')

    plt.figure(figsize=(7,5))
    quality = data[[str(i) for i in range(1,max_read_length+1)]]
    quality = pd.melt(quality).astype('int').sort_values(by=['variable'])
    sns.boxplot(x='variable', y='value',data=quality, fliersize = 0, linewidth=0.6)
    plt.xlabel('Position (bp)')
    plt.ylabel('Quality (Phred score)')
    plt.xticks([0,50,100,150,200,250,300],[0,50,100,150,200,250,300])
    plt.xlim(0,max_read_length)
    plt.savefig(fastq_file.replace('.fastq.gz','').replace('.fq.gz','') + '_quality.png', dpi=800)

    print('Figures saved!')
