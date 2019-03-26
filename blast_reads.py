import os
import gzip
import argparse
import random

from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez

import pandas as pd
import numpy as np
import operator as op

import seaborn as sns
from matplotlib import pyplot as plt

Entrez.email = "massimo.bourquin@epfl.ch"

parser = argparse.ArgumentParser()

parser.add_argument('-f', '--File', help='File, raw fastq.gz', type=str, action = 'store', required = True)
parser.add_argument('-n', '--ReadsNumber', help='Number of reads to sub-sample', type=int, action = 'store', required = True)
parser.add_argument('-e', '--EvalueThreshold', help='E-value threshold, to keep the hit', type=float, action = 'store', required = True)

args = parser.parse_args()
fastq_file = args.File
reads_number = args.ReadsNumber
e_value_threshold = args.EvalueThreshold

############ FUNCTIONS #########################################################
def GetEvalueId(xml_file):
    for i in xml_file.readlines():
        if '<Hit_id>' in i:
            GI_id = i.split('gi|')[1].split('|')[0]
        elif '<Hsp_evalue>' in i:
            evalue = i.split('<Hsp_evalue>')[1].split('</Hsp_evalue>')[0]
            return(GI_id, evalue)
            break

def GetTaxonomyGI(accession):
    handle = Entrez.efetch(db="nucleotide", id=accession,rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    return(record.annotations['taxonomy'])

def GCcontent(seq):
    g_count = seq.count('G')
    c_count = seq.count('C')
    a_count = seq.count('A')
    t_count = seq.count('T')
    try:
        gc_content = (g_count + c_count) / (g_count + c_count + a_count + t_count)
    except:
        gc_content = -1
    return(gc_content)

################################################################################


# Plot quality, n reads per file on histogram
# Plot blast matches, n reads per file
with open(fastq_file.replace('.fastq.gz','').replace('.fq.gz','') + '.tsv', 'w') as out:
    out.write('read_number\tlength\tgc_content\te_value\tmean_qual\n')
    with gzip.open(fastq_file, 'rt') as f:

        for i, sequence in enumerate(SeqIO.parse(f, 'fastq')):

            if i + 1 > reads_number:
                break
            else:
                # Create seed number, get read length
                print('Read number: ', i+1 , '/', reads_number)
                read_length = len(sequence.seq)
                gc_content = GCcontent(sequence.seq)

                # Blast sequence on NCBI
                try:
                    print('  - Blasting...')
                    f_blast = NCBIWWW.qblast('blastn', 'nt', sequence.seq)

                    # Search taxonomy using accession ID
                    print('  - Searching taxonomy...')
                    results = GetEvalueId(f_blast)

                    GI_id = results[0]
                    evalue = float(results[1])
                    f_taxonomy = GetTaxonomyGI(GI_id)
                except:
                    f_taxonomy = 'Unknown'
                    evalue = -1

                # Add output to the file
                out.write(str(i+1) + '\t' + str(read_length) + '\t' + str(gc_content) + '\t' + evalue + '\t' + mean(sequence.letter_annotations['phred_quality']) + '\n')
