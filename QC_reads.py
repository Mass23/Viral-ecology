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
    return(record.annotations['taxonomy'][0])

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


############ MAIN ##############################################################

with gzip.open(fastq_file, 'rt') as f:
    f_parsed = list(SeqIO.parse(f, 'fastq'))
    f_length = len(f_parsed)

    out_quality = []
    out_length = []
    out_taxonomy = []
    out_gc = []

    seeds = [random.randint(0, f_length-1) for i in range(0, reads_number)]

    for i in range(0, len(seeds)):
        current_seed = seeds[i]
        
        # Create seed number, get read length
        print('Read number: ', i+1 , '/', reads_number)
        read_length = len(f_parsed[current_seed].seq)
        gc_content = GCcontent(f_parsed[current_seed].seq)

        # Blast sequence on NCBI
        
        try:
            print('  - Blasting...')
            f_blast = NCBIWWW.qblast('blastn', 'nt', f_parsed[current_seed].seq)

            # Search taxonomy using accession ID
            print('  - Searching taxonomy...')
            results = GetEvalueId(f_blast)

            GI_id = results[0]
            evalue = float(results[1])
            f_taxonomy = GetTaxonomyGI(GI_id)
            
        except:
            f_taxonomy = 'Unknown'
            evalue = -1

        # Add output to the lists
        out_quality.append(f_parsed[current_seed].letter_annotations['phred_quality'])
        out_length.append(read_length)

        if evalue < e_value_threshold:
            out_taxonomy.append(f_taxonomy)
            
        if gc_content != -1:
            out_gc.append(gc_content)

    # Get quality data into the right format:
    max_length = 0
    
    for i in out_quality:
        if len(i) > max_length:
            max_length = len(i)

    qual_dict = {}
    for i in range(0, max_length):
        qual_dict[i] = []
        
        for j in out_quality:
            
            try:
                qual_dict[i].append(j[i])
                
            except:
                qual_dict[i].append(np.nan)
                
    sorted_keys = [i for i in range(1, max_length + 1)]
    sorted_vals = tuple([qual_dict[i] for i in range(0, max_length)])

    # Plot the output
    fig, ax = plt.subplots(figsize=(13,10), ncols=2, nrows=2)
    sns.set_style('whitegrid')

    print('Plotting read length...')
    ax1 = sns.distplot(out_length, ax=ax[0][0])
    ax1.set(xlabel="Read length", ylabel="Frequency")

    print('Plotting GC content...')
    ax2 = sns.distplot(out_gc, ax=ax[1][0])
    ax2.set(xlabel="GC content", ylabel="Frequency")

    print('Plotting out_taxonomy...')
    ax3 = sns.countplot(out_taxonomy, ax=ax[0][1])

    print('Plotting quality...')
    ax4 = sns.boxplot(data=sorted_vals, ax=ax[1][1], palette="Blues")
    ax4.set(xlabel="Read position", ylabel="Phred Quality Score")
    plt.xticks(plt.xticks()[0], sorted_keys)

    plt.savefig('QC_' + fastq_file + '_n' + str(reads_number) + '.png', dpi = 1000)
    print('Done!')
    plt.close()
