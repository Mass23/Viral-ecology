#!/bin/bash
import subprocess
import multiprocessing
import argparse
import glob

# python3 preprocess.py -r REF -t 8
parser = argparse.ArgumentParser()

parser.add_argument('-t', '--NumberThreads', help='Number of threads to use', type=str, action = 'store', required = True)
parser.add_argument('-r', '--ReferenceFasta', help='Fasta of the reference file used for mapping', type=int, action = 'store', required = True)

args = parser.parse_args()

n_threads = args.NumberThreads
reference = args.ReferenceFasta

fastq_files = glob.glob('*.fastq.gz')
samples = list(set([i.replace('_R1.fastq.gz','').replace('_R2.fastq.gz','') for i in fastq_files]))

def trimmomatic(sample, threads):
    r1 = ''
    r2 = ''

    if "_R1" in fastqfile:
        r1 = fastqfile
        r2 = r1.replace("_R1", "_R2")

    elif "_R2" in fastqfile:
        r2 = fastqfile
        r1 = r2.replace("_R2", "_R1")

    args = ["trimmomatic", "PE", "-threads",  threads, "-phred33",
            # Input R1, R2
            r1 , r2,
            # Output forward/reverse, paired/unpaired
            individual + "_forward_paired.fq.gz",
            individual + "_forward_unpaired.fq.gz",
            individual + "_reverse_paired.fq.gz",
            individual + "_reverse_unpaired.fq.gz",
            "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10",
            "LEADING:3",
            "TRAILING:3",
            "SLIDINGWINDOW:4:15",
            "MINLEN:36"]

            subprocess.call(' '.join(args), shell = True)

def bwa_map(sample, ref, cores):

    # Merge forward and reverse files
    args_forward = ["cat", individual + "_forward_paired.fq.gz", ">", individual + "_forward.fastq.gz"]
    args_reverse = ["cat", individual + "_reverse_paired.fq.gz", ">", individual + "_reverse.fastq.gz"]

    subprocess.call(" ".join(args_forward), shell = True)
    subprocess.call(" ".join(args_reverse), shell = True)

    # Map reads
    map_args = ["bwa", "mem", "-T", n_threads, "-M", "-R", "'" + r"@RG\tID:" + individual + r"\tSM:" + individual + r"\t'" ,ref, individual + "_forward.fastq.gz", individual + "_reverse.fastq.gz", ">", individual + "_alignment.sam"]
    subprocess.call(" ".join(map_args), shell = True)

for sample in samples_list:
    trimmomatic(sample, n_threads)
    bwa_map(sample, reference, n_threads)
