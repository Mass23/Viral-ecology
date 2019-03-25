from Bio import SeqIO
import argparse
import gzip
import glob
import time
import subprocess

# Arguments parsing
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--File', help='Raw fastq.gz file, put a all if you want all the fastq files in the directory to be processed', type=str, action = 'store', required = True)
parser.add_argument('-q', '--KmerQual', help='K-mer minimal quality average', type=int, action = 'store', required = True)
parser.add_argument('-k', '--KmerLen', help='K-mer length', type=int, action = 'store', required = True)
parser.add_argument('-m', '--MeanQual', help='Read quality mean minimal value', type=int, action = 'store', required = True)
parser.add_argument('-l', '--LowestQual', help='Read quality lowest minimal value', type=int, action = 'store', required = True)

args = parser.parse_args()

fastq_file = args.File
kmer_qual = args.KmerQual
kmer_len = args.KmerLen
mean_qual = args.MeanQual
low_qual = args.LowestQual

def GetKmerQual(kmer_len, qual_record):
    kmer_qual_list = []

    for i in range(0,len(qual_record), kmer_len):
        kmer = qual_record[i:i+kmer_len+1]
        kmer_quality = mean(kmer)
        kmer_qual_list.append(kmer_quality)

    return(min(kmer_qual_list))

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

# Function to process one file
def FilterReads(file):
    filtered_reads = []
    filtered_count = 0
    kept_count = 0

    with gzip.open(file, 'rt') as f:
        for sequence in SeqIO.parse(f, 'fastq'):
            print('Running:', str(filtered_count + kept_count + 1), 'reads processed' , end='\r')

            lowest_quality = min(sequence.letter_annotations['phred_quality'])
            mean_quality = mean(sequence.letter_annotations['phred_quality'])
            kmer_min_qual = GetKmerQual(kmer_len, sequence.letter_annotations['phred_quality'])

            if kmer_min_qual >= kmer_qual and lowest_quality > low_qual and mean_quality >= mean_qual:
                filtered_reads.append(sequence)
                kept_count += 1
            else:
                filtered_count += 1
    try:
        removed_ratio = (kept_count/(filtered_count + kept_count)) * 100
    except:
        removed_ratio = 0
    return(filtered_reads, removed_ratio)

def main():
    if fastq_file == 'all':

        print('\n', '------------ QC_trim.py ------------', '\n')
        print(time.strftime("%Y-%m-%d %H:%M"), '\n')
        print('------------------------------------')
        print(str(kmer_len) + '-mer quality', '\t', ': ', '\t', '\t', kmer_qual)
        print('Mean quality', '\t', ': ', '\t', '\t', mean_qual)
        print('Lowest quality', '\t', ': ', '\t', '\t', low_qual)
        print('All fastq files in the directory mode!')
        print('------------------------------------', '\n')

        dir_fastq1 = glob.glob('*.fastq.gz')
        dir_fastq2 = glob.glob('*.fq.gz')
        dir_fastq = dir_fastq1 + dir_fastq2
        print('Files to process: ')
        for i in dir_fastq:
            print(' -', i)
        print('\n')

        with open('All_QC_stats.txt', 'w') as out:
            for fastq in dir_fastq:
                name = fastq.replace('.fq.gz','').replace('.fastq.gz','')

                print('File: ', fastq)
                filtered_reads = FilterReads(fastq)
                out_list = filtered_reads[0]
                print(' -', str(round(filtered_reads[1],4)) + '%', 'of reads kept.')

                SeqIO.write(out_list, name + '_trimmed.fastq', 'fastq')
                subprocess.call('gzip' + name + '_trimmed.fastq', shell = True)

                print(' -', fastq, ' done!')

                out.write(fastq)
                out.write('\t')
                out.write(str(filtered_reads[1]))
                out.write('\n')

    else:
        name = fastq_file.replace('.fq.gz','').replace('.fastq.gz','')
        with open(name + '_QC_stats.txt', 'w') as out:

            print('\n', '------------ QC_trim.py ------------', '\n')
            print(time.strftime("%Y-%m-%d %H:%M"), '\n')
            print('------------------------------------')
            print(str(kmer_len) + '-mer quality', '\t', ': ', '\t', '\t', kmer_qual)
            print('Mean quality', '\t', ': ', '\t', '\t', mean_qual)
            print('Lowest quality', '\t', ': ', '\t', '\t', low_qual)
            print('File: ', fastq_file)
            print('------------------------------------', '\n')

            filtered_reads = FilterReads(fastq_file)
            out_list = filtered_reads[0]

            SeqIO.write(out_list, name + '_trimmed.fastq', 'fastq')
            subprocess.call('gzip' + name + '_trimmed.fastq', shell = True)

            print('File:', fastq_file, 'done:')
            print(' -', str(round(filtered_reads[1],4)) + '%', 'of reads kept.', '\n')

            out.write(str(fastq_file))
            out.write('\t')
            out.write(str(filtered_reads[1]))
            out.write('\n')

if __name__== "__main__":
    main()
