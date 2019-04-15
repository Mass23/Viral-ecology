import glob
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from numpy import mean as mean

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
fastq_list = glob.glob('*.fq.gz')

with open('reads_metrics.tsv', 'w') as out:
    out.write('file\tgc_content\treads_count\n')

    for fastq_file in fastq_list:
        with gzip.open(fastq_file, 'rt') as f:

            gc_content=[]
            count = 0

            for title, seq, qual in FastqGeneralIterator(f) :
                count += 1
                gc_content.append(GCcontent(seq))

            out.write(fastq_file + '\t' + str(mean(gc_content)) + '\t' + str(count) + '\n')
            print(fastq_file + 'Done!')
