import glob
import subprocess
import argpage

parser = argparse.ArgumentParser()

parser.add_argument('-r', '--ReferenceFasta', help='Reference to map the reads against, in fasta format.', type=str, action = 'store', required = True)
parser.add_argument('-n', '--Ncores', help='Number of cores for multiprocessing.', type=str, action = 'store', required = True)

args = parser.parse_args()

reference_fasta = args.ReferenceFasta
n_cores = args.Ncores

def trimmomatic(fastqfile, name, cores):
    r1 = ''
    r2 = ''

    if "_R1" in fastqfile:
        r1 = fastqfile
        r2 = r1.replace("_R1", "_R2")

    elif "_R2" in fastqfile:
        r2 = fastqfile
        r1 = r2.replace("_R2", "_R1")

    trim_args = ["trimmomatic", "PE", "-threads",  cores, "-phred33",
            # Input R1, R2
            r1 , r2,
            # Output forward/reverse, paired/unpaired
            name + "_forward_paired.fq.gz",
            name + "_forward_unpaired.fq.gz",
            name + "_reverse_paired.fq.gz",
            name + "_reverse_unpaired.fq.gz",
            "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10",
            "LEADING:10",
            "TRAILING:10",
            "SLIDINGWINDOW:4:15",
            "MINLEN:36"]

            subprocess.call(' '.join(trim_args), shell = True)

def bwa_map(name, ref, cores):

    # Merge forward and reverse files
    args_forward = ["cat", name + "_forward_paired.fq.gz", ">", name + "_forward.fastq.gz"]
    args_reverse = ["cat", name + "_reverse_paired.fq.gz", ">", name + "_reverse.fastq.gz"]

    subprocess.call(" ".join(args_forward), shell = True)
    subprocess.call(" ".join(args_reverse), shell = True)

    # Map reads
    map_args = ["bwa", "mem", "-T", cores, "-M", "-R", "'" + r"@RG\tID:" + name + r"\tSM:" + name + r"\t'" ,ref, name + "_forward.fastq.gz", name + "_reverse.fastq.gz", ">", name + "_alignment.sam"]
    subprocess.call(" ".join(map_args), shell = True)

def main():
    fastq_list = glob.glob('*_R1.fq.gz')

    for file in fastq_list:
        name = file.replace('_R1.fq.gz','')
        print('Sample:', name)
        trimmomatic(file, name, n_cores)
        bwa_map(name, reference_fasta, n_cores)
        print(' - Done!')

if __name__== "__main__":
    main()
