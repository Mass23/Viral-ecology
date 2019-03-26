# 1. OTUs calling
## 1.1 Database of assembled viromes

Two datasets were used:
- To retrieve a fasta of all genome assemblies available on NCBI:
    - List of accessions: [Assemblies_db_15-03-19.txt](https://github.com/Mass23/Viral-ecology/blob/master/Assemblies_db_15-03-19.txt)
    - Script to create fasta file: [viral_assemblies.py](https://github.com/Mass23/Viral-ecology/blob/master/viral_assemblies.py)

- Earth's virome database: https://www.nature.com/articles/nature19094#methods

## 1.2 Merge with contigs

```
cat /home/fodelian/Desktop/ViralGenomes/assembly_db/refseq_viral_genomes.fasta \
    /home/fodelian/Desktop/ViralGenomes/assembly_db/mVGs_sequences_v2.fasta  \
    /home/fodelian/Desktop/ViralGenomes/SNG/SNG_contigs.fasta  \
    /home/fodelian/Desktop/ViralGenomes/VDN/VDN_contigs.fasta  \
    /home/fodelian/Desktop/ViralGenomes/VEV/VEV_contigs.fasta  \
     > raw_db_ctgs.fasta
```

### 1.3 Find 95% identity centroids

Vsearch: https://github.com/torognes/vsearch

```
vsearch --cluster_fast raw_db_ctgs.fasta --consout 95_database.fasta --id 0.95 --iddef 0 --maxseqlength 3000000 --threads 6 --usersort
```

## 1.4 Mapping

Script: [preprocess.py](https://github.com/Mass23/Viral-ecology/blob/master/preprocess.py)

> Reads from the 214 bulk soil metagenomes were quality trimmed using Trimmomatic v0.3635 and then paired reads were mapped to 
> the viral contig database with Bowtie236, using default parameters. The output bam files were passed to BamM ‘filter’ v1.7.2 
> (http://ecogenomics.github.io/BamM/, accessed 15 December 2015) and reads that were aligned over ≥90% of their length at ≥95% > nucleic acid identity were retained.

First, we need to merge the files per sample:
```
cat SNG1_R1.fq.gz SNG2_R1.fq.gz > SNG_R1.fq.gz
cat VDN1_R1.fq.gz VDN2_R1.fq.gz > VDN_R1.fq.gz
cat VEV1_R1.fq.gz VEV2_R1.fq.gz > VEV_R1.fq.gz

cat SNG1_R2.fq.gz SNG2_R2.fq.gz > SNG_R2.fq.gz
cat VDN1_R2.fq.gz VDN2_R2.fq.gz > VDN_R2.fq.gz
cat VEV1_R2.fq.gz VEV2_R2.fq.gz > VEV_R2.fq.gz
```

1. Trimming: Trimmomatic
- http://www.usadellab.org/cms/?page=trimmomatic

2. Mapping: BWA
- http://bio-bwa.sourceforge.net/

## 1.5 Filtering

Bam filter: BamM 'filter
- http://ecogenomics.github.io/BamM/


# 2. Analysis

**References:**
- https://www.nature.com/articles/s41564-018-0190-y#ref-CR18
