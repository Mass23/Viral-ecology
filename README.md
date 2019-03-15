# 1. Pipeline

## 1.1 Sort contigs: VirSorter

## 1.2 Mapping

> Reads from the 214 bulk soil metagenomes were quality trimmed using Trimmomatic v0.3635 and then paired reads were mapped to 
> the viral contig database with Bowtie236, using default parameters. The output bam files were passed to BamM ‘filter’ v1.7.2 
> (http://ecogenomics.github.io/BamM/, accessed 15 December 2015) and reads that were aligned over ≥90% of their length at ≥95% > nucleic acid identity were retained.

1. Trimming: trimmomatic
- http://www.usadellab.org/cms/?page=trimmomatic
2. Mapping: BWA
- http://bio-bwa.sourceforge.net/
3. Bam filter: BamM 'filter
- http://ecogenomics.github.io/BamM/

## 1.3 OTUs calling
### 1.3.1 NCBI assemblies of viruses genomes
To retrieve a fasta of all genome assemblies avauilable on NCBI:
- List of accessions: [Assemblies_db_15-03-19.txt](https://github.com/Mass23/Viral-ecology/blob/master/Assemblies_db_15-03-19.txt)
- Script to create fasta file: [viral_assemblies.py](https://github.com/Mass23/Viral-ecology/blob/master/viral_assemblies.py)

A few downloads failed:
- List of accessions: [stats_db.txt](https://github.com/Mass23/Viral-ecology/blob/master/stats_db.txt)

### 1.3.2 Merge with contigs

### 1.3.3 Search OTUs

## 1.4 Mapping

## 1.5 Filtering

**References:**
- https://www.nature.com/articles/s41564-018-0190-y#ref-CR18
