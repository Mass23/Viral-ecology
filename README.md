# 1. Pipeline

## 1.1 Sort contigs: VirSorter

Reference: https://github.com/simroux/VirSorter

```
perl virsorter_db/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl -f SNG/SNG_contigs.fasta --ncpu 3 -d SNG --data-dir virsorter_db/virsorter-data --virome

perl virsorter_db/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl -f VDN/VDN_contigs.fasta --ncpu 3 -d VDN --data-dir virsorter_db/virsorter-data --virome

perl virsorter_db/VirSorter/wrapper_phage_contigs_sorter_iPlant.pl -f VEV/VEV_contigs.fasta --ncpu 3 -d VEV --data-dir virsorter_db/virsorter-data --virome
```

## 1.3 OTUs calling
### 1.3.1 NCBI assemblies of viruses genomes
To retrieve a fasta of all genome assemblies avauilable on NCBI:
- List of accessions: [Assemblies_db_15-03-19.txt](https://github.com/Mass23/Viral-ecology/blob/master/Assemblies_db_15-03-19.txt)
- Script to create fasta file: [viral_assemblies.py](https://github.com/Mass23/Viral-ecology/blob/master/viral_assemblies.py)

Earth's virome database: https://www.nature.com/articles/nature19094#methods

### 1.3.2 Merge with contigs

```
cat /home/fodelian/Desktop/ViralGenomes/assembly_db/refseq_viral_genomes.fasta \
    /home/fodelian/Desktop/ViralGenomes/assembly_db/mVGs_sequences_v2.fasta  \
    /home/fodelian/Desktop/ViralGenomes/SNG/virsorter-out/Predicted_viral_sequences/SNG_cat-1.fasta \
    /home/fodelian/Desktop/ViralGenomes/SNG/virsorter-out/Predicted_viral_sequences/SNG_cat-2.fasta \
    /home/fodelian/Desktop/ViralGenomes/SNG/virsorter-out/Predicted_viral_sequences/SNG_prophages_cat-4.fasta \
    /home/fodelian/Desktop/ViralGenomes/SNG/virsorter-out/Predicted_viral_sequences/SNG_prophages_cat-5.fasta \
    /home/fodelian/Desktop/ViralGenomes/VDN/virsorter-out/Predicted_viral_sequences/VDN_cat-1.fasta \
    /home/fodelian/Desktop/ViralGenomes/VDN/virsorter-out/Predicted_viral_sequences/VDN_cat-2.fasta \
    /home/fodelian/Desktop/ViralGenomes/VDN/virsorter-out/Predicted_viral_sequences/VDN_prophages_cat-4.fasta \
    /home/fodelian/Desktop/ViralGenomes/VDN/virsorter-out/Predicted_viral_sequences/VDN_prophages_cat-5.fasta \
    /home/fodelian/Desktop/ViralGenomes/VEV/virsorter-out/Predicted_viral_sequences/VEV_cat-1.fasta \
    /home/fodelian/Desktop/ViralGenomes/VEV/virsorter-out/Predicted_viral_sequences/VEV_cat-2.fasta \
    /home/fodelian/Desktop/ViralGenomes/VEV/virsorter-out/Predicted_viral_sequences/VEV_prophages_cat-4.fasta \
    /home/fodelian/Desktop/ViralGenomes/VEV/virsorter-out/Predicted_viral_sequences/VEV_prophages_cat-5.fasta \
     > raw_db_ctgs.fasta
```

### 1.3.3 Search 95% OTUs

```
cd-hit-est -i raw_db_ctgs.fasta -o 95otus_db_ctgs -c 0.95
```

## 1.4 Mapping

> Reads from the 214 bulk soil metagenomes were quality trimmed using Trimmomatic v0.3635 and then paired reads were mapped to 
> the viral contig database with Bowtie236, using default parameters. The output bam files were passed to BamM ‘filter’ v1.7.2 
> (http://ecogenomics.github.io/BamM/, accessed 15 December 2015) and reads that were aligned over ≥90% of their length at ≥95% > nucleic acid identity were retained.

1. Trimming: trimmomatic
- http://www.usadellab.org/cms/?page=trimmomatic
2. Mapping: BWA
- http://bio-bwa.sourceforge.net/
3. Bam filter: BamM 'filter
- http://ecogenomics.github.io/BamM/

## 1.5 Filtering

**References:**
- https://www.nature.com/articles/s41564-018-0190-y#ref-CR18
