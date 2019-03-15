from Bio import Entrez, SeqIO
import pandas as pd

Entrez.email = 'my.email@mycompany.XXX'

data_inputs = []
with open('assembly_virus_complete_genomes.txt', 'r') as f:
    for line in f.readlines():
        data_inputs.append(line.replace('\n',''))

with open('28350_viral_genomes.fasta', 'w') as f:
    for genome_id in data_inputs:
        record = Entrez.efetch(db='nucleotide', id=genome_id, rettype='fasta', retmode='text')
        filename = 'generated/genBankRecord_{}.gb'.format(genome_id)
        f.write(record.read())
