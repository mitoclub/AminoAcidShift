## read nuc and aa freqs
from collections import Counter
from Bio import SeqIO
import json

from utils import amino_acid_codes

bases = list('ACGT')
nuc_cnt_genomic = dict(zip(bases, [0,0,0,0]))

for rec in SeqIO.parse('./data/ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna', 'fasta'):
    cur_cnt = Counter(rec.seq)
    for nuc in bases:
        nuc_cnt_genomic[nuc] += cur_cnt[nuc]

with open('./data/human_genome_nuc_cnt.json', 'w') as handle:
    json.dump(nuc_cnt_genomic, handle)

print('Genome GC content:')
print((nuc_cnt_genomic['G'] + nuc_cnt_genomic['C']) / sum(nuc_cnt_genomic.values()))

bases = list('ACGT')
nuc_cnt = dict(zip(bases, [0,0,0,0]))

for rec in SeqIO.parse('./data/ncbi_dataset/data/GCF_000001405.40/cds_from_genomic.fna', 'fasta'):
    cur_cnt = Counter(rec.seq)
    for nuc in bases:
        nuc_cnt[nuc] += cur_cnt[nuc]

with open('./data/human_cds_nuc_cnt.json', 'w') as handle:
    json.dump(nuc_cnt, handle)

print('CDS GC content:')
print((nuc_cnt['G'] + nuc_cnt['C']) / sum(nuc_cnt.values()))

aa_cnt = dict(zip(amino_acid_codes, [0]*20))

for rec in SeqIO.parse('./data/ncbi_dataset/data/GCF_000001405.40/protein.faa', 'fasta'):
    cur_cnt = Counter(rec.seq)
    for aa in amino_acid_codes:
        if aa != '*':
            aa_cnt[aa] += cur_cnt[aa]

with open('./data/human_aa_cnt.json', 'w') as handle:
    json.dump(aa_cnt, handle)
