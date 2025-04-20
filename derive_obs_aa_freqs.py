
# TODO compete this script

from collections import Counter
from Bio import SeqIO
import pandas as pd

from utils import amino_acid_codes

gene = 'Cytb'
# for gene_fn in glob.glob('../../192/1data_derivation/nemu_input/*.fasta'):
gene_fn = f'../../192/1data_derivation/nemu_input/{gene}.fasta'
gene_seqs = list(SeqIO.parse(gene_fn, 'fasta'))
# sizes = [len(x) for x in gene_seqs]
gene_seqs = [x for x in gene_seqs if len(x) > 370 and len(x) < 385]

obs_aa_freqs = []
for rec in gene_seqs:
    cur_cls = rec.description.split('###')[1].split(';')[1].split('_')[0]
    row = dict(Counter(rec.seq))
    row['cls'] = cur_cls
    row['len'] = len(rec)
    obs_aa_freqs.append(row)

obs_aa_freqs = pd.DataFrame(obs_aa_freqs).fillna(0)
del obs_aa_freqs['B']
del obs_aa_freqs['X']
obs_aa_freqs = obs_aa_freqs[obs_aa_freqs.cls.isin(obs_aa_freqs.cls.value_counts().head().index)].set_index('cls')
obs_aa_freqs = obs_aa_freqs.div(obs_aa_freqs['len'], axis=0)
del obs_aa_freqs['len']
obs_aa_freqs = obs_aa_freqs.melt(ignore_index=False, var_name='aa', value_name='freq').reset_index()
obs_aa_freqs['aa'] = obs_aa_freqs.aa.map(amino_acid_codes)
obs_aa_freqs.to_csv('data/obs_aa_freqs.csv', index=False)