import os
from collections import Counter

from Bio import SeqIO
import pandas as pd
import tqdm
from pymutspec.io import read_genbank_ref
from pymutspec.annotation import CodonAnnotation

from utils import (
    amino_acid_codes, prepare_aa_subst,calc_metrics, 
    prepare_exp_aa_subst, prepare_rnd_exp_aa_subst,
)

debug = False
if debug:
    os.chdir('./sars-cov-2_aa_subst/')

coda = CodonAnnotation(1)

ref_sites_df = read_genbank_ref('data/NC_045512.2.gb')
ref_sites_df = ref_sites_df[ref_sites_df.Codon.notna()]
ref_sites_df['AA'] = ref_sites_df['Codon']\
    .apply(coda.translate_codon).map(amino_acid_codes)


def read_aa_counts_from_gb(path: str) -> dict:
    rec = next(SeqIO.parse(path, "genbank"))
    ref_df = pd.DataFrame([f.qualifiers for f in rec.features if f.type == "CDS"])
    ref_df.drop(columns=["locus_tag", "ribosomal_slippage", "codon_start", 
                        "db_xref", 'gene_synonym'], inplace=True)
    ref_df["gene"] = ref_df["gene"].apply(lambda x: x[0])
    ref_df["product"] = ref_df["product"].apply(lambda x: x[0])
    ref_df["protein_id"] = ref_df["protein_id"].apply(lambda x: x[0])
    ref_df["translation"] = ref_df["translation"].apply(lambda x: x[0])
    ref_df = ref_df[ref_df['product'] != 'ORF1a polyprotein']

    aa_counts_df = pd.DataFrame(ref_df.set_index('gene')['translation']\
                            .apply(Counter).to_dict()).T.fillna(0).astype(int)

    aa_counts_total_dct = aa_counts_df.rename(columns=amino_acid_codes).sum(0).to_dict()
    return aa_counts_total_dct


def get_site_specific_aa_counts(sites):
    aa_counts = ref_sites_df[ref_sites_df.Pos.isin(sites)]\
        .query('AA != "*"').AA.value_counts().to_dict()
    return aa_counts


def main():
    aa_freqs_total_dct = read_aa_counts_from_gb("data/NC_045512.2.gb")

    clades_spectra = pd.read_csv('data/rates_by_clade.csv')
    clades_spectra['Mut'] = clades_spectra['mut_type'].str.replace('to', '>')

    # https://media.githubusercontent.com/media/jbloomlab/SARS2-mut-spectrum/refs/heads/main/results/mutation_counts/aggregated.csv
    obs = pd.read_csv('data/aggregated.csv')
    obs = obs[(obs['subset'] == 'all') & (obs['synonymous'] == False) & (obs['exclude'] == False)]

    def _same_aa_mut(aa_mutation: str):
        """filter syn mutations? and maybe common sites for several CDS"""
        variants = aa_mutation.split(';')
        return variants.count(variants[0]) == len(variants)
    
    obs = obs[obs.aa_mutation.apply(_same_aa_mut)]
    obs['aa1'] = obs['aa_mutation'].str[0]
    obs['aa2'] = obs['aa_mutation'].str[-1]

    metrics_total = []
    for clade in tqdm.tqdm(clades_spectra.clade.unique(), desc='Clades'):
        # Select clade OBS substitutions
        obs_clade = obs[obs['clade'] == clade].copy()
        obs_clade = obs_clade[(obs_clade.aa1 != '*') & (obs_clade.aa2 != '*')]

        # neutral model
        spectrum_clade = clades_spectra[clades_spectra['clade'] == clade]
        exp_aa_subst, _ = prepare_exp_aa_subst(spectrum_clade, 'rate', 1)
        aa_subst = prepare_aa_subst(obs_clade, exp_aa_subst, aa_freqs_total_dct)
        cur_metrics = calc_metrics(aa_subst)
        cur_metrics['model'] = 'neutral'
        cur_metrics['clade'] = clade
        cur_metrics['replica'] = 1
        metrics_total.append(cur_metrics)

        # random model
        for i in range(1, 21):
            exp_aa_subst_rnd, _ = prepare_rnd_exp_aa_subst(1)
            aa_subst = prepare_aa_subst(obs_clade, exp_aa_subst_rnd, aa_freqs_total_dct)
            cur_metrics = calc_metrics(aa_subst)
            cur_metrics['model'] = 'random'
            cur_metrics['clade'] = clade
            cur_metrics['replica'] = i
            metrics_total.append(cur_metrics)

    metrics_total_df = pd.DataFrame(metrics_total)\
        .set_index(['model', 'clade', 'replica'])
    metrics_total_df.to_csv('data/rnd_fit_metrics.csv', float_format='%g')


if __name__ == "__main__":
    main()