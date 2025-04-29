import os
from collections import Counter, defaultdict
import json

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import tqdm
from pymutspec.annotation import transcriptor
from scipy.spatial.distance import pdist, squareform


from utils import (
    amino_acid_codes, prepare_aa_subst,
    calc_metrics, prepare_exp_aa_subst, prepare_rnd_exp_aa_subst,
    plot_obs_vs_exp,
)

color_mapping12 = {
    "C>A": "deepskyblue",
    "G>T": "deepskyblue",
    "C>G": "black",
    "G>C": "black",
    "C>T": "red",
    "G>A": "red",
    "T>A": "silver",
    "A>T": "silver",
    "T>C": "yellowgreen",
    "A>G": "yellowgreen",
    "T>G": "pink",
    "A>C": "pink",
}

debug = True
if debug:
    os.chdir('./PCAWG')


# with open('./data/human_genome_nuc_cnt.json') as f:
#     nuc_cnt_genomic = json.load(f)
with open('./data/human_cds_nuc_cnt.json') as f:
    nuc_cnt_cds = json.load(f)


def main():
    # read mutations
    obs = pd.read_csv('data/missense_mutations.csv').assign(count=1)
    obs['aa1'] = obs.aa_change.str[0]
    obs['aa2'] = obs.aa_change.str[-1]

    available_cancers = obs[obs.dataset == 'WES']['disease'].value_counts()

    # read amino acid freqs from protein files
    with open('./data/human_aa_cnt.json') as f:
        aa_freqs_dct_raw = json.load(f)
        aa_freqs_dct = {amino_acid_codes[k]: v for k,v in aa_freqs_dct_raw.items()}

    # read samples spectra
    spectra_wes192 = pd.read_csv(
        './data/Input_Data_PCAWG7_23K_Spectra_DB/Mutation_Catalogs_--_Spectra_of_Individual_Tumours/WES_Other.192.csv')
    spectra_wes192.loc[spectra_wes192.Strand == 'T', 'Mutation type'] = \
        spectra_wes192.loc[spectra_wes192.Strand == 'T', 'Mutation type'].str.translate(transcriptor)
    
    spectra_wes12 = spectra_wes192.groupby('Mutation type')\
        .sum().reset_index().rename(columns={'Mutation type': 'Mut'})
    wes_exp12 = spectra_wes12.Mut.apply(lambda x: nuc_cnt_cds[x[0]])
    wes_exp12.index = spectra_wes12.Mut
    wes_exp12 /= wes_exp12.sum() / 3
    spectra_wes12 = (spectra_wes12.set_index('Mut').T / wes_exp12.to_dict()).T
    spectra_wes12 = spectra_wes12 / spectra_wes12.sum()

    metrics_total = []
    i = 0
    nrows, ncols = 5, 5
    fig0, axs_sp = plt.subplots(nrows, ncols, figsize=(24, 15))
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*4-1.5))
    for disease in available_cancers[available_cancers>1000].index:
        cur_samples = [x for x in spectra_wes12.columns if x.startswith(disease)]
        if len(cur_samples) < 50:
            continue
        cur_spectra = spectra_wes12[cur_samples].T
        cur_spectra = cur_spectra[(cur_spectra > 0).sum(1) >= 4]
        cossim = 1 - squareform(pdist(cur_spectra, 'cosine'))
        cur_spectra = cur_spectra[(cossim < 0.7).sum(1) < len(cur_spectra) / 5]
        if len(cur_spectra) < 15:
            continue
        mean_cossim = np.mean(1-pdist(cur_spectra, 'cosine'))
        median_cossim = np.median(1-pdist(cur_spectra, 'cosine'))
        print(f"{disease[:15]}\t{len(cur_samples)}\t{len(cur_spectra)}\t{mean_cossim:.2f}\t{median_cossim:.2f}")
        
        cur_obs = obs[(obs['sample'].isin(cur_spectra.index.str.replace(disease+'::', ''))) & 
                      (obs.dataset == 'WES')]
        cur_obs = cur_obs[(cur_obs.aa1 != '*') & (cur_obs.aa2 != '*')]

        cur_spectrum = cur_spectra.mean().rename('MutSpec').reset_index()
        exp_aa_subst, _ = prepare_exp_aa_subst(cur_spectrum, 'MutSpec', 1)

        aa_subst = prepare_aa_subst(cur_obs, exp_aa_subst, aa_freqs_dct)
        cur_metrics = calc_metrics(aa_subst)
        cur_metrics['cancer'] = disease
        cur_metrics['model'] = 'neutral'
        metrics_total.append(cur_metrics)
        
        rnd_metrics = defaultdict(int)
        nreplics = 20
        for _ in range(nreplics):
            rnd_exp_aa_subst, _ = prepare_rnd_exp_aa_subst(gc=1)
            rnd_aa_subst = prepare_aa_subst(cur_obs, rnd_exp_aa_subst, aa_freqs_dct)
            cur_rnd_metrics = calc_metrics(rnd_aa_subst)
            for k,v in cur_rnd_metrics.items():
                rnd_metrics[k] += v / nreplics
        
        rnd_metrics = dict(rnd_metrics)
        rnd_metrics['cancer'] = disease
        rnd_metrics['model'] = 'random'
        metrics_total.append(rnd_metrics)

        print(disease, len(cur_obs), '##################')
        
        if i < 25:
            ax_sp = axs_sp[i // 5, i % 5]
            sns.barplot(cur_spectrum, x='Mut', y='MutSpec', palette=color_mapping12, ax=ax_sp)
            nmuts = cur_obs.shape[0]
            ax_sp.set_title(f'{disease} ({nmuts:.1f} mutations)')

            ax = axs[i // ncols, i % ncols]
            plot_obs_vs_exp(
                aa_subst, ax=ax, show=False, 
                text=f"{disease.split('-')[:15]}", text_x=-2.2, text_y=-4.)
        i += 1

    fig0.tight_layout()
    fig0.savefig('./figures/ms12cancers.pdf')

    fig.tight_layout()
    fig.savefig('./figures/neutral_model_fit_nuclear_somatic_cancers.pdf')
    plt.close()

    metrics_total_df = pd.DataFrame(metrics_total)\
        .set_index(['cancer', 'model'])
    metrics_total_df.to_csv('./data/nuclear_somatic_fit_metrics_cancers.csv', float_format='%g')
    print(metrics_total_df)


if __name__ == "__main__":
    main()