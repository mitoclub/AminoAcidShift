import os
from collections import Counter, defaultdict
import json

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import tqdm
from pymutspec.annotation import transcriptor


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


with open('./data/human_genome_nuc_cnt.json') as f:
    nuc_cnt_genomic = json.load(f)
with open('./data/human_cds_nuc_cnt.json') as f:
    nuc_cnt_cds = json.load(f)


def main():
    # read mutations
    obs = pd.read_csv('data/missense_mutations.csv').assign(count=1)
    obs['aa1'] = obs.aa_change.str[0]
    obs['aa2'] = obs.aa_change.str[-1]

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
    spectra_wes12 = spectra_wes12.reset_index()

    spectra_wgs192 = pd.read_csv(
        './data/Input_Data_PCAWG7_23K_Spectra_DB/Mutation_Catalogs_--_Spectra_of_Individual_Tumours/WGS_Other.192.csv')
    spectra_wgs192.loc[spectra_wgs192.Strand == 'T', 'Mutation type'] = \
        spectra_wgs192.loc[spectra_wgs192.Strand == 'T', 'Mutation type'].str.translate(transcriptor)
    
    spectra_wgs12 = spectra_wgs192.groupby('Mutation type')\
        .sum().reset_index().rename(columns={'Mutation type': 'Mut'})
    wgs_exp12 = spectra_wgs12.Mut.apply(lambda x: nuc_cnt_genomic[x[0]])
    wgs_exp12.index = spectra_wgs12.Mut
    wgs_exp12 /= wgs_exp12.sum() / 3
    spectra_wgs12 = (spectra_wgs12.set_index('Mut').T / wes_exp12.to_dict()).T
    spectra_wgs12 = spectra_wgs12 / spectra_wgs12.sum()
    spectra_wgs12 = spectra_wgs12.reset_index()

    metrics_total = []
    i = 0
    nrows, ncols = 4, 5
    fig0, axs_sp = plt.subplots(nrows, ncols, figsize=(24, 15))
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*4-1.5))
    for sample in obs.groupby('sample').chr.count().sort_values(ascending=False).head(150).index:
        cur_obs = obs[obs['sample'] == sample]
        cur_obs = cur_obs[(cur_obs.aa1 != '*') & (cur_obs.aa2 != '*')]

        disease, dataset = cur_obs[['disease','dataset']].iloc[0]

        try:
            if dataset == 'WES':
                cur_spectrum = spectra_wes12[['Mut', f'{disease}::{sample}']]
            elif dataset == 'WGS':
                cur_spectrum = spectra_wgs12[['Mut', f'{disease}::{sample}']]
            else:
                raise
        except Exception as e:
            print(e)
            continue

        cur_spectrum.columns = ['Mut', 'MutSpec']
        exp_aa_subst, _ = prepare_exp_aa_subst(cur_spectrum, 'MutSpec', 1)

        aa_subst = prepare_aa_subst(cur_obs, exp_aa_subst, aa_freqs_dct)
        cur_metrics = calc_metrics(aa_subst)
        cur_metrics['disease'] = disease
        cur_metrics['sample'] = sample
        cur_metrics['model'] = 'neutral'
        metrics_total.append(cur_metrics)
        
        rnd_metrics = defaultdict(int)
        nreplics = 1
        for _ in range (nreplics):
            rnd_exp_aa_subst, _ = prepare_rnd_exp_aa_subst(gc=1)
            rnd_aa_subst = prepare_aa_subst(cur_obs, rnd_exp_aa_subst, aa_freqs_dct)
            cur_rnd_metrics = calc_metrics(rnd_aa_subst)
            for k,v in cur_rnd_metrics.items():
                rnd_metrics[k] += v / nreplics
        
        rnd_metrics = dict(rnd_metrics)
        rnd_metrics['disease'] = disease
        rnd_metrics['sample'] = sample
        rnd_metrics['model'] = 'random'
        metrics_total.append(rnd_metrics)

        print(disease, sample, len(cur_obs), '##################')
        
        if i < 20:
            ax_sp = axs_sp[i // 5, i % 5]
            sns.barplot(cur_spectrum, x='Mut', y='MutSpec', palette=color_mapping12, ax=ax_sp)
            nmuts = cur_obs.shape[0]
            ax_sp.set_title(f'{disease}::{sample}::{dataset}\n({nmuts:.1f} mutations)')

            ax = axs[i // ncols, i % ncols]
            plot_obs_vs_exp(
                aa_subst, ax=ax, show=False, 
                text=f"{disease.split('-')[0]}::{sample}", text_x=-2.2, text_y=-4.)
        i += 1
    
    # # unshow empty plots
    # axs[i // ncols, i % ncols].axis('off')
    # for x in range(2):
    #     i += 1
    #     axs[i // ncols, i % ncols].axis('off')
        

    fig0.tight_layout()
    fig0.savefig('./figures/ms12.pdf')

    fig.tight_layout()
    fig.savefig('./figures/neutral_model_fit_nuclear_somatic.pdf')
    plt.close()

    metrics_total_df = pd.DataFrame(metrics_total)\
        .set_index(['disease', 'sample', 'model'])
    metrics_total_df.to_csv('./data/nuclear_somatic_fit_metrics.csv', float_format='%g')
    print(metrics_total_df)


if __name__ == "__main__":
    main()