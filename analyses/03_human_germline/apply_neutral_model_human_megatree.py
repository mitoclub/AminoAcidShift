from collections import Counter
from pathlib import Path

from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd


from aminoacidshift import (
    amino_acid_codes, prepare_aa_subst,
    calc_metrics, prepare_exp_aa_subst,
    plot_obs_vs_exp,
)
from aminoacidshift.paths import EXTERNAL, REFERENCE

HERE = Path(__file__).resolve().parent


def read_aa_counts_from_gb(path=None):
    if path is None:
        path = REFERENCE / "NC_012920.1.gb"
    refseq = SeqIO.parse(path, 'genbank')
    refseq = list(refseq)[0]

    aa_freqs = dict()
    for fea in refseq.features:
        if fea.type == 'CDS':
            # print(fea.qualifiers['gene'])
            # print(fea.qualifiers['translation'])
            # print('---')
            aa_freqs[fea.qualifiers['gene'][0]] = dict(Counter(fea.qualifiers['translation'][0]))
    df = pd.DataFrame(aa_freqs).T.rename(columns=amino_acid_codes)
    return df


def main():
    spectra = pd.read_csv(EXTERNAL / "human_megatree_genes_spectra.csv", index_col=0)

    # read amino acid freqs from protein files
    aa_freqs = read_aa_counts_from_gb()

    obs = pd.read_csv(EXTERNAL / "raw_human_megatree.tsv", sep="\t")
    obs = obs[(obs.TypeRef == 'CDS') & (obs.Label == 0)]\
        .rename(columns={'Aa1': 'aa1', 'Aa2': 'aa2', 'ProbaFull': 'count'})

    metrics_total = []
    i = 0
    nrows, ncols = 4, 4
    fig, axs = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*4-1.5))
    for gene in spectra.index.unique():
        cur_spectrum = spectra.loc[gene]
        exp_aa_subst, _ = prepare_exp_aa_subst(cur_spectrum, 'MutSpec', 2)
        
        cur_obs = obs[obs['GeneRef'] == gene].copy()
        print(gene, len(cur_obs), '##################')

        # for total sites set
        cur_aa_freqs_dct = aa_freqs.loc[gene].to_dict()
        aa_subst = prepare_aa_subst(cur_obs, exp_aa_subst, cur_aa_freqs_dct)
        cur_metrics = calc_metrics(aa_subst)
        cur_metrics['gene'] = gene
        metrics_total.append(cur_metrics)

        ax = axs[i // ncols, i % ncols]
        plot_obs_vs_exp(
            aa_subst, ax=ax, show=False, 
            text=f"{gene}", text_x=-2.2, text_y=-4.)
        i += 1
    
    # unshow empty plots
    axs[i // ncols, i % ncols].axis('off')
    for x in range(2):
        i += 1
        axs[i // ncols, i % ncols].axis('off')
        

    plt.tight_layout()
    plt.savefig(HERE / "neutral_model_fit_megatree.pdf", dpi=300)
    plt.close()

    metrics_total_df = pd.DataFrame(metrics_total)\
        .set_index('gene')
    metrics_total_df.to_csv(HERE / "megatree_fit_metrics.csv", float_format="%g")
    print(metrics_total_df)


if __name__ == "__main__":
    main()