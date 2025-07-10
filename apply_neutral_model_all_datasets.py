import os
from collections import Counter

from Bio import SeqIO
import pandas as pd
import tqdm
from pymutspec.io import read_genbank_ref
from pymutspec.annotation import CodonAnnotation, transcriptor

from utils import (
    amino_acid_codes, prepare_aa_subst,
    calc_metrics, prepare_exp_aa_subst, 
)

##  inputs:
# - aa frequencies
# - spectra
# - mutations


coda = CodonAnnotation(2)
path_to_genbank = 'NC_012920.1.gb'

# ref_sites_df = read_genbank_ref(path_to_genbank)
# ref_sites_df = ref_sites_df[ref_sites_df.Codon.notna()]
# ref_sites_df['AA'] = ref_sites_df['Codon']\
#     .apply(coda.translate_codon).map(amino_acid_codes)

def generate_alt_codon(ref_codon, pos_in_codon: int, alt_nuc: str):
    assert pos_in_codon in [1, 2, 3]
    assert ref_codon[pos_in_codon-1] != alt_nuc
    ref_codon = list(ref_codon)
    ref_codon[pos_in_codon-1] = alt_nuc
    return ''.join(ref_codon)


def read_aa_counts_from_gb(path: str) -> dict:
    rec = next(SeqIO.parse(path, "genbank"))
    ref_df = pd.DataFrame([f.qualifiers for f in rec.features if f.type == "CDS"])
    # ref_df.drop(columns=["locus_tag", "ribosomal_slippage", "codon_start", 
    #                     "db_xref", 'gene_synonym'], inplace=True)
    ref_df["gene"] = ref_df["gene"].apply(lambda x: x[0])
    ref_df["product"] = ref_df["product"].apply(lambda x: x[0])
    ref_df["protein_id"] = ref_df["protein_id"].apply(lambda x: x[0])
    ref_df["translation"] = ref_df["translation"].apply(lambda x: x[0])

    genes_aa_counts_df = pd.DataFrame(ref_df.set_index('gene')['translation']\
                            .apply(Counter).to_dict()).T.fillna(0).astype(int)

    aa_counts_total = genes_aa_counts_df[genes_aa_counts_df.index != 'ND6']\
        .sum(0).rename('TOTALH').to_frame().T
    genes_aa_counts_df = pd.concat([genes_aa_counts_df, aa_counts_total], axis=0)\
        .rename(columns=amino_acid_codes)
    return genes_aa_counts_df.T.to_dict()


# def get_site_specific_aa_counts(sites):
#     aa_counts = ref_sites_df[ref_sites_df.Pos.isin(sites)]\
#         .query('AA != "*"').AA.value_counts().to_dict()
#     return aa_counts


def main():
    # read aa freqs for adjustment
    genes_aa_freqs = read_aa_counts_from_gb(path_to_genbank)

    # read mutations from all datasets
    chordates_species_mut = pd.read_csv('./vertebrates_aa_subst/dataset/obs_muts.csv')\
        .rename(columns={'RefAa': 'aa1', 'AltAa': 'aa2', 'ProbaMut': 'count'})\
            .query('gene != "ND6"')
    # print(chordates_species_mut.groupby('gene')['count'].sum())
    chordates_species_mut_Cytb = chordates_species_mut.query('gene == "Cytb"')
    chordates_species_mut_ND2 = chordates_species_mut.query('gene == "ND2"')

    megatree_mut = pd.read_csv('external_datasets/raw_human_megatree.csv')\
        .rename(columns={'Aa1': 'aa1', 'Aa2': 'aa2', 'ProbaFull': 'count'})
    megatree_mut = megatree_mut[(megatree_mut.TypeRef == 'CDS') & 
                                (megatree_mut.Label == 0) & 
                                (megatree_mut.GeneRef != 'ND6')]
    
    cancer_mut = pd.read_csv('https://raw.githubusercontent.com/mitoclub/mtdna-192component-mutspec-chordata/refs/heads/main/0cancer/data/mutations.csv')
    cancer_mut = cancer_mut[(cancer_mut.Type == 'CDS') & 
                            (cancer_mut.Label == 0) & 
                            (cancer_mut.GeneName != 'ND6')].assign(count=1)
    cancer_mut['aa1'] = cancer_mut.apply(lambda x: coda.translate_codon(x.Codon), axis=1)
    cancer_mut['aa2'] = cancer_mut.apply(lambda x: coda.translate_codon(x.AltCodon), axis=1)

    mitomap_mut = pd.read_csv('external_datasets/MutationsSomatic_MITOMAP_Foswiki.csv')
    mitomap_mut = mitomap_mut[mitomap_mut['Amino Acid Change'].str.match('[A-Z]-[A-Z]$')]
    mitomap_mut = mitomap_mut.query('Locus != "MT-ND6"').assign(count=1)
    mitomap_mut['aa1'] = mitomap_mut['Amino Acid Change'].str.split('-').str[0]
    mitomap_mut['aa2'] = mitomap_mut['Amino Acid Change'].str.split('-').str[1]

    mitomap_pathogenic = pd.read_csv('external_datasets/ConfirmedMutations_MITOMAP_Foswiki.csv')
    mitomap_pathogenic = mitomap_pathogenic[(mitomap_pathogenic['Locus Type'] == "Coding") & 
                                            (mitomap_pathogenic['Locus'] != "MT-ND6") &
                                            (mitomap_pathogenic['aaΔ-or-RNA'].str.match('[A-Z][0-9]{1,4}[A-Z]$'))]
    mitomap_pathogenic['aa1'] = mitomap_pathogenic['aaΔ-or-RNA'].str[0]
    mitomap_pathogenic['aa2'] = mitomap_pathogenic['aaΔ-or-RNA'].str[-1]
    mitomap_pathogenic['count'] = 1

    gtex_mut = pd.read_csv('./external_datasets/gtex_annotated.csv')
    gtex_mut = gtex_mut[(gtex_mut.Type == 'CDS') & (gtex_mut.PosInCodon > 0)].assign(count=1)
    gtex_mut['AltCodon'] = gtex_mut.apply(lambda x: generate_alt_codon(x.Codon, x.PosInCodon, x.DerivedAllele), axis=1)
    gtex_mut['aa1'] = gtex_mut.apply(lambda x: coda.translate_codon(x.Codon), axis=1)
    gtex_mut['aa2'] = gtex_mut.apply(lambda x: coda.translate_codon(x.AltCodon), axis=1)

    global_chordates_mut = pd.read_csv(
        './external_datasets/raw_global_vert_cytb_mutations.tsv', sep='\t'
    ).query('Label == 0').rename(columns={'RefAa': 'aa1', 'AltAa': 'aa2', 'ProbaFull': 'count'})

    # read spectra
    chordates_species_ms12_raw = pd.read_csv('../192/1data_derivation/dataset/MutSpecVertebrates12.csv.gz')
    chordates_species_ms12_Cytb = chordates_species_ms12_raw.query('Gene == "Cytb"')\
        .groupby(['Mut']).MutSpec.mean().reset_index()
    chordates_species_ms12_ND2 = chordates_species_ms12_raw.query('Gene == "ND2"')\
        .groupby(['Mut']).MutSpec.mean().reset_index()
    
    # # MBE spectrum for all datasets
    # global_chordates_ms12 = megatree_ms12 = cancer_ms12 = chordates_species_ms12
    
    # custom spectra
    megatree_ms12 = pd.read_csv('../human_mtDNA_megatree/data/spectra/ms12syn.csv')
    megatree_ms12['Mut'] = megatree_ms12['Mut'].str.translate(transcriptor)
    cancer_ms12 = pd.read_csv('../mtdnaMutSpecOfCancers/data/mutspecs/cancer_mutspec12custom.csv')\
        .rename(columns={'MutSpec_AllWithoutDloop': 'MutSpec'})
    cancer_ms12['Mut'] = cancer_ms12['Mut'].str.translate(transcriptor)

    # global_chordates_ms12 = pd.read_csv("external_datasets/global_cytb_chordates_ms12syn.tsv", sep='\t')

    # use mean chordates ms for global dataset
    total_spectra = [chordates_species_ms12_Cytb, chordates_species_ms12_Cytb, chordates_species_ms12_ND2, 
                     megatree_ms12, cancer_ms12, cancer_ms12, cancer_ms12, cancer_ms12]
    total_obs = [global_chordates_mut, chordates_species_mut_Cytb, chordates_species_mut_ND2, 
                 megatree_mut, cancer_mut, gtex_mut, mitomap_mut, mitomap_pathogenic]
    total_freqs_genes = ['CYTB', 'CYTB', 'ND2', 'TOTALH', 'TOTALH', 'TOTALH', 'TOTALH', 'TOTALH']
    total_labels = ['global_chordates_cytb', 'chordates_species_Cytb', 'chordates_species_ND2', 
                    'megatree', 'cancer', 'gtex', 'mitomap', 'mitomap_pathogenic']
    
    metrics_total = []
    for ms12, obs, freqs_gene, dataset in zip(
        total_spectra, total_obs, total_freqs_genes, total_labels):

        print(dataset, len(obs), '##################')

        exp_aa_subst, _exp_aa_subst_matrix = prepare_exp_aa_subst(ms12, 'MutSpec', 2)

        # Select clade OBS AA substitutions
        obs = obs[(obs.aa1 != '*') & (obs.aa2 != '*')]

        aa_subst = prepare_aa_subst(obs, exp_aa_subst, genes_aa_freqs[freqs_gene])
        cur_metrics = calc_metrics(aa_subst)
        cur_metrics['dataset'] = dataset
        metrics_total.append(cur_metrics)

    metrics_total_df = pd.DataFrame(metrics_total).set_index(['dataset'])
    metrics_total_df.to_csv('datasets_fit_metrics.csv', float_format='%g')
    print(metrics_total_df['r2,spearman_corr,mut_count'.split(',')])


if __name__ == "__main__":
    main()