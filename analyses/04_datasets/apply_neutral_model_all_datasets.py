from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from pymutspec.annotation import CodonAnnotation, transcriptor

from aminoacidshift import (
    amino_acid_codes, prepare_aa_subst,
    calc_metrics, prepare_exp_aa_subst,
    get_random_spectra,
)
from aminoacidshift.paths import CANCER_MUTSPEC, EXTERNAL, MEGATREE, NEMU_192, REFERENCE, VERTEBRATES_OBS

coda = CodonAnnotation(2)
path_to_genbank = REFERENCE / "NC_012920.1.gb"

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

    ## Mutations from comparative-species level Cytb and ND2

    # global_chordates_mut = pd.read_csv(
    #     './external_datasets/raw_global_vert_cytb_mutations.tsv', sep='\t'
    # ).query('Label == 0').rename(columns={'RefAa': 'aa1', 'AltAa': 'aa2', 'ProbaFull': 'count'})

    mammals_dir = EXTERNAL / "mammals_mutations"
    mammals_cytb_mut_gtr = pd.read_csv(mammals_dir / "cytb_gtr.tsv", sep="\t")
    mammals_cytb_mut_unr = pd.read_csv(mammals_dir / "cytb_unrest.tsv", sep="\t")
    mammals_cytb_mut_gtr = mammals_cytb_mut_gtr.query('Label == 0 & ProbaFull > 0.3')\
        .rename(columns={'RefAa': 'aa1', 'AltAa': 'aa2', 'ProbaFull': 'count'})
    mammals_cytb_mut_unr = mammals_cytb_mut_unr.query('Label == 0 & ProbaFull > 0.3')\
        .rename(columns={'RefAa': 'aa1', 'AltAa': 'aa2', 'ProbaFull': 'count'})

    mammals_nd1_mut_gtr = pd.read_csv(mammals_dir / "nd1_gtr.tsv", sep="\t")
    mammals_nd1_mut_unr = pd.read_csv(mammals_dir / "nd1_unrest.tsv", sep="\t")
    mammals_nd1_mut_gtr = mammals_nd1_mut_gtr.query('Label == 0 & ProbaFull > 0.3')\
        .rename(columns={'RefAa': 'aa1', 'AltAa': 'aa2', 'ProbaFull': 'count'})
    mammals_nd1_mut_unr = mammals_nd1_mut_unr.query('Label == 0 & ProbaFull > 0.3')\
        .rename(columns={'RefAa': 'aa1', 'AltAa': 'aa2', 'ProbaFull': 'count'})


    ## VERTEBRATES Cytb and ND2 germline mutations
    chordates_species_mut = pd.read_csv(VERTEBRATES_OBS)\
        .rename(columns={'RefAa': 'aa1', 'AltAa': 'aa2', 'ProbaMut': 'count'})
        
    chordates_species_mut_Cytb = chordates_species_mut.query('gene == "Cytb"')
    chordates_species_mut_ND2 = chordates_species_mut.query('gene == "ND2"')
    del chordates_species_mut

    ## Human mtDNA germline mutations from megatree
    megatree_mut = pd.read_csv(EXTERNAL / "raw_human_megatree.tsv", sep="\t")\
        .rename(columns={'Aa1': 'aa1', 'Aa2': 'aa2', 'ProbaFull': 'count'})
    megatree_mut = megatree_mut[(megatree_mut.TypeRef == 'CDS') & 
                                (megatree_mut.Label == 0) & 
                                (megatree_mut.GeneRef != 'ND6')]
    
    ## Human mtDNA somatic mutations
    cancer_mut = pd.read_csv('https://raw.githubusercontent.com/mitoclub/mtdna-192component-mutspec-chordata/refs/heads/main/0cancer/data/mutations.csv')
    cancer_mut = cancer_mut[(cancer_mut.Type == 'CDS') & 
                            (cancer_mut.Label == 0) & 
                            (cancer_mut.GeneName != 'ND6')].assign(count=1)
    cancer_mut['aa1'] = cancer_mut.apply(lambda x: coda.translate_codon(x.Codon), axis=1)
    cancer_mut['aa2'] = cancer_mut.apply(lambda x: coda.translate_codon(x.AltCodon), axis=1)

    ## Human mtDNA somatic and pathogenic mutations from MITOMAP
    mitomap_mut = pd.read_csv(EXTERNAL / "MutationsSomatic_MITOMAP_Foswiki.csv")
    mitomap_mut = mitomap_mut[mitomap_mut['Amino Acid Change'].str.match('[A-Z]-[A-Z]$')]
    mitomap_mut = mitomap_mut.query('Locus != "MT-ND6"').assign(count=1)
    mitomap_mut['aa1'] = mitomap_mut['Amino Acid Change'].str.split('-').str[0]
    mitomap_mut['aa2'] = mitomap_mut['Amino Acid Change'].str.split('-').str[1]

    mitomap_pathogenic = pd.read_csv(EXTERNAL / "ConfirmedMutations_MITOMAP_Foswiki.csv")
    mitomap_pathogenic = mitomap_pathogenic[(mitomap_pathogenic['Locus Type'] == "Coding") & 
                                            (mitomap_pathogenic['Locus'] != "MT-ND6") &
                                            (mitomap_pathogenic['aaΔ-or-RNA'].str.match('[A-Z][0-9]{1,4}[A-Z]$'))]
    mitomap_pathogenic['aa1'] = mitomap_pathogenic['aaΔ-or-RNA'].str[0]
    mitomap_pathogenic['aa2'] = mitomap_pathogenic['aaΔ-or-RNA'].str[-1]
    mitomap_pathogenic['count'] = 1

    ## Human mtDNA nearly somatic mutations from gTEX
    gtex_mut = pd.read_csv(EXTERNAL / "gtex_annotated.csv")
    gtex_mut = gtex_mut[(gtex_mut.Type == 'CDS') & (gtex_mut.PosInCodon > 0)].assign(count=1)
    gtex_mut['AltCodon'] = gtex_mut.apply(lambda x: generate_alt_codon(x.Codon, x.PosInCodon, x.DerivedAllele), axis=1)
    gtex_mut['aa1'] = gtex_mut.apply(lambda x: coda.translate_codon(x.Codon), axis=1)
    gtex_mut['aa2'] = gtex_mut.apply(lambda x: coda.translate_codon(x.AltCodon), axis=1)

    # read spectra
    chordates_species_ms12_raw = pd.read_csv(
        NEMU_192 / "1data_derivation/dataset/MutSpecVertebrates12.csv.gz")

    def _class_spectrum(gene, cls="Mammalia"):
        return (
            chordates_species_ms12_raw.query("Gene == @gene & Class == @cls")
            .groupby(["Mut"]).MutSpec.mean()
            .reset_index()
        )

    chordates_species_ms12_Cytb = _class_spectrum("Cytb")
    chordates_species_ms12_ND2 = _class_spectrum("ND2")
    chordates_species_ms12_ND1 = _class_spectrum("ND1")
    assert len(chordates_species_ms12_Cytb) == 12
    assert len(chordates_species_ms12_ND2) == 12
    if len(chordates_species_ms12_ND1) != 12:
        print(
            "Skipping mam_nd1_gtr and mam_nd1_unr: MutSpecVertebrates12 has no "
            f"12-class ND1 Mammalia spectrum ({len(chordates_species_ms12_ND1)} mutation types)."
        )
        chordates_species_ms12_ND1 = None

    megatree_ms12 = pd.read_csv(MEGATREE / "data/spectra/ms12syn.csv")
    megatree_ms12["Mut"] = megatree_ms12["Mut"].str.translate(transcriptor)
    cancer_ms12 = pd.read_csv(
        CANCER_MUTSPEC / "data/mutspecs/cancer_mutspec12custom.csv"
    ).rename(columns={"MutSpec_AllWithoutDloop": "MutSpec"})
    cancer_ms12["Mut"] = cancer_ms12["Mut"].str.translate(transcriptor)

    # global_chordates_ms12 = pd.read_csv("external_datasets/global_cytb_chordates_ms12syn.tsv", sep='\t')

    # use mean chordates ms for global dataset
    total_spectra = [
        chordates_species_ms12_Cytb, chordates_species_ms12_Cytb, 
        chordates_species_ms12_ND1, chordates_species_ms12_ND1, 
        chordates_species_ms12_Cytb, chordates_species_ms12_ND2,
        megatree_ms12, cancer_ms12, cancer_ms12, cancer_ms12, cancer_ms12]
    total_obs = [
        mammals_cytb_mut_gtr, mammals_cytb_mut_unr, 
        mammals_nd1_mut_gtr, mammals_nd1_mut_unr, 
        chordates_species_mut_Cytb, chordates_species_mut_ND2, 
        megatree_mut, cancer_mut, gtex_mut, mitomap_mut, mitomap_pathogenic]
    total_freqs_genes = ['CYTB', 'CYTB', 'ND1', 'ND1', 'CYTB', 'ND2', 
                         'TOTALH', 'TOTALH', 'TOTALH', 'TOTALH', 'TOTALH']
    total_labels = [
        'mam_cytb_gtr', 'mam_cytb_unr',
        'mam_nd1_gtr', 'mam_nd1_unr',
        'chordates_species_Cytb', 'chordates_species_ND2', 
        'megatree', 'cancer', 'gtex', 'mitomap', 'mitomap_pathogenic']

    if chordates_species_ms12_ND1 is None:
        keep = [label not in {"mam_nd1_gtr", "mam_nd1_unr"} for label in total_labels]
        total_spectra = [item for item, ok in zip(total_spectra, keep) if ok]
        total_obs = [item for item, ok in zip(total_obs, keep) if ok]
        total_freqs_genes = [item for item, ok in zip(total_freqs_genes, keep) if ok]
        total_labels = [item for item, ok in zip(total_labels, keep) if ok]

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

        # get metrics on random spectra
        r2_rnd = []
        pearson_corr_rnd = []
        rnd_spectra = get_random_spectra(100)
        for _, ms_rnd in rnd_spectra.groupby('replica'):
            exp_aa_subst, _1 = prepare_exp_aa_subst(ms_rnd, 'MutSpec', 2)
            aa_subst_rnd = prepare_aa_subst(obs, exp_aa_subst, genes_aa_freqs[freqs_gene])
            cur_metrics_rnd = calc_metrics(aa_subst_rnd)
            r2_rnd.append(cur_metrics_rnd['r2'])
            pearson_corr_rnd.append(cur_metrics_rnd['pearson_corr_clr'])
        cur_metrics['r2_rnd'] = np.mean(r2_rnd)
        cur_metrics['pearson_corr_rnd'] = np.mean(pearson_corr_rnd)

        metrics_total.append(cur_metrics)

    out_path = Path(__file__).resolve().parent / "datasets_fit_metrics.csv"
    metrics_total_df = pd.DataFrame(metrics_total).set_index(['dataset'])
    metrics_total_df['pearson_corr_clr_squared'] = metrics_total_df['pearson_corr_clr'] ** 2
    metrics_total_df['pearson_corr_clr_squared_rnd'] = metrics_total_df['pearson_corr_rnd'] ** 2
    metrics_total_df.to_csv(out_path, float_format='%g')

    print(metrics_total_df['pearson_corr_clr_squared,pearson_corr_clr_squared_rnd,mut_count'.split(',')])


if __name__ == "__main__":
    main()