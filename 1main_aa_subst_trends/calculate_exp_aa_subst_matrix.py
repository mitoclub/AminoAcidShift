import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pymutspec.annotation import CodonAnnotation
from pymutspec.constants import possible_codons

alphabet = 'ACGT'


def collect_possible_changes(gc=2):
    coda = CodonAnnotation(gc)
    nucls = alphabet
    i = 1
    data = []
    for cdn1 in possible_codons:
        aa1 = coda.translate_codon(cdn1)
        for pic in range(3):
            nuc1 = cdn1[pic]
            for nuc2 in nucls:
                if nuc1 == nuc2:
                    continue
                cdn2 = list(cdn1)
                cdn2[pic] = nuc2
                cdn2 = ''.join(cdn2)
                aa2 = coda.translate_codon(cdn2)
                is_syn = aa1 == aa2
                sbs = f'{nuc1}>{nuc2}'
                data.append((pic, cdn1, cdn2, aa1, aa2, is_syn, sbs))
                i += 1

    df_changes = pd.DataFrame(data, columns=['pic', 'cdn1', 'cdn2', 'aa1', 'aa2', 'is_syn', 'sbs'])
    return df_changes


def nuc_spectrum_to_matrix(spec):
    '''
    convert dictionary of mutation counts to mutation matrix
    '''
    M = np.zeros((4,4))
    for i1,n1 in enumerate(alphabet):
        for i2,n2 in enumerate(alphabet):
            if n1!=n2:
                M[i2,i1] = spec[f"{n1}>{n2}"]
    # normalize off-diagonal rates (just for standardization, doesn't affect the results)
    M /= M.sum()
    # # will the diagonal with 'outflow' term to guarantee conservation of probability
    d = M.sum(axis=0)
    np.fill_diagonal(M,-d)

    return M


def cdn_spectrum_to_matrix(cdn_sbs):
    '''
    convert dictionary of mutation counts to mutation matrix
    '''
    n = len(possible_codons)
    M = np.zeros((n, n))
    for i1,cdn1 in enumerate(possible_codons):
        for i2,cdn2 in enumerate(possible_codons):
            if cdn1!=cdn2:
                val = cdn_sbs[(cdn1, cdn2)] if (cdn1, cdn2) in cdn_sbs.index else 0.
                M[i2,i1] = val
    # normalize off-diagonal rates (just for standardization, doesn't affect the results)
    M /= M.sum()
    # will the diagonal with 'outflow' term to guarantee conservation of probability
    d = M.sum(axis=0)
    np.fill_diagonal(M,-d)
    return M


def get_equilibrium_probabilities(M):
    evals, evecs = np.linalg.eig(M)
    # find zero eigenvalue
    ii = np.argmin(np.abs(evals))
    assert np.abs(evals[ii])<1e-10
    # pull out corresponding eigenvector, return normalized to sum_i p_i = 1
    p = evecs[:,ii]
    return p/p.sum()


amino_acid_codes = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "*",
}


def main():
    ms12 = pd.read_csv('../192/1data_derivation/dataset/MutSpecVertebrates12.csv.gz')
    ms12cytb = ms12[ms12.Gene == 'Cytb']
    mean_vert_ms12 = ms12cytb.groupby(['Mut']).MutSpec.mean().to_dict()

    df_changes = collect_possible_changes(2)
    df_changes['rate'] = df_changes['sbs'].map(mean_vert_ms12)
    df_changes['aa1'] = df_changes['aa1'].map(amino_acid_codes)
    df_changes['aa2'] = df_changes['aa2'].map(amino_acid_codes)
    aa_subst_vert = df_changes[(df_changes.aa1 != '*')&(df_changes.aa2 != '*')]\
        .groupby(['aa1', 'aa2'])['rate'].sum()
    aa_subst_vert.to_csv('eq_freqs/data/exp_aa_subst_vert.csv', float_format='%g')

    aa_subst_5cls = []
    classes_spectra = ms12cytb.groupby(['Class', 'Mut']).MutSpec.mean().unstack()
    for _cls in classes_spectra.index:
        mean_cls_ms12 = classes_spectra.loc[_cls].to_dict()

        df_changes = collect_possible_changes(2)
        df_changes['rate'] = df_changes['sbs'].map(mean_cls_ms12)
        df_changes['aa1'] = df_changes['aa1'].map(amino_acid_codes)
        df_changes['aa2'] = df_changes['aa2'].map(amino_acid_codes)
        aa_subst_vert = df_changes[(df_changes.aa1 != '*')&(df_changes.aa2 != '*')]\
            .groupby(['aa1', 'aa2'])['rate'].sum().reset_index()
        aa_subst_5cls.append(aa_subst_vert.assign(cls=_cls))
    
    aa_subst_5cls = pd.concat(aa_subst_5cls)
    aa_subst_5cls.to_csv('eq_freqs/data/exp_aa_subst_5cls.csv', index=False, float_format='%g')
    

if __name__ == "__main__":
    main()
