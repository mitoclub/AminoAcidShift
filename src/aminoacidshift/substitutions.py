import numpy as np
import pandas as pd
from pymutspec.annotation import CodonAnnotation
from pymutspec.constants import possible_codons
from scipy.stats import uniform

from aminoacidshift.constants import (
    alphabet,
    amino_acid_codes,
    grantham_long,
    transitions,
    transversions,
)


def collect_possible_changes(gc, spectrum_dict=None):
    coda = CodonAnnotation(gc)
    nucls = alphabet
    data = []
    for cdn1 in possible_codons:
        aa1 = amino_acid_codes[coda.translate_codon(cdn1)]
        for pic in range(3):
            nuc1 = cdn1[pic]
            for nuc2 in nucls:
                if nuc1 == nuc2:
                    continue
                cdn2 = list(cdn1)
                cdn2[pic] = nuc2
                cdn2 = "".join(cdn2)
                aa2 = amino_acid_codes[coda.translate_codon(cdn2)]
                is_syn = aa1 == aa2
                sbs = f"{nuc1}>{nuc2}"
                data.append((pic, cdn1, cdn2, aa1, aa2, is_syn, sbs))

    df_changes = pd.DataFrame(
        data, columns=["pic", "cdn1", "cdn2", "aa1", "aa2", "is_syn", "sbs"]
    )

    if spectrum_dict is not None:
        df_changes["rate"] = df_changes["sbs"].map(spectrum_dict)

    return df_changes


def nuc_spectrum_to_matrix(spec):
    """Convert a dictionary of mutation counts to a mutation matrix."""
    matrix = np.zeros((4, 4))
    for i1, n1 in enumerate(alphabet):
        for i2, n2 in enumerate(alphabet):
            if n1 != n2:
                matrix[i1, i2] = spec[f"{n1}>{n2}"]
    matrix /= matrix.sum()
    return matrix


def cdn_spectrum_to_matrix(df_changes: pd.DataFrame):
    cdn_sbs = df_changes.groupby(["cdn1", "cdn2"])["rate"].sum().reset_index()
    pcdn = cdn_sbs.pivot(index="cdn1", columns="cdn2", values="rate").fillna(0.0)
    pcdn = (pcdn.T / pcdn.sum(axis=1)).T
    return pcdn


def get_equilibrium_probabilities(pmatrix: np.ndarray):
    """
    Equilibrium probabilities from the eigenvector of eigenvalue ~ 0.

    Solves pi @ Q = 0, equivalently Q.T @ pi.T = 0.
    """
    generator = pmatrix.copy()
    generator = generator - np.diag(np.sum(generator, axis=1))

    values, vectors = np.linalg.eig(generator.T)
    pi = vectors[:, np.isclose(values, 0)].real.flatten()
    pi = pi / pi.sum()
    return pi


def get_equilibrium_freqs(spectrum: pd.DataFrame, rate_col="MutSpec", gc=1):
    coda = CodonAnnotation(gc)
    spectrum_dict = spectrum.set_index("Mut")[rate_col].to_dict()
    df_changes = collect_possible_changes(gc, spectrum_dict)

    matrix = cdn_spectrum_to_matrix(df_changes)
    eq_prob = get_equilibrium_probabilities(matrix.values)

    eq_freqs_cdn = pd.Series(dict(zip(possible_codons, eq_prob)))
    eq_freqs_cdn.name = "eq_freq"
    eq_freqs_cdn.index.name = "cdn"
    eq_freqs_cdn = eq_freqs_cdn.reset_index()
    eq_freqs_cdn["aa"] = eq_freqs_cdn["cdn"].map(coda.translate_codon).map(amino_acid_codes)

    eq_freqs_aa = eq_freqs_cdn[eq_freqs_cdn.aa != "*"].groupby("aa")["eq_freq"].sum()
    eq_freqs_aa /= eq_freqs_aa.sum()
    eq_freqs_aa = eq_freqs_aa.sort_values(ascending=False).reset_index()

    return eq_freqs_cdn, eq_freqs_aa


def get_random_spectra(n=10, tstv_sampling="poisson"):
    """Random 12-component spectra.

    For a typical spectrum the transition/transversion ratio is between 2 and 10.
    ``tstv_sampling`` is ``poisson`` or ``random``.
    """
    columns = [f"{n1}>{n2}" for n1 in alphabet for n2 in alphabet if n1 != n2]
    data = uniform.rvs(size=(n, 12))
    rnd_spectra = pd.DataFrame(data, columns=columns)
    sampled_tstv_ratios = (
        rnd_spectra[transitions].sum(axis=1) / rnd_spectra[transversions].sum(axis=1)
    ).values

    if tstv_sampling == "poisson":
        tstv_ratios = np.random.poisson(3.5, n) + np.random.random(n)
        multiplier = tstv_ratios / sampled_tstv_ratios
        rnd_spectra[transitions] *= multiplier[:, None]
    elif tstv_sampling == "random":
        tstv_ratios = sampled_tstv_ratios
    else:
        raise ValueError(f"Unknown tstv_sampling: {tstv_sampling}")

    rnd_spectra = (rnd_spectra.T / rnd_spectra.sum(axis=1)).T
    rnd_spectra["tstv_ratio"] = tstv_ratios
    rnd_spectra.reset_index(names="replica", inplace=True)
    rnd_spectra_long = rnd_spectra.melt(
        ["replica", "tstv_ratio"], var_name="Mut", value_name="MutSpec"
    )
    return rnd_spectra_long


def get_random_spectrum(tstv_ratio=None):
    """One random 12-component spectrum. Typical ts/tv ratios are between 2 and 20."""
    rnd_spectrum = pd.DataFrame(
        {
            "Mut": [f"{n1}>{n2}" for n1 in alphabet for n2 in alphabet if n1 != n2],
            "MutSpec": uniform.rvs(size=12),
        }
    )
    if tstv_ratio is not None:
        ts_frac = rnd_spectrum[rnd_spectrum["Mut"].isin(transitions)].MutSpec.sum()
        tv_frac = rnd_spectrum[~rnd_spectrum["Mut"].isin(transitions)].MutSpec.sum()
        multiplier = tstv_ratio / (ts_frac / tv_frac)
        rnd_spectrum.loc[rnd_spectrum["Mut"].isin(transitions), "MutSpec"] *= multiplier

    rnd_spectrum["MutSpec"] = rnd_spectrum["MutSpec"] / rnd_spectrum["MutSpec"].sum()
    return rnd_spectrum


def get_equilibrium_freqs_rnd(gc=1, tstv_ratio=2.0):
    """Equilibrium frequencies from one random spectrum."""
    rnd_spectrum = get_random_spectrum(tstv_ratio=tstv_ratio)
    res = get_equilibrium_freqs(rnd_spectrum, rate_col="MutSpec", gc=gc)
    return rnd_spectrum, res


def prepare_exp_aa_subst(spectrum: pd.DataFrame, rate_col="rate", gc=1, save_path=None):
    df_changes = collect_possible_changes(gc=gc)
    spectrum_dict = spectrum.set_index("Mut")[rate_col].to_dict()
    df_changes["rate"] = df_changes["sbs"].map(spectrum_dict)

    per_sbs = (
        df_changes.query('aa1 != aa2 & aa1 != "*" & aa2 != "*"')
        .groupby(["aa1", "aa2", "sbs"], as_index=False)["rate"]
        .sum()
    )
    # Join mutation labels explicitly. Summing the string column used to
    # concatenate them, which current pandas drops when numeric_only is set.
    exp_aa_subst = per_sbs.groupby(["aa1", "aa2"], as_index=False).agg(
        rate=("rate", "sum"),
        sbs=("sbs", lambda values: "|".join(values)),
    )

    if save_path:
        exp_aa_subst.to_csv(save_path, float_format="%g", index=False)
    exp_aa_subst_matrix = exp_aa_subst.pivot(index="aa1", columns="aa2", values="rate").fillna(0.0)
    return exp_aa_subst, exp_aa_subst_matrix


def prepare_exp_cdn_subst(spectrum: pd.DataFrame, rate_col="rate", gc=1, save_path=None):
    df_changes = collect_possible_changes(gc=gc)
    spectrum_dict = spectrum.set_index("Mut")[rate_col].to_dict()
    df_changes["rate"] = df_changes["sbs"].map(spectrum_dict)
    exp_cdn_subst = df_changes.groupby(["cdn1", "cdn2"])["rate"].sum().reset_index()

    if save_path:
        exp_cdn_subst.to_csv(save_path, float_format="%g", index=False)
    exp_cdn_subst_matrix = exp_cdn_subst.pivot(index="cdn1", columns="cdn2", values="rate").fillna(0.0)
    return exp_cdn_subst, exp_cdn_subst_matrix


def prepare_exp_aa_subst_codons(
    spectrum: pd.DataFrame, rate_col="rate", gc=1, save_path=None, codon_pos="all"
):
    df_changes = collect_possible_changes(gc=gc)
    spectrum_dict = spectrum.set_index("Mut")[rate_col].to_dict()
    df_changes["rate"] = df_changes["sbs"].map(spectrum_dict)

    if codon_pos == "all":
        df_changes_flt = df_changes.query('aa1 != "*" & aa2 != "*"')
    elif codon_pos == "12":
        df_changes_flt = df_changes.query('aa1 != "*" & aa2 != "*" & pic in [0, 1]')
    elif isinstance(codon_pos, int):
        pic = codon_pos - 1
        df_changes_flt = df_changes.query('aa1 != "*" & aa2 != "*" & pic == @pic')
    else:
        raise ValueError(
            f"Invalid codon_pos value: {codon_pos}. Allowed values are 'all', '12', 1, 2, or 3."
        )

    exp_aa_subst = df_changes_flt.groupby(["aa1", "aa2"])["rate"].sum().reset_index()

    if save_path:
        exp_aa_subst.to_csv(save_path, float_format="%g", index=False)
    exp_aa_subst_matrix = exp_aa_subst.pivot(index="aa1", columns="aa2", values="rate").fillna(0.0)
    return exp_aa_subst, exp_aa_subst_matrix


def prepare_rnd_exp_aa_subst(gc=1, tstv_ratio=None, save_path=None):
    """Expected amino-acid substitutions from one random spectrum."""
    rnd_spectrum = get_random_spectrum(tstv_ratio=tstv_ratio)
    return prepare_exp_aa_subst(rnd_spectrum, rate_col="MutSpec", gc=gc, save_path=save_path)


def prepare_aa_subst(obs_df: pd.DataFrame, exp_aa_subst: pd.DataFrame, ref_aa_freqs: dict):
    aa_subst = obs_df.groupby(["aa1", "aa2"])["count"].sum().rename("nobs").reset_index()
    aa_subst = aa_subst[
        (aa_subst["aa1"] != aa_subst["aa2"])
        & (aa_subst["aa1"] != "*")
        & (aa_subst["aa2"] != "*")
    ]
    aa_subst["aa1"] = aa_subst["aa1"].map(amino_acid_codes)
    aa_subst["aa2"] = aa_subst["aa2"].map(amino_acid_codes)

    aa_subst = aa_subst.merge(exp_aa_subst.rename(columns={"rate": "rate_exp"}), "right").fillna(0)
    aa_subst = aa_subst[aa_subst["aa1"] != aa_subst["aa2"]]

    ref_aa_total_cnt = sum(x for x in ref_aa_freqs.values() if x > 0)
    aa_subst["ref_aa1_freq"] = aa_subst["aa1"].map(ref_aa_freqs) / ref_aa_total_cnt
    aa_subst["nobs_scaled"] = (aa_subst["nobs"] / aa_subst["ref_aa1_freq"]).replace(np.inf, np.nan)
    aa_subst["nobs_scaled"] = aa_subst["nobs_scaled"] / aa_subst["nobs_scaled"].sum() * aa_subst["nobs"].sum()

    aa_subst["nexp"] = aa_subst["rate_exp"] / aa_subst["rate_exp"].sum() * aa_subst["nobs_scaled"].sum()
    aa_subst["diff"] = aa_subst["nobs_scaled"] - aa_subst["nexp"]
    aa_subst["pe"] = aa_subst["diff"] / aa_subst["nobs_scaled"] * 100
    aa_subst["nobs_freqs"] = aa_subst["nobs_scaled"] / aa_subst["nobs_scaled"].sum()
    aa_subst["nexp_freqs"] = aa_subst["nexp"] / aa_subst["nexp"].sum()
    aa_subst["obs_relative_freq"] = aa_subst["nobs_freqs"] / aa_subst["nexp_freqs"]
    aa_subst = aa_subst.merge(grantham_long, "left")
    return aa_subst
