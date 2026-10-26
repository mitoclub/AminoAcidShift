"""Derive observed and equilibrium frequencies for Cytb and ND2.

Inputs: data/raw/MutSpecVertebrates12.csv.gz and
data/raw/output/{Gene}/{Gene}__{Species}/seqs_unique.fasta.
Outputs: data/nuc_freqs.csv, data/cdn_freqs.csv, and data/aa_freqs.csv.
Each table has columns Class, Gene, Species, Type, then one column per state.
Type is obs or eq, so each species has two rows. Nucleotides also have subset.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO
from pymutspec.annotation import CodonAnnotation
from tqdm import tqdm

from aminoacidshift import (
    alphabet,
    amino_acid_codes,
    cdn_spectrum_to_matrix,
    collect_possible_changes,
    get_equilibrium_probabilities,
    nuc_spectrum_to_matrix,
    prepare_exp_aa_subst,
)

HERE = Path(__file__).resolve().parent
SPECTRUM_PATH = HERE / "data" / "raw" / "MutSpecVertebrates12.csv.gz"
FASTA_ROOT = HERE / "data" / "raw" / "output"
DATA_DIR = HERE / "data"
FIGURE_DIR = HERE / "figures"
GENETIC_CODE = 2

GENE_ORDER = ["Cytb", "ND2"]
CLASS_ORDER = ["Mammalia", "Aves", "Lepidosauria", "Amphibia", "Actinopteri"]
SUBSET_ORDER = ["all", "pos1", "pos2", "pos3", "fourfold"]
SUBSET_LABELS = {
    "all": "CDS",
    "pos1": "1st",
    "pos2": "2nd",
    "pos3": "3rd",
    "fourfold": "Four-fold",
}
CLASS_PALETTE = {
    "Mammalia": "#4C78A8",
    "Aves": "#F58518",
    "Lepidosauria": "#54A24B",
    "Amphibia": "#E45756",
    "Actinopteri": "#72B7B2",
}
BASE_INDEX = {base: index for index, base in enumerate(alphabet)}
AA_ORDER = [name for letter, name in amino_acid_codes.items() if letter != "*"]


def stationary(matrix: np.ndarray) -> np.ndarray | None:
    """Equilibrium of a rate matrix, or None when it is not a distribution."""
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        return None
    if not np.isfinite(matrix).all():
        return None
    try:
        pi = get_equilibrium_probabilities(matrix.astype(float, copy=False))
    except (ValueError, np.linalg.LinAlgError):
        return None
    if pi.shape != (matrix.shape[0],) or not np.isfinite(pi).all():
        return None
    if np.any(pi < -1e-8):
        return None
    pi = np.clip(pi, 0.0, None)
    total = pi.sum()
    if total <= 0:
        return None
    return pi / total


def nucleotide_frequency(bases: str) -> np.ndarray | None:
    counts = np.zeros(len(alphabet), dtype=float)
    for base in bases:
        index = BASE_INDEX.get(base)
        if index is not None:
            counts[index] += 1
    total = counts.sum()
    if total == 0:
        return None
    return counts / total


def average(vectors: list[np.ndarray | None]) -> np.ndarray | None:
    usable = [vector for vector in vectors if vector is not None]
    if not usable:
        return None
    return np.mean(np.stack(usable), axis=0)


def sequence_content(
    seq: str,
    codon_index: dict[str, int],
    aa_index: dict[str, int],
    fourfold_codons: set[str],
    coda: CodonAnnotation,
) -> tuple[dict[str, np.ndarray | None], np.ndarray | None, np.ndarray | None]:
    seq = str(seq).upper()
    nucleotides = {name: None for name in SUBSET_ORDER}
    nucleotides["all"] = nucleotide_frequency(seq)
    framed = seq[: len(seq) - (len(seq) % 3)]
    codon_frequency = None
    aa_frequency = None
    if framed:
        nucleotides["pos1"] = nucleotide_frequency(framed[0::3])
        nucleotides["pos2"] = nucleotide_frequency(framed[1::3])
        nucleotides["pos3"] = nucleotide_frequency(framed[2::3])
        codon_counts = np.zeros(len(codon_index), dtype=float)
        aa_counts = np.zeros(len(aa_index), dtype=float)
        fourfold_bases = []
        for start in range(0, len(framed), 3):
            codon = framed[start : start + 3]
            codon_id = codon_index.get(codon)
            if codon_id is not None:
                codon_counts[codon_id] += 1
                amino_acid = amino_acid_codes[coda.translate_codon(codon)]
                if amino_acid != "*":
                    aa_counts[aa_index[amino_acid]] += 1
            if codon in fourfold_codons:
                fourfold_bases.append(codon[2])
        nucleotides["fourfold"] = nucleotide_frequency("".join(fourfold_bases))
        codon_total = codon_counts.sum()
        if codon_total:
            codon_frequency = codon_counts / codon_total
        aa_total = aa_counts.sum()
        if aa_total:
            aa_frequency = aa_counts / aa_total
    return nucleotides, codon_frequency, aa_frequency


def observed_content(path: Path, codon_order: list[str], aa_order: list[str], coda: CodonAnnotation):
    codon_index = {codon: index for index, codon in enumerate(codon_order)}
    aa_index = {amino_acid: index for index, amino_acid in enumerate(aa_order)}
    fourfold_codons = {codon for codon in codon_order if coda.is_fourfold(codon)}
    nucleotide_rows = {name: [] for name in SUBSET_ORDER}
    codon_rows = []
    aa_rows = []
    n_sequences = 0
    for record in SeqIO.parse(path, "fasta"):
        if record.id.upper() == "OUTGRP":
            continue
        n_sequences += 1
        nucleotides, codon_frequency, aa_frequency = sequence_content(
            record.seq, codon_index, aa_index, fourfold_codons, coda
        )
        for name, frequency in nucleotides.items():
            nucleotide_rows[name].append(frequency)
        codon_rows.append(codon_frequency)
        aa_rows.append(aa_frequency)
    observed_nucleotides = {
        name: average(vectors) for name, vectors in nucleotide_rows.items()
    }
    return observed_nucleotides, average(codon_rows), average(aa_rows), n_sequences


def codon_order_from_changes(changes: pd.DataFrame) -> list[str]:
    rated = changes.copy()
    rated["rate"] = 1.0
    return list(cdn_spectrum_to_matrix(rated).index)


def codon_equilibrium(changes: pd.DataFrame, spectrum: dict[str, float], codon_order: list[str]):
    rated = changes.copy()
    rated["rate"] = rated["sbs"].map(spectrum)
    matrix = cdn_spectrum_to_matrix(rated)
    aligned = matrix.reindex(index=codon_order, columns=codon_order)
    return stationary(aligned.to_numpy(dtype=float))


def amino_acid_equilibrium(spectrum: pd.DataFrame, aa_order: list[str]):
    _, matrix = prepare_exp_aa_subst(spectrum, rate_col="MutSpec", gc=GENETIC_CODE)
    aligned = matrix.reindex(index=aa_order, columns=aa_order, fill_value=0.0)
    return stationary(aligned.to_numpy(dtype=float))


def append_paired(rows, meta, state_name, states, observed, equilibrium, extra=None):
    if observed is None and equilibrium is None:
        return
    observed_values = observed if observed is not None else np.full(len(states), np.nan)
    equilibrium_values = equilibrium if equilibrium is not None else np.full(len(states), np.nan)
    for state, obs_freq, eq_freq in zip(states, observed_values, equilibrium_values):
        row = dict(meta)
        if extra:
            row.update(extra)
        row[state_name] = state
        row["obs_freq"] = obs_freq
        row["eq_freq"] = eq_freq
        rows.append(row)


def write_wide(rows, state_col: str, states: list[str], id_cols: list[str], path: Path) -> None:
    long = pd.DataFrame(rows)
    parts = []
    for kind, value_col in (("obs", "obs_freq"), ("eq", "eq_freq")):
        part = long.pivot(index=id_cols, columns=state_col, values=value_col).reindex(columns=states)
        part = part.reset_index()
        part.insert(len(id_cols), "Type", kind)
        parts.append(part)
    wide = pd.concat(parts, ignore_index=True)
    wide["_gene"] = wide["Gene"].map({gene: index for index, gene in enumerate(GENE_ORDER)})
    wide["_type"] = wide["Type"].map({"obs": 0, "eq": 1})
    sort_cols = ["_gene", "Class", "Species", "_type"]
    if "subset" in wide.columns:
        wide["_subset"] = wide["subset"].map({name: index for index, name in enumerate(SUBSET_ORDER)})
        sort_cols = ["_gene", "Class", "Species", "_subset", "_type"]
    wide = wide.sort_values(sort_cols).drop(columns=["_gene", "_type"] + (["_subset"] if "_subset" in wide.columns else []))
    wide.to_csv(path, index=False, float_format="%.6g")


def read_frequencies(path: Path, state_col: str) -> pd.DataFrame:
    """Melt a Type=obs/eq frequency table back to one row per state."""
    table = pd.read_csv(path)
    id_cols = [column for column in ("Class", "Gene", "Species", "subset") if column in table.columns]
    state_cols = [column for column in table.columns if column not in id_cols + ["Type"]]
    observed = table.loc[table["Type"] == "obs", id_cols + state_cols]
    equilibrium = table.loc[table["Type"] == "eq", id_cols + state_cols]
    observed = observed.melt(id_vars=id_cols, value_vars=state_cols, var_name=state_col, value_name="obs_freq")
    equilibrium = equilibrium.melt(id_vars=id_cols, value_vars=state_cols, var_name=state_col, value_name="eq_freq")
    return observed.merge(equilibrium, on=id_cols + [state_col], how="left")


def main() -> None:
    DATA_DIR.mkdir(exist_ok=True)
    coda = CodonAnnotation(GENETIC_CODE)
    changes = collect_possible_changes(GENETIC_CODE)
    states = {
        "nucleotide": list(alphabet),
        "codon": codon_order_from_changes(changes),
        "aa": AA_ORDER,
    }
    spectra = pd.read_csv(SPECTRUM_PATH)
    spectra = spectra[spectra["Gene"].isin(GENE_ORDER)]
    groups = list(spectra.groupby(["Gene", "Species"], sort=False))

    nuc_rows, cdn_rows, aa_rows = [], [], []
    n_missing_fasta = 0
    n_failed = {"nucleotide": 0, "codon": 0, "aa": 0}

    for (gene, species), group in tqdm(groups, desc="Species-gene"):
        meta = {
            "Class": group["Class"].iloc[0],
            "Gene": gene,
            "Species": species,
        }
        spectrum = group.set_index("Mut")["MutSpec"].to_dict()
        nuc_eq = stationary(nuc_spectrum_to_matrix(spectrum))
        cdn_eq = codon_equilibrium(changes, spectrum, states["codon"])
        aa_eq = amino_acid_equilibrium(group[["Mut", "MutSpec"]], states["aa"])
        for level, equilibrium in (
            ("nucleotide", nuc_eq),
            ("codon", cdn_eq),
            ("aa", aa_eq),
        ):
            if equilibrium is None:
                n_failed[level] += 1

        fasta_path = FASTA_ROOT / gene / f"{gene}__{species}" / "seqs_unique.fasta"
        if fasta_path.exists():
            nuc_obs, cdn_obs, aa_obs, _ = observed_content(
                fasta_path, states["codon"], states["aa"], coda
            )
        else:
            n_missing_fasta += 1
            nuc_obs = {name: None for name in SUBSET_ORDER}
            cdn_obs = None
            aa_obs = None

        for subset, observed in nuc_obs.items():
            append_paired(
                nuc_rows,
                meta,
                "nucleotide",
                states["nucleotide"],
                observed,
                nuc_eq,
                {"subset": subset},
            )
        append_paired(cdn_rows, meta, "cdn", states["codon"], cdn_obs, cdn_eq)
        append_paired(aa_rows, meta, "aa", states["aa"], aa_obs, aa_eq)

    write_wide(
        nuc_rows,
        "nucleotide",
        states["nucleotide"],
        ["Class", "Gene", "Species", "subset"],
        DATA_DIR / "nuc_freqs.csv",
    )
    write_wide(
        cdn_rows,
        "cdn",
        states["codon"],
        ["Class", "Gene", "Species"],
        DATA_DIR / "cdn_freqs.csv",
    )
    write_wide(
        aa_rows,
        "aa",
        states["aa"],
        ["Class", "Gene", "Species"],
        DATA_DIR / "aa_freqs.csv",
    )
    print(f"pairs: {len(groups)}")
    print(f"missing fasta: {n_missing_fasta}")
    print(f"failed equilibria: {n_failed}")


if __name__ == "__main__":
    main()
