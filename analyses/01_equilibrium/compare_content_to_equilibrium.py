"""Compare species- and gene-specific content with mutational equilibrium.

For each vertebrate mitochondrial spectrum, estimate nucleotide, codon, and
amino-acid equilibrium and compare it with the coding sequences used to build
that gene tree. The outgroup record is excluded. Observed frequencies are the
unweighted mean across the remaining sequences.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
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

GENE_ORDER = ["CO1", "CO3", "Cytb", "ND2"]
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
    pi = stationary(aligned.to_numpy(dtype=float))
    return pi


def amino_acid_equilibrium(spectrum: pd.DataFrame, aa_order: list[str]):
    _, matrix = prepare_exp_aa_subst(spectrum, rate_col="MutSpec", gc=GENETIC_CODE)
    aligned = matrix.reindex(index=aa_order, columns=aa_order, fill_value=0.0)
    return stationary(aligned.to_numpy(dtype=float))


def manhattan(observed: np.ndarray | None, equilibrium: np.ndarray | None) -> float:
    if observed is None or equilibrium is None:
        return np.nan
    return float(np.abs(observed - equilibrium).sum())


def append_frequencies(rows, meta, state_name, states, values, value_name, extra=None):
    if values is None:
        return
    for state, value in zip(states, values):
        row = dict(meta)
        if extra:
            row.update(extra)
        row[state_name] = state
        row[value_name] = value
        rows.append(row)


def plot_distances(frame: pd.DataFrame, x: str, x_order: list[str], title: str, path: Path) -> None:
    plot_frame = frame.replace([np.inf, -np.inf], np.nan).dropna(subset=["manhattan"])
    if plot_frame.empty:
        return
    legend = x != "Class"
    grid = sns.catplot(
        data=plot_frame,
        x=x,
        y="manhattan",
        hue="Class",
        col="Gene",
        col_order=[gene for gene in GENE_ORDER if gene in set(plot_frame["Gene"])],
        order=x_order,
        hue_order=[name for name in CLASS_ORDER if name in set(plot_frame["Class"])],
        kind="box",
        palette=CLASS_PALETTE,
        height=3.6,
        aspect=1.05,
        fliersize=1.5,
        linewidth=0.6,
        sharey=True,
        legend=legend,
    )
    grid.set_axis_labels("", "Manhattan distance")
    grid.set_titles("{col_name}")
    grid.figure.suptitle(title, y=1.03)
    for axis in grid.axes.flat:
        axis.tick_params(axis="x", rotation=30)
    grid.savefig(path)
    plt.close(grid.figure)


def main() -> None:
    DATA_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)
    coda = CodonAnnotation(GENETIC_CODE)
    changes = collect_possible_changes(GENETIC_CODE)
    states = {
        "nucleotide": list(alphabet),
        "codon": codon_order_from_changes(changes),
        "aa": AA_ORDER,
    }
    spectra = pd.read_csv(SPECTRUM_PATH)
    groups = list(spectra.groupby(["Gene", "Species"], sort=False))

    nuc_eq_rows, nuc_obs_rows, nuc_dist_rows = [], [], []
    cdn_eq_rows, cdn_obs_rows, cdn_dist_rows = [], [], []
    aa_eq_rows, aa_obs_rows, aa_dist_rows = [], [], []
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

        append_frequencies(nuc_eq_rows, meta, "nucleotide", states["nucleotide"], nuc_eq, "eq_freq")
        append_frequencies(cdn_eq_rows, meta, "cdn", states["codon"], cdn_eq, "eq_freq")
        append_frequencies(aa_eq_rows, meta, "aa", states["aa"], aa_eq, "eq_freq")
        append_frequencies(cdn_obs_rows, meta, "cdn", states["codon"], cdn_obs, "obs_freq")
        append_frequencies(aa_obs_rows, meta, "aa", states["aa"], aa_obs, "obs_freq")
        for subset, observed in nuc_obs.items():
            append_frequencies(
                nuc_obs_rows,
                meta,
                "nucleotide",
                states["nucleotide"],
                observed,
                "obs_freq",
                {"subset": subset},
            )
            nuc_dist_rows.append({**meta, "subset": subset, "manhattan": manhattan(observed, nuc_eq)})
        cdn_dist_rows.append({**meta, "manhattan": manhattan(cdn_obs, cdn_eq)})
        aa_dist_rows.append({**meta, "manhattan": manhattan(aa_obs, aa_eq)})

    tables = {
        "nuc_eq_freqs.csv": nuc_eq_rows,
        "nuc_obs_freqs.csv": nuc_obs_rows,
        "nuc_distances.csv": nuc_dist_rows,
        "cdn_eq_freqs.csv": cdn_eq_rows,
        "cdn_obs_freqs.csv": cdn_obs_rows,
        "cdn_distances.csv": cdn_dist_rows,
        "aa_eq_freqs.csv": aa_eq_rows,
        "aa_obs_freqs.csv": aa_obs_rows,
        "aa_distances.csv": aa_dist_rows,
    }
    frames = {}
    for name, rows in tables.items():
        frame = pd.DataFrame(rows)
        frame.to_csv(DATA_DIR / name, index=False, float_format="%.6g")
        frames[name] = frame

    nuc_plot = frames["nuc_distances.csv"].copy()
    nuc_plot["subset"] = nuc_plot["subset"].map(SUBSET_LABELS)
    plot_distances(
        nuc_plot,
        "subset",
        [SUBSET_LABELS[name] for name in SUBSET_ORDER],
        "Nucleotide content versus mutational equilibrium",
        FIGURE_DIR / "nuc_eq_distance_distributions.pdf",
    )
    plot_distances(
        frames["cdn_distances.csv"],
        "Class",
        CLASS_ORDER,
        "Codon content versus mutational equilibrium",
        FIGURE_DIR / "cdn_eq_distance_distributions.pdf",
    )
    plot_distances(
        frames["aa_distances.csv"],
        "Class",
        CLASS_ORDER,
        "Amino-acid content versus mutational equilibrium",
        FIGURE_DIR / "aa_eq_distance_distributions.pdf",
    )
    print(f"pairs: {len(groups)}")
    print(f"missing fasta: {n_missing_fasta}")
    print(f"failed equilibria: {n_failed}")


if __name__ == "__main__":
    main()
