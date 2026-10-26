"""Compare observed content with each species' own equilibrium.

Inputs: data/nuc_freqs.csv, data/cdn_freqs.csv, and data/aa_freqs.csv.
Outputs: data/distances.csv and the distance and frequency figures in figures/.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

from derive_freqs import (
    AA_ORDER,
    CLASS_ORDER,
    CLASS_PALETTE,
    DATA_DIR,
    FIGURE_DIR,
    GENE_ORDER,
    SUBSET_LABELS,
    SUBSET_ORDER,
    alphabet,
    read_frequencies,
)

NUC_DISTANCE_COLUMNS = {
    "all": "nuc_all",
    "pos1": "nuc_pos1",
    "pos2": "nuc_pos2",
    "pos3": "nuc_pos3",
    "fourfold": "nuc_fourfold",
}
DISTANCE_COLUMNS = list(NUC_DISTANCE_COLUMNS.values()) + ["cdn", "aa"]
SOURCE_ORDER = ["Observed", "Equilibrium"]
SOURCE_PALETTE = {"Observed": "#4C78A8", "Equilibrium": "#F58518"}


def vector_distance(frame: pd.DataFrame, state_col: str, id_cols: list[str]) -> pd.Series:
    observed = frame.pivot(index=id_cols, columns=state_col, values="obs_freq")
    equilibrium = frame.pivot(index=id_cols, columns=state_col, values="eq_freq")
    equilibrium = equilibrium.reindex(index=observed.index, columns=observed.columns)
    values = np.abs(observed.to_numpy(dtype=float) - equilibrium.to_numpy(dtype=float))
    ok = np.isfinite(values).all(axis=1)
    distance = np.full(len(observed), np.nan)
    distance[ok] = values[ok].sum(axis=1)
    return pd.Series(distance, index=observed.index, name="manhattan")


def distance_table(nucleotides: pd.DataFrame, codons: pd.DataFrame, amino_acids: pd.DataFrame) -> pd.DataFrame:
    keys = ["Class", "Gene", "Species"]
    nucleotide_distance = vector_distance(
        nucleotides, "nucleotide", [*keys, "subset"]
    ).unstack("subset")
    nucleotide_distance = nucleotide_distance.reindex(columns=SUBSET_ORDER)
    nucleotide_distance.columns = [NUC_DISTANCE_COLUMNS[name] for name in nucleotide_distance.columns]
    codon_distance = vector_distance(codons, "cdn", keys).rename("cdn")
    amino_acid_distance = vector_distance(amino_acids, "aa", keys).rename("aa")
    distances = nucleotide_distance.join(codon_distance).join(amino_acid_distance).reset_index()
    return distances[["Class", "Gene", "Species", *DISTANCE_COLUMNS]]


def plot_frequency_pages(frame, state_col, state_order, path: Path, subset_col=None) -> None:
    pages = []
    classes = [name for name in CLASS_ORDER if name in set(frame["Class"])]
    genes = [name for name in GENE_ORDER if name in set(frame["Gene"])]
    subsets = [None]
    if subset_col is not None:
        subsets = [name for name in SUBSET_ORDER if name in set(frame[subset_col])]
    for class_name in classes:
        for gene in genes:
            for subset in subsets:
                page = frame[(frame["Class"] == class_name) & (frame["Gene"] == gene)]
                if subset is not None:
                    page = page[page[subset_col] == subset]
                if page.empty:
                    continue
                pages.append((class_name, gene, subset, page))
    if not pages:
        return
    with PdfPages(path) as pdf:
        for class_name, gene, subset, page in pages:
            long = page.melt(
                id_vars=[state_col],
                value_vars=["obs_freq", "eq_freq"],
                var_name="source",
                value_name="freq",
            )
            long["source"] = long["source"].map(
                {"obs_freq": "Observed", "eq_freq": "Equilibrium"}
            )
            long = long.dropna(subset=["freq"])
            if long.empty:
                continue
            width = max(6, 0.28 * len(state_order))
            fig, axis = plt.subplots(figsize=(width, 4.5))
            sns.boxplot(
                data=long,
                x=state_col,
                y="freq",
                hue="source",
                order=state_order,
                hue_order=SOURCE_ORDER,
                palette=SOURCE_PALETTE,
                fliersize=1.5,
                ax=axis,
            )
            subset_label = "" if subset is None else f", {SUBSET_LABELS[subset]}"
            axis.set_title(f"{class_name} {gene}{subset_label}")
            axis.set_xlabel("")
            axis.set_ylabel("Frequency")
            axis.tick_params(axis="x", rotation=90)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def plot_distance_pairs(distances: pd.DataFrame, path: Path) -> None:
    pairs = [
        ("nuc_all", "cdn", "CDS", "Codon"),
        ("nuc_all", "aa", "CDS", "Amino acid"),
        ("cdn", "aa", "Codon", "Amino acid"),
        ("nuc_pos1", "nuc_all", "1st", "CDS"),
        ("nuc_pos2", "nuc_all", "2nd", "CDS"),
        ("nuc_pos3", "nuc_all", "3rd", "CDS"),
        ("nuc_fourfold", "nuc_all", "Four-fold", "CDS"),
    ]
    genes = [gene for gene in GENE_ORDER if gene in set(distances["Gene"])]
    fig, axes = plt.subplots(
        len(genes),
        len(pairs),
        figsize=(3.1 * len(pairs), 3.1 * len(genes)),
        squeeze=False,
        sharex="col",
        sharey="col",
    )
    legend_axis = None
    for row, gene in enumerate(genes):
        subset = distances[distances["Gene"] == gene]
        for col, (x_col, y_col, x_label, y_label) in enumerate(pairs):
            axis = axes[row, col]
            sns.scatterplot(
                data=subset,
                x=x_col,
                y=y_col,
                hue="Class",
                hue_order=CLASS_ORDER,
                palette=CLASS_PALETTE,
                ax=axis,
                s=12,
                linewidth=0,
                legend=row == 0 and col == 0,
            )
            axis.set_xlabel(x_label if row == len(genes) - 1 else "")
            axis.set_ylabel(f"{gene}\n{y_label}" if col == 0 else y_label)
            if row == 0 and col == 0:
                legend_axis = axis
    if legend_axis is not None and legend_axis.get_legend() is not None:
        sns.move_legend(legend_axis, "upper left", bbox_to_anchor=(1.02, 1), title="Class")
    fig.suptitle("Pairwise Manhattan distances to each species' own equilibrium", y=1.01)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


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
    FIGURE_DIR.mkdir(exist_ok=True)
    nucleotides = read_frequencies(DATA_DIR / "nuc_freqs.csv", "nucleotide")
    codons = read_frequencies(DATA_DIR / "cdn_freqs.csv", "cdn")
    amino_acids = read_frequencies(DATA_DIR / "aa_freqs.csv", "aa")
    distances = distance_table(nucleotides, codons, amino_acids)
    distances.to_csv(DATA_DIR / "distances.csv", index=False, float_format="%.6g")

    nuc_plot = distances.melt(
        id_vars=["Class", "Gene", "Species"],
        value_vars=list(NUC_DISTANCE_COLUMNS.values()),
        var_name="subset",
        value_name="manhattan",
    )
    nuc_plot["subset"] = nuc_plot["subset"].map(
        {column: SUBSET_LABELS[subset] for subset, column in NUC_DISTANCE_COLUMNS.items()}
    )
    plot_distances(
        nuc_plot,
        "subset",
        [SUBSET_LABELS[name] for name in SUBSET_ORDER],
        "Nucleotide content versus mutational equilibrium",
        FIGURE_DIR / "nuc_eq_distance_distributions.pdf",
    )
    for column, title, filename in (
        ("cdn", "Codon content versus mutational equilibrium", "cdn_eq_distance_distributions.pdf"),
        ("aa", "Amino-acid content versus mutational equilibrium", "aa_eq_distance_distributions.pdf"),
    ):
        level_plot = distances.melt(
            id_vars=["Class", "Gene", "Species"],
            value_vars=[column],
            value_name="manhattan",
        )
        plot_distances(level_plot, "Class", CLASS_ORDER, title, FIGURE_DIR / filename)
    plot_distance_pairs(distances, FIGURE_DIR / "distance_pairwise_dotplots.pdf")
    plot_frequency_pages(
        nucleotides,
        "nucleotide",
        list(alphabet),
        FIGURE_DIR / "nuc_freq_by_class_gene.pdf",
        subset_col="subset",
    )
    plot_frequency_pages(
        codons,
        "cdn",
        list(dict.fromkeys(codons["cdn"])),
        FIGURE_DIR / "cdn_freq_by_class_gene.pdf",
    )
    plot_frequency_pages(
        amino_acids,
        "aa",
        AA_ORDER,
        FIGURE_DIR / "aa_freq_by_class_gene.pdf",
    )


if __name__ == "__main__":
    main()
