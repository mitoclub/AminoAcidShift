"""PCA of Cytb and ND2 frequencies for birds and fish.

Inputs: data/nuc_freqs.csv, data/cdn_freqs.csv, and data/aa_freqs.csv.
Outputs: figures/pca_aa_freqs.pdf, figures/pca_cdn_freqs.pdf, and figures/pca_nuc_freqs.pdf.
The decomposition is fit on observed frequencies; equilibrium is projected into it.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA

from derive_freqs import CLASS_PALETTE, DATA_DIR, FIGURE_DIR, read_frequencies

GENES = ["Cytb", "ND2"]
CLASSES = ["Aves", "Actinopteri"]
LEVELS = (
    ("aa", "aa", "aa_freqs.csv", None, "Amino acid"),
    ("cdn", "cdn", "cdn_freqs.csv", None, "Codon"),
    ("nuc", "nucleotide", "nuc_freqs.csv", "all", "CDS nucleotide"),
)


def frequency_matrix(frame: pd.DataFrame, state_col: str, value_col: str) -> pd.DataFrame:
    states = list(dict.fromkeys(frame[state_col]))
    wide = frame.pivot(index=["Class", "Gene", "Species"], columns=state_col, values=value_col)
    wide = wide.reindex(columns=states)
    complete = np.isfinite(wide.to_numpy(dtype=float)).all(axis=1)
    return wide.loc[complete]


def scores_for(observed: pd.DataFrame, equilibrium: pd.DataFrame, genes: list[str]):
    observed = observed[observed.index.get_level_values("Gene").isin(genes)]
    equilibrium = equilibrium[equilibrium.index.get_level_values("Gene").isin(genes)]
    pca = PCA(n_components=2, random_state=0)
    observed_scores = pca.fit_transform(observed.to_numpy(dtype=float))
    shared = equilibrium.index.intersection(observed.index)
    equilibrium_scores = pca.transform(equilibrium.loc[shared].to_numpy(dtype=float))
    observed_frame = observed.index.to_frame(index=False)
    observed_frame["PC1"] = observed_scores[:, 0]
    observed_frame["PC2"] = observed_scores[:, 1]
    observed_frame["source"] = "Observed"
    equilibrium_frame = shared.to_frame(index=False)
    equilibrium_frame["PC1"] = equilibrium_scores[:, 0]
    equilibrium_frame["PC2"] = equilibrium_scores[:, 1]
    equilibrium_frame["source"] = "Equilibrium"
    return pca, pd.concat([observed_frame, equilibrium_frame], ignore_index=True)


def draw_panel(axis, frame: pd.DataFrame, title: str, explained: np.ndarray) -> None:
    for class_name in CLASSES:
        class_frame = frame[frame["Class"] == class_name]
        for source, marker, size in (
            ("Observed", "o", 16),
            ("Equilibrium", "X", 22),
        ):
            points = class_frame[class_frame["source"] == source]
            if points.empty:
                continue
            axis.scatter(
                points["PC1"],
                points["PC2"],
                c=CLASS_PALETTE[class_name],
                marker=marker,
                s=size,
                alpha=0.65 if source == "Observed" else 0.85,
                linewidths=0,
                rasterized=True,
            )
    axis.set_title(title)
    axis.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
    axis.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
    axis.axhline(0, color="0.85", lw=0.6, zorder=0)
    axis.axvline(0, color="0.85", lw=0.6, zorder=0)


def legend_handles():
    class_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=CLASS_PALETTE[name],
            markeredgecolor="none",
            markersize=6,
            label=name,
        )
        for name in CLASSES
    ]
    source_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="0.25",
            markeredgecolor="none",
            markersize=6,
            label="Observed",
        ),
        Line2D(
            [0],
            [0],
            marker="X",
            color="none",
            markerfacecolor="0.25",
            markeredgecolor="none",
            markersize=7,
            label="Equilibrium",
        ),
    ]
    return class_handles + source_handles


def plot_level(observed: pd.DataFrame, equilibrium: pd.DataFrame, title: str, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.1), constrained_layout=True)
    for axis, gene in zip(axes, GENES):
        pca, scores = scores_for(observed, equilibrium, [gene])
        draw_panel(axis, scores, gene, pca.explained_variance_ratio_)
        print(
            f"{title} {gene}: "
            f"PC1 {pca.explained_variance_ratio_[0] * 100:.1f}%, "
            f"PC2 {pca.explained_variance_ratio_[1] * 100:.1f}%, "
            f"observed {int((scores.source == 'Observed').sum())}, "
            f"equilibrium {int((scores.source == 'Equilibrium').sum())}",
            flush=True,
        )
    fig.legend(
        handles=legend_handles(),
        loc="center left",
        bbox_to_anchor=(1.005, 0.5),
        frameon=False,
        title=title,
    )
    fig.suptitle(f"{title} frequencies, PCA fit on observed frequencies", fontsize=13)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    for level, state_col, filename, subset, title in LEVELS:
        frame = read_frequencies(DATA_DIR / filename, state_col)
        frame = frame[frame["Gene"].isin(GENES) & frame["Class"].isin(CLASSES)]
        if subset is not None:
            frame = frame[frame["subset"] == subset]
        observed = frequency_matrix(frame, state_col, "obs_freq")
        equilibrium = frequency_matrix(frame, state_col, "eq_freq")
        plot_level(observed, equilibrium, title, FIGURE_DIR / f"pca_{level}_freqs.pdf")


if __name__ == "__main__":
    main()
