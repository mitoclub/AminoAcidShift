"""Compare fish and bird content with equilibria from both classes.

Inputs: data/nuc_freqs.csv, data/cdn_freqs.csv, and data/aa_freqs.csv.
Outputs: data/fish_bird_distances.csv and figures/fish_bird_{nuc,cdn,aa}.pdf.
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
    DATA_DIR,
    FIGURE_DIR,
    GENE_ORDER,
    SUBSET_LABELS,
    SUBSET_ORDER,
    read_frequencies,
)

FOCAL_CLASSES = ["Actinopteri", "Aves"]
COMPARISON_PALETTE = {
    "own": "#4C78A8",
    "within": "#54A24B",
    "other": "#E45756",
}
KEYS = ["Gene", "focal_class", "Species"]


def state_order(frame: pd.DataFrame, state_col: str) -> list[str]:
    return list(dict.fromkeys(frame[state_col]))


def species_matrices(frame: pd.DataFrame, state_col: str, states: list[str]):
    species = list(dict.fromkeys(frame["Species"]))
    observed = (
        frame.pivot(index="Species", columns=state_col, values="obs_freq")
        .reindex(index=species, columns=states)
        .to_numpy(dtype=float)
    )
    equilibrium = (
        frame.pivot(index="Species", columns=state_col, values="eq_freq")
        .reindex(index=species, columns=states)
        .to_numpy(dtype=float)
    )
    return species, observed, equilibrium


def pair_table(
    gene: str,
    level: str,
    subset: str,
    focal_class: str,
    focal_species: list[str],
    focal_obs: np.ndarray,
    focal_eq: np.ndarray,
    reference_class: str,
    reference_species: list[str],
    reference_eq: np.ndarray,
    comparison: str,
) -> pd.DataFrame | None:
    focal_ok = np.isfinite(focal_obs).all(axis=1) & np.isfinite(focal_eq).all(axis=1)
    reference_ok = np.isfinite(reference_eq).all(axis=1)
    focal_index = np.flatnonzero(focal_ok)
    reference_index = np.flatnonzero(reference_ok)
    if focal_index.size == 0 or reference_index.size == 0:
        return None
    distances = np.abs(
        focal_obs[focal_index, None, :] - reference_eq[None, reference_index, :]
    ).sum(axis=2)
    if comparison == "own":
        keep = focal_index[:, None] == reference_index[None, :]
    elif comparison == "within":
        keep = focal_index[:, None] != reference_index[None, :]
    else:
        keep = np.ones(distances.shape, dtype=bool)
    focal_pairs, reference_pairs = np.nonzero(keep)
    if focal_pairs.size == 0:
        return None
    focal_names = np.asarray(focal_species, dtype=object)[focal_index[focal_pairs]]
    reference_names = np.asarray(reference_species, dtype=object)[reference_index[reference_pairs]]
    return pd.DataFrame(
        {
            "Gene": gene,
            "level": level,
            "subset": subset,
            "focal_class": focal_class,
            "Species": focal_names,
            "reference_class": reference_class,
            "reference_species": reference_names,
            "comparison": comparison,
            "distance": distances[focal_pairs, reference_pairs],
        }
    )


def contrasts_for_slice(
    gene: str,
    level: str,
    subset: str,
    by_class: dict[str, tuple[list[str], np.ndarray, np.ndarray]],
) -> list[pd.DataFrame]:
    blocks = []
    for focal_class in FOCAL_CLASSES:
        other_class = next(name for name in FOCAL_CLASSES if name != focal_class)
        if focal_class not in by_class or other_class not in by_class:
            continue
        focal_species, focal_obs, focal_eq = by_class[focal_class]
        other_species, _, other_eq = by_class[other_class]
        for comparison, reference_class, reference_species, reference_eq in (
            ("own", focal_class, focal_species, focal_eq),
            ("within", focal_class, focal_species, focal_eq),
            ("other", other_class, other_species, other_eq),
        ):
            block = pair_table(
                gene,
                level,
                subset,
                focal_class,
                focal_species,
                focal_obs,
                focal_eq,
                reference_class,
                reference_species,
                reference_eq,
                comparison,
            )
            if block is not None:
                blocks.append(block)
    return blocks


def build_contrasts(frame: pd.DataFrame, state_col: str, level: str, subset: str) -> pd.DataFrame:
    states = state_order(frame, state_col)
    blocks = []
    for gene in GENE_ORDER:
        gene_frame = frame[frame["Gene"] == gene]
        if subset:
            gene_frame = gene_frame[gene_frame["subset"] == subset]
        by_class = {}
        for class_name in FOCAL_CLASSES:
            class_frame = gene_frame[gene_frame["Class"] == class_name]
            if class_frame.empty:
                continue
            by_class[class_name] = species_matrices(class_frame, state_col, states)
        blocks.extend(contrasts_for_slice(gene, level, subset, by_class))
    if not blocks:
        return pd.DataFrame(
            columns=[
                "Gene",
                "level",
                "subset",
                "focal_class",
                "Species",
                "reference_class",
                "reference_species",
                "comparison",
                "distance",
            ]
        )
    return pd.concat(blocks, ignore_index=True)


def species_medians(frame: pd.DataFrame) -> pd.DataFrame:
    own = frame.loc[frame["comparison"] == "own", [*KEYS, "distance"]].rename(
        columns={"distance": "own"}
    )
    medians = (
        frame.loc[frame["comparison"] != "own"]
        .groupby([*KEYS, "comparison"], sort=False)["distance"]
        .median()
        .reset_index()
    )
    return own.merge(medians, on=KEYS, how="inner")


def report(frame: pd.DataFrame, level: str, subset: str) -> None:
    points = species_medians(frame)
    label = level if not subset else f"{level} {subset}"
    for class_name in FOCAL_CLASSES:
        class_points = points[points["focal_class"] == class_name]
        if class_points.empty:
            print(f"{label} {class_name}: no species")
            continue
        wide = class_points.pivot(index=["Gene", "Species"], columns="comparison", values="distance")
        if not {"within", "other"}.issubset(wide.columns):
            print(f"{label} {class_name}: no paired medians")
            continue
        usable = wide.dropna(subset=["within", "other"])
        if usable.empty:
            print(f"{label} {class_name}: no paired medians")
            continue
        own = class_points.drop_duplicates(["Gene", "Species"]).set_index(["Gene", "Species"])["own"]
        closer_within = own.loc[usable.index] < usable["within"]
        closer_other = own.loc[usable.index] < usable["other"]
        closer_both = closer_within & closer_other
        print(
            f"{label} {class_name}: "
            f"own < within {closer_within.mean():.3f}, "
            f"own < other {closer_other.mean():.3f}, "
            f"own < both {closer_both.mean():.3f} "
            f"({int(closer_both.sum())}/{len(usable)})"
        )


def plot_identity(frame: pd.DataFrame, title: str):
    points = species_medians(frame)
    if points.empty:
        return None
    grid = sns.relplot(
        data=points,
        x="own",
        y="distance",
        hue="comparison",
        col="Gene",
        row="focal_class",
        col_order=[gene for gene in GENE_ORDER if gene in set(points["Gene"])],
        row_order=[name for name in FOCAL_CLASSES if name in set(points["focal_class"])],
        hue_order=["within", "other"],
        palette=COMPARISON_PALETTE,
        kind="scatter",
        height=3.2,
        aspect=0.95,
        s=18,
        linewidth=0,
        facet_kws={"sharex": False, "sharey": False},
    )
    for axis in grid.axes.flat:
        if not axis.collections and not axis.lines:
            continue
        values = []
        for collection in axis.collections:
            offsets = collection.get_offsets()
            if len(offsets):
                values.append(offsets)
        if not values:
            continue
        stacked = np.vstack(values)
        finite = stacked[np.isfinite(stacked).all(axis=1)]
        if finite.size == 0:
            continue
        limit = float(np.max(finite)) * 1.05
        axis.plot([0, limit], [0, limit], color="0.45", lw=0.8, zorder=0)
        axis.set_xlim(0, limit)
        axis.set_ylim(0, limit)
    grid.set_axis_labels("Distance to own equilibrium", "Median distance")
    grid.set_titles("{row_name} | {col_name}")
    grid.figure.suptitle(title, y=1.03)
    if grid.legend is not None:
        grid.legend.set_title("Comparison")
    grid.tight_layout()
    return grid


def add_page(grid, pdf) -> None:
    if grid is None:
        return
    pdf.savefig(grid.figure, bbox_inches="tight")
    plt.close(grid.figure)


def write_table(frame: pd.DataFrame, path: Path, first: bool) -> None:
    frame.to_csv(
        path,
        mode="w" if first else "a",
        header=first,
        index=False,
        float_format="%.6g",
    )


def main() -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    output = DATA_DIR / "fish_bird_distances.csv"
    first = True
    nucleotides = read_frequencies(DATA_DIR / "nuc_freqs.csv", "nucleotide")
    with PdfPages(FIGURE_DIR / "fish_bird_nuc.pdf") as pdf:
        for subset in SUBSET_ORDER:
            print(f"nucleotide {subset}", flush=True)
            frame = build_contrasts(nucleotides, "nucleotide", "nuc", subset)
            write_table(frame, output, first)
            first = False
            add_page(plot_identity(frame, f"Nucleotide distances, {SUBSET_LABELS[subset]}"), pdf)
            report(frame, "nuc", subset)
            del frame

    for level, state_col, title in (
        ("cdn", "cdn", "Codon"),
        ("aa", "aa", "Amino acid"),
    ):
        print(level, flush=True)
        frame = build_contrasts(read_frequencies(DATA_DIR / f"{level}_freqs.csv", state_col), state_col, level, "")
        write_table(frame, output, first)
        first = False
        with PdfPages(FIGURE_DIR / f"fish_bird_{level}.pdf") as pdf:
            add_page(plot_identity(frame, f"{title} distances"), pdf)
        report(frame, level, "")


if __name__ == "__main__":
    main()
