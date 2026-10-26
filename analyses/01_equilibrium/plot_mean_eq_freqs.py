"""Redraw the Cytb equilibrium summaries.

Inputs: data/aa_freqs.csv, data/nuc_freqs.csv, and
data/raw/MutSpecVertebrates12.csv.gz.
Outputs: figures/aa_eq_freqs_vert_cytb.pdf, figures/aa_eq_freqs_cls_heatmap.pdf,
figures/nuc_eq_freqs_cls.pdf, figures/aa_eq_obs_freqs_5cls.pdf, and
figures/exp_aa_subst_vert_matrix.pdf.
With --average spectra the same stems gain a _spectra suffix.

--average equilibrium (default): unweighted mean of per-species equilibria.
The substitution matrix is the mean of per-species matrices.
--average spectra: unweighted mean of mutation spectra, then the equilibrium
and substitution matrix of that mean spectrum.
"""

import argparse

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import seaborn.objects as so

from aminoacidshift import alphabet, nuc_spectrum_to_matrix, plot_exp_heatmap, prepare_exp_aa_subst
from derive_freqs import (
    AA_ORDER,
    DATA_DIR,
    FIGURE_DIR,
    GENETIC_CODE,
    SPECTRUM_PATH,
    amino_acid_equilibrium,
    stationary,
)

GENE = "Cytb"
CLASS_PLOT_ORDER = ["Actinopteri", "Amphibia", "Lepidosauria", "Mammalia", "Aves"]
ID_COLUMNS = {"Class", "Gene", "Species", "subset", "Type"}
MUTATIONS = [f"{source}>{target}" for source in alphabet for target in alphabet if source != target]
AVERAGE_EQUILIBRIUM = "equilibrium"
AVERAGE_SPECTRA = "spectra"


def state_columns(table: pd.DataFrame) -> list[str]:
    return [column for column in table.columns if column not in ID_COLUMNS]


def cytb_rows(table: pd.DataFrame) -> pd.DataFrame:
    rows = table[table["Gene"] == GENE]
    if "subset" in rows.columns:
        rows = rows[rows["subset"] == "all"]
    return rows


def species_equilibria(table: pd.DataFrame, states: list[str]) -> pd.DataFrame:
    rows = cytb_rows(table)
    equilibrium = rows.loc[rows["Type"] == "eq", ["Class", "Species", *states]]
    return equilibrium.dropna(subset=states)


def observations(table: pd.DataFrame, states: list[str], state_name: str) -> pd.DataFrame:
    rows = cytb_rows(table)
    observed = rows.loc[rows["Type"] == "obs", ["Class", "Species", *states]]
    return observed.melt(
        id_vars=["Class", "Species"],
        value_vars=states,
        var_name=state_name,
        value_name="freq",
    )


def output_path(filename: str, average: str):
    path = FIGURE_DIR / filename
    if average == AVERAGE_EQUILIBRIUM:
        return path
    return path.with_stem(f"{path.stem}_{average}")


def cytb_spectra() -> pd.DataFrame:
    spectra = pd.read_csv(SPECTRUM_PATH)
    return spectra[spectra["Gene"] == GENE]


def mean_spectrum(spectra: pd.DataFrame) -> pd.Series:
    mean = spectra.groupby("Mut")["MutSpec"].mean().reindex(MUTATIONS)
    if mean.isna().any():
        missing = ", ".join(mean.index[mean.isna()])
        raise RuntimeError(f"Mean spectrum is missing {missing}")
    return mean


def equilibrium_from_spectrum(spectrum: pd.Series) -> tuple[pd.Series, pd.Series]:
    spec = spectrum.to_dict()
    nucleotide = stationary(nuc_spectrum_to_matrix(spec))
    frame = spectrum.rename_axis("Mut").reset_index(name="MutSpec")
    amino_acid = amino_acid_equilibrium(frame, AA_ORDER)
    if nucleotide is None or amino_acid is None:
        raise RuntimeError("Mean spectrum has no equilibrium")
    return pd.Series(nucleotide, index=list(alphabet)), pd.Series(amino_acid, index=AA_ORDER)


def class_means_from_spectra(spectra: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series]:
    nucleotide_rows = []
    amino_acid_rows = []
    for class_name in CLASS_PLOT_ORDER:
        nucleotide, amino_acid = equilibrium_from_spectrum(mean_spectrum(spectra[spectra["Class"] == class_name]))
        nucleotide_rows.append(nucleotide.rename(class_name))
        amino_acid_rows.append(amino_acid.rename(class_name))
    nucleotide_vertebrate, amino_acid_vertebrate = equilibrium_from_spectrum(mean_spectrum(spectra))
    return (
        pd.DataFrame(nucleotide_rows),
        pd.DataFrame(amino_acid_rows),
        nucleotide_vertebrate,
        amino_acid_vertebrate,
    )


def class_means_from_equilibria(nucleotides: pd.DataFrame, amino_acids: pd.DataFrame):
    nucleotide_states = state_columns(nucleotides)
    amino_acid_states = state_columns(amino_acids)
    nucleotide_species = species_equilibria(nucleotides, nucleotide_states)
    amino_acid_species = species_equilibria(amino_acids, amino_acid_states)
    nucleotide_classes = nucleotide_species.groupby("Class")[nucleotide_states].mean().loc[CLASS_PLOT_ORDER]
    amino_acid_classes = amino_acid_species.groupby("Class")[amino_acid_states].mean().loc[CLASS_PLOT_ORDER]
    print(f"vertebrate amino acids: {len(amino_acid_species)} species")
    return (
        nucleotide_classes,
        amino_acid_classes,
        nucleotide_species[nucleotide_states].mean(),
        amino_acid_species[amino_acid_states].mean(),
    )


def substitution_matrix(spectra: pd.DataFrame, average: str) -> pd.DataFrame:
    if average == AVERAGE_SPECTRA:
        frame = mean_spectrum(spectra).rename_axis("Mut").reset_index(name="MutSpec")
        _, matrix = prepare_exp_aa_subst(frame, rate_col="MutSpec", gc=GENETIC_CODE)
        aligned = matrix.reindex(index=AA_ORDER, columns=AA_ORDER).fillna(0.0)
        print("substitution matrix: mean spectrum")
    else:
        matrices = []
        for _, spectrum in spectra.groupby("Species", sort=False):
            if spectrum["MutSpec"].notna().sum() < 12 or spectrum["Mut"].nunique() < 12:
                continue
            _, matrix = prepare_exp_aa_subst(spectrum[["Mut", "MutSpec"]], rate_col="MutSpec", gc=GENETIC_CODE)
            if not np.isfinite(matrix.to_numpy(dtype=float)).all():
                continue
            aligned = matrix.reindex(index=AA_ORDER, columns=AA_ORDER).fillna(0.0)
            matrices.append(aligned.to_numpy(dtype=float))
        if not matrices:
            raise RuntimeError("No complete Cytb spectra were available for the substitution matrix")
        aligned = pd.DataFrame(np.mean(matrices, axis=0), index=AA_ORDER, columns=AA_ORDER)
        print(f"substitution matrix: {len(matrices)} species")
    aligned.index.name = "aa1"
    aligned.columns.name = "aa2"
    return aligned


def plot_vertebrate_amino_acids(amino_acid_mean: pd.Series, observed: pd.DataFrame, average: str) -> None:
    bars = amino_acid_mean.sort_values(ascending=False).rename_axis("aa").reset_index(name="freq")
    plt.figure(figsize=(14, 4))
    ax = sns.barplot(bars, x="aa", y="freq", color="lightgray", edgecolor=".55", label="Equilibrium")
    sns.boxplot(
        data=observed,
        x="aa",
        y="freq",
        order=bars["aa"],
        showfliers=True,
        fliersize=1,
        ax=ax,
        color="lightgreen",
        label="Observed",
    )
    ax.set_title("Vertebrates Cytb")
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.tight_layout()
    plt.savefig(output_path("aa_eq_freqs_vert_cytb.pdf", average))
    plt.close()


def plot_class_heatmap(amino_acid_classes: pd.DataFrame, average: str) -> None:
    means = amino_acid_classes * 100
    column_order = means.mean(axis=0).sort_values(ascending=False).index
    plt.figure(figsize=(16, 4))
    sns.heatmap(
        means.loc[:, column_order],
        cmap="coolwarm",
        annot=True,
        cbar_kws={"label": "Equilibrium frequency (%)"},
    )
    plt.ylabel("")
    plt.xlabel("")
    plt.tight_layout()
    plt.savefig(output_path("aa_eq_freqs_cls_heatmap.pdf", average), bbox_inches="tight")
    plt.close()


def plot_class_nucleotides(nucleotide_classes: pd.DataFrame, average: str) -> None:
    long = nucleotide_classes.reset_index(names="Class").melt(id_vars="Class", var_name="nucl", value_name="freq")
    (
        so.Plot(long, y="Class", x="freq", color="nucl")
        .add(so.Bar(), so.Stack())
        .scale(
            y=so.Nominal(order=CLASS_PLOT_ORDER),
            color=so.Nominal(order=list(nucleotide_classes.columns)),
        )
        .layout(size=(6, 4))
        .save(output_path("nuc_eq_freqs_cls.pdf", average))
    )
    plt.close()


def plot_class_amino_acids(amino_acid_classes: pd.DataFrame, observed: pd.DataFrame, average: str) -> None:
    amino_acid_order = amino_acid_classes.mean(axis=0).sort_values(ascending=False).index
    bars = amino_acid_classes.reset_index(names="Class").melt(id_vars="Class", var_name="aa", value_name="freq")
    panel_order = list(reversed(CLASS_PLOT_ORDER))

    fig, axes = plt.subplots(len(panel_order), 1, figsize=(10, 3.5 * len(panel_order)), sharey=True)
    for axis, class_name in zip(axes, panel_order):
        sns.barplot(
            data=bars[bars["Class"] == class_name],
            x="aa",
            y="freq",
            order=amino_acid_order,
            color="lightgray",
            edgecolor=".55",
            ax=axis,
        )
        sns.boxplot(
            data=observed[observed["Class"] == class_name],
            x="aa",
            y="freq",
            order=amino_acid_order,
            showfliers=True,
            fliersize=1,
            ax=axis,
            color="lightgreen",
            label="Observed",
        )
        axis.set_title(class_name)
        axis.set_xlabel("")
    fig.tight_layout()
    fig.savefig(output_path("aa_eq_obs_freqs_5cls.pdf", average), bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot Cytb equilibria averaged across species.")
    parser.add_argument(
        "--average",
        choices=(AVERAGE_EQUILIBRIUM, AVERAGE_SPECTRA),
        default=AVERAGE_EQUILIBRIUM,
        help=(
            "equilibrium: mean of per-species equilibrium frequencies (default). "
            "spectra: mean of mutation spectra, then the equilibrium of that mean spectrum."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    FIGURE_DIR.mkdir(exist_ok=True)
    amino_acids = pd.read_csv(DATA_DIR / "aa_freqs.csv")
    nucleotides = pd.read_csv(DATA_DIR / "nuc_freqs.csv")
    observed = observations(amino_acids, state_columns(amino_acids), "aa")
    spectra = cytb_spectra()
    if args.average == AVERAGE_SPECTRA:
        nucleotide_classes, amino_acid_classes, _, amino_acid_vertebrate = class_means_from_spectra(spectra)
    else:
        nucleotide_classes, amino_acid_classes, _, amino_acid_vertebrate = class_means_from_equilibria(
            nucleotides, amino_acids
        )
    plot_vertebrate_amino_acids(amino_acid_vertebrate, observed, args.average)
    plot_class_heatmap(amino_acid_classes, args.average)
    plot_class_nucleotides(nucleotide_classes, args.average)
    plot_class_amino_acids(amino_acid_classes, observed, args.average)
    plot_exp_heatmap(
        substitution_matrix(spectra, args.average),
        output_path("exp_aa_subst_vert_matrix.pdf", args.average),
        show=False,
    )


if __name__ == "__main__":
    main()
