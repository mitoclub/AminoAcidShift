"""Write Plotly pages from the committed tables. Does not refit models."""

from pathlib import Path

import pandas as pd
import plotly.express as px

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent / "interactive"


def write(fig, name: str) -> None:
    fig.write_html(OUT / name, include_plotlyjs="cdn", full_html=True)


def main() -> None:
    OUT.mkdir(exist_ok=True)

    expected = pd.read_csv(ROOT / "analyses/01_equilibrium/data/exp_aa_subst_vert.csv")
    expected = expected[expected["aa1"] != expected["aa2"]]
    matrix = expected.pivot(index="aa1", columns="aa2", values="rate").fillna(0.0)
    heatmap = px.imshow(
        matrix,
        color_continuous_scale="Reds",
        title="Expected amino-acid substitution rates, vertebrate Cytb",
        aspect="equal",
    )
    heatmap.update_layout(xaxis_title="To", yaxis_title="From")
    write(heatmap, "expected_aa.html")

    datasets = pd.read_csv(ROOT / "analyses/04_datasets/datasets_fit_metrics.csv")
    fit = px.scatter(
        datasets,
        x="mut_count",
        y="pearson_corr_clr",
        text="dataset",
        hover_data=["r2", "spearman_corr", "mut_type_count"],
        title="Neutral-model fit across mtDNA datasets",
    )
    fit.update_traces(textposition="top center")
    write(fit, "datasets_fit.html")

    megatree = pd.read_csv(ROOT / "analyses/03_human_germline/megatree_fit_metrics.csv")
    genes = px.bar(
        megatree,
        x="gene",
        y="r2",
        hover_data=["mut_count", "pearson_corr_clr"],
        title="Human germline fit by gene",
    )
    write(genes, "megatree_r2.html")

    pcawg = pd.read_csv(ROOT / "analyses/05_pcawg/data/nuclear_somatic_fit_metrics.csv")
    neutral = pcawg.loc[pcawg["model"] == "neutral"]
    samples = px.scatter(
        neutral,
        x="mut_count",
        y="r2",
        color="disease",
        hover_data=["sample", "pearson_corr"],
        title="PCAWG neutral-model fit, saved per sample",
    )
    write(samples, "pcawg_fit.html")


if __name__ == "__main__":
    main()
