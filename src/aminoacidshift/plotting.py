import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import r2_score

from aminoacidshift.metrics import calc_slope


def plot_exp_heatmap(exp_aa_subst_matrix: pd.DataFrame, save_path: str, show=True, annot=False):
    """Plot an expected amino-acid substitution matrix."""
    freqs_from = (exp_aa_subst_matrix.sum(1) - exp_aa_subst_matrix.values.diagonal()).copy()
    freqs_to = (exp_aa_subst_matrix.sum(0) - exp_aa_subst_matrix.values.diagonal()).copy()
    max_flow_value = max(freqs_to.max(), freqs_from.max()) + 0.1

    fig, axs = plt.subplots(
        2, 3, figsize=(11, 10), width_ratios=[0.1, 1, 0.05], height_ratios=[1, 0.1]
    )
    sns.heatmap(
        exp_aa_subst_matrix,
        annot=annot,
        fmt=".2f",
        ax=axs[0, 1],
        cbar_ax=axs[0, 2],
        cmap="light:r",
        cbar_kws={"label": "Substitution rate"},
        mask=exp_aa_subst_matrix == 0,
    )
    axs[0, 1].set_ylabel("")
    axs[0, 1].set_xlabel("")
    axs[0, 1].set_xticks([])
    axs[0, 1].set_yticks([])

    sns.barplot(freqs_from.reset_index(), y="aa1", x=0, ax=axs[0, 0], color="lightgray", edgecolor="black")
    axs[0, 0].set_ylabel("From", fontsize=14)
    axs[0, 0].set_xlabel("")
    axs[0, 0].spines["top"].set_visible(False)
    axs[0, 0].spines["right"].set_visible(False)
    axs[0, 0].spines["bottom"].set_visible(False)
    axs[0, 0].spines["left"].set_visible(False)
    axs[0, 0].invert_xaxis()
    axs[1, 1].set_ylim(0, max_flow_value)
    axs[0, 0].set_xticks([])
    axs[0, 0].set_yticklabels(axs[0, 0].get_yticklabels(), fontsize=12)

    sns.barplot(freqs_to.reset_index(), x="aa2", y=0, ax=axs[1, 1], color="lightgray", edgecolor="black")
    axs[1, 1].set_ylabel("")
    axs[1, 1].set_xlabel("To", fontsize=14)
    axs[1, 1].spines["top"].set_visible(False)
    axs[1, 1].spines["right"].set_visible(False)
    axs[1, 1].spines["bottom"].set_visible(False)
    axs[1, 1].spines["left"].set_visible(False)
    axs[1, 1].invert_yaxis()
    axs[1, 1].set_ylim(max_flow_value, 0)
    axs[1, 1].set_yticks([])
    axs[1, 1].set_xticklabels(axs[1, 1].get_xticklabels(), fontsize=12)

    axs[1, 0].remove()
    axs[1, 2].remove()

    plt.tight_layout()
    plt.savefig(save_path)
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_aa_eq_freqs(exp_aa_subst_matrix, save_path: str, show=True, figsize=(12, 3)):
    freqs_from = (exp_aa_subst_matrix.sum(1) - exp_aa_subst_matrix.values.diagonal()).copy()
    freqs_to = (exp_aa_subst_matrix.sum(0) - exp_aa_subst_matrix.values.diagonal()).copy()

    _freqs_to = (freqs_to / freqs_to.sum()).copy()
    _freqs_from = (freqs_from / freqs_from.sum()).copy()
    _freqs_to.index.name = "aa"
    _freqs_from.index.name = "aa"

    plt.figure(figsize=figsize)
    ax = sns.barplot(
        (_freqs_to - _freqs_from).sort_values(ascending=False).reset_index(),
        x="aa",
        y=0,
        color="lightgray",
        edgecolor=".65",
    )

    for bar in ax.patches:
        if bar.get_height() > 0.05:
            bar.set_color("orangered")
            bar.set_edgecolor("0.02")
        elif bar.get_height() > 0.02:
            bar.set_color("lightcoral")
            bar.set_edgecolor("0.05")
        elif bar.get_height() < -0.05:
            bar.set_color("royalblue")
            bar.set_edgecolor("0.02")
        elif bar.get_height() < -0.02:
            bar.set_color("lightblue")
            bar.set_edgecolor("0.45")

    plt.legend(
        [
            plt.Rectangle((0, 0), 1, 1, fc="orangered", edgecolor="0.45"),
            plt.Rectangle((0, 0), 1, 1, fc="lightcoral", edgecolor="0.45"),
            plt.Rectangle((0, 0), 1, 1, fc="lightblue", edgecolor="0.45"),
            plt.Rectangle((0, 0), 1, 1, fc="royalblue", edgecolor="0.45"),
        ],
        ["Top Gainers", "Gainers", "Loosers", "Top Loosers"],
        loc="upper right",
    )

    plt.ylabel("Rates [To - From]", fontsize=13)
    plt.xlabel("")
    plt.xticks(fontsize=12)
    plt.savefig(save_path)
    if show:
        plt.show()
    else:
        plt.close()


def plot_obs_vs_exp(
    aa_subst, save_path=None, ax=None, show=True, text="", text_x=-2.2, text_y=-4.0
):
    aa_subst = aa_subst.copy()
    aa_subst["nobs_freqs_log"] = np.log10(aa_subst["nobs_freqs"])
    aa_subst["nexp_freqs_log"] = np.log10(aa_subst["nexp_freqs"])

    rho, pval = spearmanr(aa_subst["nobs_freqs"], aa_subst["nexp_freqs"])
    k, b = calc_slope(aa_subst["nobs_freqs"], aa_subst["nexp_freqs"])

    if ax is None:
        plt.figure(figsize=(5, 5))
        ax = plt.gca()

    ax = sns.scatterplot(
        aa_subst, ax=ax, x="nobs_freqs_log", y="nexp_freqs_log", color="blue", alpha=0.5
    )
    ax = sns.regplot(
        aa_subst[aa_subst["nobs_freqs"] > 0],
        color="blue",
        scatter=False,
        x="nobs_freqs_log",
        y="nexp_freqs_log",
        ax=ax,
        label=f"y={k:.2f}x+{b:.1g}",
    )

    ticks = np.linspace(-4, -1, 4)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    formatter = lambda x, pos: "10$^{" + str(int(x)) + "}$"
    ax.get_xaxis().set_major_formatter(formatter)
    ax.get_yaxis().set_major_formatter(formatter)

    ax.text(text_x, text_y, f"{text}\nSpearman: {rho:.2f}\np-value: {pval:.2g}", fontsize=10)
    ax.plot([-4, -1], [-4, -1], color="black", linestyle="--", label="y=x")
    ax.set_xlabel("Observed AA frequencies")
    ax.set_ylabel("Expected AA frequencies")
    ax.legend(loc="upper left", fontsize=10)
    if save_path:
        plt.savefig(save_path)
    if show:
        plt.show()


def plot_subst_freqs(aa_subst, title=""):
    aa_subst = aa_subst.copy()
    aa_subst["nobs_freqs_log"] = np.log10(aa_subst["nobs_freqs"])
    aa_subst["nexp_freqs_log"] = np.log10(aa_subst["nexp_freqs"])

    y_true, y_pred = aa_subst["nobs_freqs"], aa_subst["nexp_freqs"]
    cor_res = pearsonr(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    print(r2)
    print(f"Pearson correlation: {cor_res.correlation:.3f} (p-value: {cor_res.pvalue:.3g})")

    plt.figure(figsize=(8, 8))
    ax = sns.regplot(
        aa_subst[aa_subst["nobs_freqs"] > 0],
        color="blue",
        scatter_kws={"alpha": 0.5, "s": 50},
        y="nobs_freqs_log",
        x="nexp_freqs_log",
    )

    ticks = np.linspace(-4, -1, 4)
    ticks_minor = np.log10(
        np.concat(
            [
                np.linspace(10**-4, 10**-3, 10),
                np.linspace(10**-3, 10**-2, 10)[1:],
                np.linspace(10**-2, 10**-1, 10)[1:],
            ]
        )
    )
    ax.set_xticks(ticks, ticks, size=14)
    ax.set_yticks(ticks, ticks, size=14)
    ax.set_xticks(ticks_minor, minor=True)
    ax.set_yticks(ticks_minor, minor=True)
    formatter = lambda x, pos: "10$^{" + str(int(x)) + "}$"
    ax.get_xaxis().set_major_formatter(formatter)
    ax.get_yaxis().set_major_formatter(formatter)

    plt.plot([-4, -1], [-4, -1], color="black", linestyle="--")
    plt.ylabel("Observed AA substitution frequencies", fontsize=16)
    plt.xlabel("Predicted AA substitution frequencies", fontsize=16)
    plt.title(title, fontsize=16)
    plt.legend(title=f"$R^2$ = {r2:.2f}", title_fontsize=16)
    return ax
