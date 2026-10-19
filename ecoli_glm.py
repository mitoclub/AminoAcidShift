"""Python rewrite of the E. coli negative-binomial analysis.

This version uses the saved ``dataEcoli.txt`` table directly. That table
already contains the engineered predictors and the original R negative-binomial
predictions in ``predevents``. The script refits the models in Python and
compares the new predictions against the saved ones.

https://doi.org/10.1073/pnas.2119720119
https://github.com/alejvcano/Mutbias2022
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import pearsonr, spearmanr
from statsmodels.discrete.count_model import ZeroInflatedNegativeBinomialP, ZeroInflatedPoisson
from statsmodels.tools import add_constant

def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def fit_poisson(data: pd.DataFrame):
    exog = add_constant(data[["log_mutAll"]], has_constant="add")
    model = sm.GLM(data["events"], exog, family=sm.families.Poisson(), offset=data["log_frequency"])
    return model.fit()


def fit_negative_binomial(data: pd.DataFrame, predictors: Sequence[str], maxiter: int = 50):
    exog = add_constant(data[list(predictors)], has_constant="add")
    model = sm.NegativeBinomial(data["events"], exog, offset=data["log_frequency"])
    return model.fit(maxiter=maxiter, disp=False)


def fit_zero_inflated_poisson(data: pd.DataFrame):
    exog = add_constant(data[["log_mutAll"]], has_constant="add")
    exog_infl = np.ones((len(data), 1))
    model = ZeroInflatedPoisson(
        data["events"],
        exog,
        exog_infl=exog_infl,
        offset=data["log_frequency"],
        inflation="logit",
    )
    return model.fit(method="bfgs", maxiter=200, disp=False)


def fit_zero_inflated_nb(data: pd.DataFrame, predictors: Sequence[str]):
    exog = add_constant(data[list(predictors)], has_constant="add")
    exog_infl = np.ones((len(data), 1))
    model = ZeroInflatedNegativeBinomialP(
        data["events"],
        exog,
        exog_infl=exog_infl,
        offset=data["log_frequency"],
        inflation="logit",
    )
    return model.fit(method="bfgs", maxiter=200, disp=False)


def predicted_response(result, data: pd.DataFrame, predictors: Sequence[str]) -> np.ndarray:
    exog = add_constant(data[list(predictors)], has_constant="add")
    return np.asarray(result.predict(exog=exog, offset=data["log_frequency"]))


def compare_aic(*results) -> pd.Series:
    return pd.Series({f"model_{i+1}": res.aic for i, res in enumerate(results)})

def run_resampling(data: pd.DataFrame, n_iter: int, seed: int = 12345) -> np.ndarray:
    rng = np.random.default_rng(seed)
    corrdist: List[float] = []
    predictors = ["log_mutAll"]

    while len(corrdist) < n_iter:
        resample_data = data.copy()
        resample_data["events"] = rng.permutation(resample_data["events"].to_numpy())
        try:
            temp = fit_negative_binomial(resample_data, predictors, maxiter=50)
            fit = predicted_response(temp, resample_data, predictors)
            corrdist.append(float(pearsonr(fit, resample_data["events"])[0]))
        except Exception:
            continue

    return np.asarray(corrdist)


def load_saved_table(base_dir: Path, table_path: Path) -> pd.DataFrame:
    data = read_table(table_path)

    required = ["events", "frequency", "log_frequency", "log_avgMutRate", "log_pathsCount", "log_mutAll"]
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise ValueError(
            "dataEcoli.txt is missing required columns: " + ", ".join(missing)
        )

    return data


def print_model_summary(name: str, result) -> None:
    print(f"\n=== {name} ===")
    try:
        print(result.summary())
    except Exception:
        print(result)


def main() -> None:
    parser = argparse.ArgumentParser(description="Python rewrite of ecoli_glm.Rmd")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory containing dataEcoli.txt",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Path to dataEcoli.txt",
    )
    parser.add_argument(
        "--resamples",
        type=int,
        default=0,
        help="Number of permutations for the resampling correlation test",
    )
    parser.add_argument(
        "--skip-resampling",
        action="store_true",
        help="Skip the permutation test section",
    )
    args = parser.parse_args()

    table_path = args.input or (args.base_dir / "dataEcoli.txt")
    if not table_path.exists():
        raise FileNotFoundError(f"Missing input table: {table_path}")

    data = load_saved_table(args.base_dir, table_path)

    # print(data.head(20).to_string(index=False))

    poisson = fit_poisson(data)
    nb_1a = fit_negative_binomial(data, ["log_mutAll"], maxiter=50)
    nb_1b = fit_negative_binomial(data, ["log_avgMutRate", "log_pathsCount"], maxiter=50)
    nb_paths = fit_negative_binomial(data, ["log_pathsCount"], maxiter=50)
    nb_avg = fit_negative_binomial(data, ["log_avgMutRate"], maxiter=50)

    print_model_summary("Poisson", poisson)
    print_model_summary("Negative binomial model 4a", nb_1a)

    print(nb_1a._results.summary2().tables[1].to_string())

    print_model_summary("Negative binomial model 4b", nb_1b)
    print_model_summary("Negative binomial paths only", nb_paths)
    print_model_summary("Negative binomial avg rate only", nb_avg)

    try:
        zip_poisson = fit_zero_inflated_poisson(data)
        print_model_summary("Zero-inflated Poisson", zip_poisson)
        print(compare_aic(poisson, zip_poisson).to_string())
    except Exception as exc:
        print(f"Zero-inflated Poisson failed: {exc}")

    try:
        zinb_1 = fit_zero_inflated_nb(data, ["log_mutAll"])
        zinb_2 = fit_zero_inflated_nb(data, ["log_avgMutRate", "log_pathsCount"])
        print_model_summary("Zero-inflated NB model 4a", zinb_1)
        print_model_summary("Zero-inflated NB model 4b", zinb_2)
        print(compare_aic(nb_1a, zinb_1, zinb_2).to_string())
    except Exception as exc:
        print(f"Zero-inflated negative binomial failed: {exc}")

    data["predevents_python"] = predicted_response(nb_1a, data, ["log_mutAll"])
    pearson_corr = pearsonr(data["predevents_python"], data["events"])
    spearman_corr = spearmanr(data["predevents_python"], data["events"])
    print(f"Pearson correlation between predicted and real events: {pearson_corr.statistic:.6f}")
    print(f"Pearson p-value: {pearson_corr.pvalue:.6g}")
    print(f"Spearman correlation between predicted and real events: {spearman_corr.correlation:.6f}")
    print(f"Spearman p-value: {spearman_corr.pvalue:.6g}")

    if "predevents" in data.columns:
        saved_pearson = pearsonr(data["predevents_python"], data["predevents"])
        saved_spearman = spearmanr(data["predevents_python"], data["predevents"])
        print()
        print(f"Pearson correlation between Python and saved predictions: {saved_pearson.statistic:.6f}")
        print(f"Pearson p-value vs saved predictions: {saved_pearson.pvalue:.6g}")
        print(f"Spearman correlation between Python and saved predictions: {saved_spearman.correlation:.6f}")
        print(f"Spearman p-value vs saved predictions: {saved_spearman.pvalue:.6g}")

    data["log_events"] = data["events"].astype(float) + 0.5
    data["log_offset"] = data["log_events"] - data["log_frequency"]

    data["predevents_freq_only"] = predicted_response(
        fit_negative_binomial(data, ["log_frequency"], maxiter=500), data, ["log_frequency"]
    )

    spearman_pair = spearmanr(data["log_mutAll"], data["log_offset"])
    print(f"Spearman(log_mutAll, log_offset) = {spearman_pair.correlation:.6f}")

    pearson_fit = pearsonr(data["predevents_python"], data["events"])
    print(f"Pearson(predicted, real) = {pearson_fit.statistic:.6f}")

    if not args.skip_resampling and args.resamples > 0:
        corrdist = run_resampling(data, args.resamples)
        out_path = args.base_dir / "corrdistEcoli.txt"
        pd.Series(corrdist).to_csv(out_path, sep="\t", index=False, header=False)
        print(f"Saved {len(corrdist)} resampling correlations to corrdistEcoli.txt")


if __name__ == "__main__":
    main()