import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_percentage_error, mean_squared_error, r2_score

from aminoacidshift.glm import (
    fit_negative_binomial,
    fit_poisson,
    fit_zero_inflated_nb,
    fit_zero_inflated_poisson,
    predicted_response,
)


def transform_clr(v, func_log=np.log):
    """Centered log-ratio transformation for compositional data."""
    v = np.asarray(v)
    if np.any(v <= 0) and func_log == np.log:
        raise ValueError("All elements of input vector must be positive for CLR transformation.")
    log_v = func_log(v)
    mean_log_v = np.mean(log_v)
    return log_v - mean_log_v


def calc_accuracy(y_true, y_pred):
    """Legacy accuracy: sum of per-class minima divided by the observed total."""
    true_positive = np.min([y_true, y_pred], axis=0)
    return true_positive.sum() / y_true.sum()


def calc_f1(y_true, y_pred):
    """Legacy macro and weighted F1 on compositional counts."""
    f1_weighted = 0.0
    f1_macro = 0.0
    total = y_true.sum()
    n = len(y_true)

    for observed, expected in zip(y_true, y_pred):
        if observed < expected:
            tp = observed
            fp = expected - observed
            fn = 0
        elif observed > expected:
            tp = expected
            fp = 0
            fn = observed - expected
        else:
            tp = observed
            fp = fn = 0

        f1 = (2 * tp) / (2 * tp + fp + fn)
        f1_macro += f1 / n
        f1_weighted += f1 * (observed / total)

    return f1_macro, f1_weighted


def calc_slope(y_true, y_pred):
    """Slope and intercept of y_pred ~ y_true."""
    try:
        slope, intercept = np.polyfit(y_true, y_pred, 1)
    except Exception:
        slope, intercept = np.nan, np.nan
    return slope, intercept


def weighted_average_percentage_error(y_true, y_pred):
    return mean_absolute_percentage_error(y_true, y_pred, sample_weight=y_true)


def _try_fit(name, fitter):
    try:
        return fitter()
    except Exception as exc:
        print(f"{name} failed: {exc}")
        return None


def fit_glm(aa_subst: pd.DataFrame):
    """Fit count GLMs of observed substitutions on the expected mutational rate.

    Rows with a non-positive expected rate are dropped before the logarithm.
    """
    usable = aa_subst.loc[aa_subst["rate_exp"] > 0].copy()
    data = usable[["aa1", "aa2", "sbs"]].copy()
    data["events"] = usable["nobs"].astype(int)
    data["log_mutAll"] = np.log(usable["rate_exp"].to_numpy())
    data["log_frequency"] = np.log(usable["ref_aa1_freq"].to_numpy())

    poisson = _try_fit(
        "Poisson",
        lambda: fit_poisson(data),
    )
    nb = _try_fit(
        "Negative Binomial",
        lambda: fit_negative_binomial(data, ["log_mutAll"], maxiter=50),
    )
    zip_poisson = _try_fit("Zero-inflated Poisson", lambda: fit_zero_inflated_poisson(data))
    zinb = _try_fit(
        "Zero-inflated Negative Binomial",
        lambda: fit_zero_inflated_nb(data, ["log_mutAll"]),
    )

    results = []
    for model, name in zip(
        [poisson, nb, zip_poisson, zinb],
        [
            "Poisson",
            "Negative Binomial",
            "Zero-inflated Poisson",
            "Zero-inflated Negative Binomial",
        ],
    ):
        if model is None:
            continue
        try:
            predevents = predicted_response(model, data, ["log_mutAll"])
            pearson_corr, pearson_p = pearsonr(data["events"], predevents)
            spearman_corr, spearman_p = spearmanr(data["events"], predevents)
            r2 = r2_score(data["events"], predevents)
        except Exception as exc:
            print(f"Failed to calculate metrics for {name}: {exc}")
            continue
        results.append(
            {
                "model": name,
                "r2": r2,
                "pearson_corr": pearson_corr,
                "pearson_p": pearson_p,
                "spearman_corr": spearman_corr,
                "spearman_p": spearman_p,
                "log_likelihood": model._results.llf,
                "beta": model.params["log_mutAll"],
                "beta_se": model.bse["log_mutAll"],
            }
        )

    if nb is not None:
        try:
            data["predevents"] = predicted_response(nb, data, ["log_mutAll"])
        except Exception as exc:
            print(f"Negative binomial prediction failed: {exc}")
    return results, data


def calc_metrics(aa_subst: pd.DataFrame):
    aa_subst = aa_subst.dropna().query("rate_exp > 0")
    nobs_clr = transform_clr(aa_subst["nobs_freqs"], np.log1p)
    nexp_clr = transform_clr(aa_subst["nexp_freqs"], np.log1p)

    # Sum of Bernoulli terms on compositional frequencies. This is not a
    # multinomial likelihood; see docs/open_questions.md.
    log_likelihood = np.sum(
        aa_subst["nobs_freqs"] * np.log(aa_subst["nexp_freqs"])
        + (1 - aa_subst["nobs_freqs"]) * np.log(1 - aa_subst["nexp_freqs"])
    )

    try:
        rmse_clr = mean_squared_error(nobs_clr, nexp_clr) ** 0.5
    except Exception:
        rmse_clr = np.nan

    mask = aa_subst["nobs"] > 0
    mape = mean_absolute_percentage_error(aa_subst["nobs_freqs"][mask], aa_subst["nexp_freqs"][mask])
    wape = weighted_average_percentage_error(aa_subst["nobs_freqs"][mask], aa_subst["nexp_freqs"][mask])

    spearman_corr, spearman_p = spearmanr(aa_subst["nobs_freqs"], aa_subst["nexp_freqs"])
    pearson_corr_clr, pearson_p_clr = pearsonr(nobs_clr, nexp_clr)

    mut_count = np.sum(aa_subst["nobs"])
    r2 = r2_score(aa_subst["nobs_freqs"], aa_subst["nexp_freqs"])
    r2_clr = r2_score(nobs_clr, nexp_clr)

    corr_chem_vs_rel_freq = pearsonr(
        aa_subst["obs_relative_freq"],
        aa_subst["grantham_distance"],
    )

    metrics = {
        "r2": r2,
        "r2_clr": r2_clr,
        "mape": mape,
        "wape": wape,
        "spearman_corr": spearman_corr,
        "spearman_p": spearman_p,
        "pearson_corr_clr": pearson_corr_clr,
        "pearson_p_clr": pearson_p_clr,
        "rmse_clr": rmse_clr,
        "log_likelihood": log_likelihood,
        "mut_count": mut_count,
        "mut_type_count": aa_subst.nobs.ne(0).sum(),
        "corr_chem_vs_rel_freq": corr_chem_vs_rel_freq.correlation,
        "corr_chem_vs_rel_freq_pval": corr_chem_vs_rel_freq.pvalue,
    }
    try:
        res_glm, _ = fit_glm(aa_subst)
    except Exception as exc:
        print(f"GLM fit failed: {exc}")
        res_glm = []
    nb_metrics = next((row for row in res_glm if row.get("model") == "Negative Binomial"), None)
    if nb_metrics is not None:
        metrics.update({f"nb_{key}": value for key, value in nb_metrics.items()})

    return metrics
