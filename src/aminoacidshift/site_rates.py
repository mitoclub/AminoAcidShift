import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma


def calculate_truncated_mean(alpha, scale, lower, upper):
    """Conditional expectation E[X | lower < X < upper] for a Gamma distribution."""
    prob_mass = gamma.cdf(upper, a=alpha, scale=scale) - gamma.cdf(lower, a=alpha, scale=scale)
    if prob_mass == 0:
        return 0.0

    num_upper = gamma.cdf(upper, a=alpha + 1, scale=scale)
    num_lower = gamma.cdf(lower, a=alpha + 1, scale=scale)
    return (alpha * scale) * (num_upper - num_lower) / prob_mass


def categorize_site_rates_robust_plus_invariant(
    mut_counts, n_categories=6, hotspot_percentile=99.0, plot=False, verbose=True
):
    """
    Categorize sites into invariant sites, gamma-distributed rates, and hotspots.

    Category 0 is sites with zero observed mutations. Categories 1 to N-1 are
    gamma bins. Category N is the hotspot tail.
    """
    rates_observed = np.array(mut_counts, dtype=float)
    global_mean = np.mean(rates_observed)
    if global_mean == 0:
        site_categories = np.zeros(len(rates_observed), dtype=int)
        category_rates = np.zeros(n_categories)
        return site_categories, category_rates, np.nan

    rates_observed = rates_observed / global_mean
    mask_zero = rates_observed == 0
    rates_zero = rates_observed[mask_zero]
    rates_nonzero = rates_observed[~mask_zero]

    cutoff_value = np.percentile(rates_nonzero, hotspot_percentile)
    mask_normal = rates_nonzero <= cutoff_value
    mask_hotspot = rates_nonzero > cutoff_value
    rates_normal = rates_nonzero[mask_normal]
    rates_hotspot = rates_nonzero[mask_hotspot]

    if verbose:
        print(
            f"Summary: {len(rates_zero)} zero sites, "
            f"{len(rates_normal)} gamma sites, {len(rates_hotspot)} hotspots."
        )

    alpha, _loc, scale = gamma.fit(rates_normal, floc=0)
    n_gamma_cats = n_categories - 2
    quantiles = np.linspace(0, 1, n_gamma_cats + 1)
    gamma_thresholds = gamma.ppf(quantiles, a=alpha, scale=scale)
    gamma_thresholds[-1] = cutoff_value

    category_rates = [0.0]
    for i in range(n_gamma_cats):
        lower, upper = gamma_thresholds[i], gamma_thresholds[i + 1]
        category_rates.append(calculate_truncated_mean(alpha, scale, lower, upper))

    if len(rates_hotspot) > 0:
        hotspot_rate = np.mean(rates_hotspot)
    else:
        hotspot_rate = cutoff_value * 1.5
    category_rates.append(hotspot_rate)
    category_rates = np.array(category_rates)

    site_categories = np.zeros(len(rates_observed), dtype=int)
    nonzero_indices = np.where(~mask_zero)[0]
    digitization_bins = np.concatenate([gamma_thresholds, [np.inf]])
    gamma_hotspot_cats = np.digitize(rates_nonzero, digitization_bins[1:]) + 1
    site_categories[nonzero_indices] = gamma_hotspot_cats

    if plot:
        plt.figure(figsize=(12, 6))
        max_plot_x = np.percentile(rates_observed, 99.5) * 1.2
        bins = np.linspace(0, max_plot_x, 100)
        plt.hist(rates_observed, bins=bins, density=True, alpha=0.5, color="gray", label="Observed Rates")
        x = np.linspace(0, max_plot_x, 500)
        prop_normal = len(rates_normal) / len(rates_observed)
        plt.plot(
            x,
            gamma.pdf(x, a=alpha, scale=scale) * prop_normal,
            "r-",
            lw=2,
            label="Gamma Fit (Normal Data Only)",
        )
        for threshold in gamma_thresholds[1:-1]:
            plt.axvline(threshold, color="blue", linestyle="--", alpha=0.6)
        plt.axvline(cutoff_value, color="red", linestyle="-", lw=2, label="Hotspot Cutoff")
        for rate in category_rates:
            plt.scatter(rate, 0, color="black", zorder=5, s=50, marker="x")
        plt.xlabel("Relative Site Rate")
        plt.ylabel("Density")
        plt.title(f"Site Categorization: {n_gamma_cats} Gamma + 1 Hotspot")
        plt.legend()
        plt.show()

    return site_categories, category_rates, alpha
