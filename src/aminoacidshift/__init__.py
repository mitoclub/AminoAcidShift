"""Public API for the amino-acid substitution model."""

from aminoacidshift.constants import (
    alphabet,
    amino_acid_codes,
    color_mapping12,
    grantham_long,
    transitions,
    transversions,
)
from aminoacidshift.metrics import (
    calc_accuracy,
    calc_f1,
    calc_metrics,
    calc_slope,
    fit_glm,
    transform_clr,
    weighted_average_percentage_error,
)
from aminoacidshift.plotting import (
    plot_aa_eq_freqs,
    plot_exp_heatmap,
    plot_obs_vs_exp,
    plot_subst_freqs,
)
from aminoacidshift.site_rates import (
    calculate_truncated_mean,
    categorize_site_rates_robust_plus_invariant,
)
from aminoacidshift.substitutions import (
    cdn_spectrum_to_matrix,
    collect_possible_changes,
    get_equilibrium_freqs,
    get_equilibrium_freqs_rnd,
    get_equilibrium_probabilities,
    get_random_spectra,
    get_random_spectrum,
    nuc_spectrum_to_matrix,
    prepare_aa_subst,
    prepare_exp_aa_subst,
    prepare_exp_aa_subst_codons,
    prepare_exp_cdn_subst,
    prepare_rnd_exp_aa_subst,
)

__all__ = [
    "alphabet",
    "amino_acid_codes",
    "calc_accuracy",
    "calc_f1",
    "calc_metrics",
    "calc_slope",
    "calculate_truncated_mean",
    "categorize_site_rates_robust_plus_invariant",
    "cdn_spectrum_to_matrix",
    "collect_possible_changes",
    "color_mapping12",
    "fit_glm",
    "get_equilibrium_freqs",
    "get_equilibrium_freqs_rnd",
    "get_equilibrium_probabilities",
    "get_random_spectra",
    "get_random_spectrum",
    "grantham_long",
    "nuc_spectrum_to_matrix",
    "plot_aa_eq_freqs",
    "plot_exp_heatmap",
    "plot_obs_vs_exp",
    "plot_subst_freqs",
    "prepare_aa_subst",
    "prepare_exp_aa_subst",
    "prepare_exp_aa_subst_codons",
    "prepare_exp_cdn_subst",
    "prepare_rnd_exp_aa_subst",
    "transform_clr",
    "transitions",
    "transversions",
    "weighted_average_percentage_error",
]
