# Open questions

These items can change a scientific conclusion. The code changes below are in place. The committed CSVs and PDFs were not regenerated.

## Fixed in code, results not regenerated

- Whole-genome PCAWG spectra in `apply_neutral_model_nuclear_somatic.py` were divided by the exome nucleotide frequencies. They now use the genomic frequencies. The saved `nuclear_somatic_fit_metrics.csv` still reflects the old normalization.
- `MutSpecVertebrates12.csv.gz` contains CO1, CO3, Cytb, and ND2, not ND1. Mammal ND1 rows were previously scored with the ND2 spectrum. The batch script now skips `mam_nd1_gtr` and `mam_nd1_unr` until an ND1 spectrum exists. The committed `datasets_fit_metrics.csv` still contains those two rows.
- `calc_metrics` returns `pearson_corr_clr`, the Pearson correlation of centered log-ratios. Callers that read `pearson_corr` from a fresh metrics dict now use that key. Saved tables still have the old columns (`pearson_corr`, `slope`, `ks_stat`, and others).
- Negative-binomial summaries are taken from the GLM result whose model name is `Negative Binomial`. If that fit fails, the `nb_*` columns are omitted instead of reading whichever model finished second.
- Expected substitution labels (`sbs`) are joined with `|`. They are no longer produced by summing a string column.
- `fit_glm` drops rows with `rate_exp <= 0` before taking the logarithm.
- Human megatree mutations are read from `data/external/raw_human_megatree.tsv`. The `.csv` path did not match the file in the tree.

## Left unchanged on purpose

- `log_likelihood` in `calc_metrics` sums Bernoulli terms on compositional frequencies. It is not a multinomial likelihood, so it should not be compared across datasets as one.
- In `human_germline.ipynb`, each row's `MutRate` divides by the table-wide `ObsNum` and `ExpNum` totals, not by a gene-specific denominator. The same notebook keeps only `StrandRef == 1`, which drops L-strand coding substitutions.
- PCAWG missense extraction takes amino-acid change, codon change, and strand from fixed VEP `INFO` fields 15, 16, and 19, then reads the amino acids with `str[0]` and `str[-1]`. That matches one-letter codes such as `A123V` only when the field is already a one-letter pair. Confirm the field order against one annotated VCF before the next refit.
- `analyses/01_equilibrium/calculate_exp_aa_subst_matrix.py` indexes mutation matrices as `M[destination, source]` and stores outflow on the diagonal. `aminoacidshift.substitutions` uses the other off-diagonal index and a row-stochastic codon matrix. The equilibrium notebook still imports the local module. Compare the two on `data/vert_aa_eq_freqs.csv` before deleting the duplicate.
