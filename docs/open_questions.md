# Open questions

These items can change a scientific conclusion. The code changes below are in place. Megatree, cross-dataset, and equilibrium tables were regenerated with the current metrics. PCAWG tables were not, because the missense table and tumour-spectra catalogs are not in the tree.

## Fixed in code

- Whole-genome PCAWG spectra in `apply_neutral_model_nuclear_somatic.py` were divided by the exome nucleotide frequencies. They now use the genomic frequencies. `nuclear_somatic_fit_metrics.csv` still reflects the old normalization.
- `MutSpecVertebrates12.csv.gz` contains CO1, CO3, Cytb, and ND2, not ND1. Mammal ND1 rows were previously scored with the ND2 spectrum. The batch script now skips a run when the observed table or the 12-class spectrum is missing. The regenerated `datasets_fit_metrics.csv` has no mammal rows and no ND1 rows.
- `calc_metrics` returns `pearson_corr_clr`, the Pearson correlation of centered log-ratios. Regenerated megatree and cross-dataset tables use that column. The PCAWG table still has `pearson_corr`, `slope`, and `ks_stat`.
- Negative-binomial summaries are taken from the GLM result whose model name is `Negative Binomial`. If that fit fails, the `nb_*` columns are omitted instead of reading whichever model finished second.
- Expected substitution labels (`sbs`) are joined with `|`. They are no longer produced by summing a string column.
- `fit_glm` drops rows with `rate_exp <= 0` before taking the logarithm.
- Human megatree mutations are read from `data/external/raw_human_megatree.tsv`. The `.csv` path did not match the file in the tree.

## Left unchanged on purpose

- `log_likelihood` in `calc_metrics` sums Bernoulli terms on compositional frequencies. It is not a multinomial likelihood, so it should not be compared across datasets as one.
- In `human_germline.ipynb`, each row's `MutRate` divides by the table-wide `ObsNum` and `ExpNum` totals, not by a gene-specific denominator. The same notebook keeps only `StrandRef == 1`, which drops L-strand coding substitutions.
- PCAWG missense extraction takes amino-acid change, codon change, and strand from fixed VEP `INFO` fields 15, 16, and 19, then reads the amino acids with `str[0]` and `str[-1]`. That matches one-letter codes such as `A123V` only when the field is already a one-letter pair. Confirm the field order against one annotated VCF before the next refit.
- `analyses/01_equilibrium/calculate_exp_aa_subst_matrix.py` indexes mutation matrices as `M[destination, source]` and stores outflow on the diagonal. `aminoacidshift.substitutions` uses the other off-diagonal index and a row-stochastic codon matrix. The equilibrium notebook still imports the local module. Compare the two on `data/vert_aa_eq_freqs.csv` before deleting the duplicate.
