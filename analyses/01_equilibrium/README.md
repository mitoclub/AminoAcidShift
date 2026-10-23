# Equilibrium frequencies

Question: how far is the observed nucleotide, codon, and amino-acid content of each species and gene from the equilibrium implied by its own mitochondrial mutation spectrum?

Inputs

1. [`MutSpecVertebrates12.csv.gz`](./data/raw/MutSpecVertebrates12.csv.gz): 12-class spectra for CO1, CO3, Cytb, and ND2.
2. Coding sequences in `data/raw/output/{Gene}/{Gene}__{Species}/seqs_unique.fasta`. The outgroup record (`OUTGRP`) is excluded. Frequencies are the mean across the remaining sequences.

Run from this directory:

```bash
uv run python compare_content_to_equilibrium.py
```

The script uses vertebrate mitochondrial code 2 and the equilibrium functions in `aminoacidshift`. It writes nine tables in `data/`:

- `nuc_eq_freqs.csv`, `nuc_obs_freqs.csv`, `nuc_distances.csv`
- `cdn_eq_freqs.csv`, `cdn_obs_freqs.csv`, `cdn_distances.csv`
- `aa_eq_freqs.csv`, `aa_obs_freqs.csv`, `aa_distances.csv`

Nucleotide distances are reported for the whole CDS, each codon position, and four-fold third positions. Distance distributions are in `figures/nuc_eq_distance_distributions.pdf`, `figures/cdn_eq_distance_distributions.pdf`, and `figures/aa_eq_distance_distributions.pdf`.

`prepare_data_for_aa_shift.ipynb` and `gene_strand_forecast.ipynb` are the earlier class-level and strand explorations. `calculate_exp_aa_subst_matrix.py` keeps its own matrix orientation; see [`docs/open_questions.md`](../../docs/open_questions.md).
