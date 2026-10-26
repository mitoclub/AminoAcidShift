# Equilibrium frequencies

Question: how far is the observed nucleotide, codon, and amino-acid content of each species and gene from the equilibrium implied by its own mitochondrial mutation spectrum?

Inputs

1. [`MutSpecVertebrates12.csv.gz`](./data/raw/MutSpecVertebrates12.csv.gz): 12-class spectra. This analysis uses Cytb and ND2.
2. Coding sequences in `data/raw/output/{Gene}/{Gene}__{Species}/seqs_unique.fasta`. The outgroup record (`OUTGRP`) is excluded. Frequencies are the mean across the remaining sequences.

Run from this directory:

```bash
uv run python derive_freqs.py
uv run python compare_content_to_equilibrium.py
```

`derive_freqs.py` uses vertebrate mitochondrial code 2 and the equilibrium functions in `aminoacidshift`. It writes three tables in `data/`. Each has columns `Class`, `Gene`, `Species`, `Type`, then one column per state. `Type` is `obs` or `eq`, so each species has two rows. Missing equilibria are left blank.

- `nuc_freqs.csv` also has `subset` (`all`, `pos1`, `pos2`, `pos3`, `fourfold`) before `Type`, and columns `A`, `C`, `G`, `T`. The equilibrium does not depend on the subset, so it is repeated on each subset.
- `cdn_freqs.csv` has one column per codon.
- `aa_freqs.csv` has one column per amino acid.

`compare_content_to_equilibrium.py` reads these tables and writes `distances.csv`, which is not committed. Distance boxplots are `figures/nuc_eq_distance_distributions.pdf`, `figures/cdn_eq_distance_distributions.pdf`, and `figures/aa_eq_distance_distributions.pdf`. `figures/distance_pairwise_dotplots.pdf` scatters those distances against each other, with one point per species. `figures/nuc_freq_by_class_gene.pdf`, `figures/cdn_freq_by_class_gene.pdf`, and `figures/aa_freq_by_class_gene.pdf` are boxplots of observed and equilibrium frequencies for each class and gene, and for each nucleotide subset.

Fish and birds are compared with:

```bash
uv run python compare_fish_bird_equilibrium.py
```

Within each gene, every Actinopteri or Aves species that has an equilibrium is a focal species. The script records the Manhattan distance from its observed frequencies to its own equilibrium (`own`), to the equilibrium of every other species in the same class (`within`), and to the equilibrium of every species in the other class (`other`). Pairs are in `data/fish_bird_distances.csv`, which is not committed. `figures/fish_bird_nuc.pdf`, `figures/fish_bird_cdn.pdf`, and `figures/fish_bird_aa.pdf` plot, for each focal species, the own distance against its median within-class distance and its median other-class distance. Nucleotide pages are split by subset. The medians are used only in those figures and in the printed fractions.

Frequency PCA, fit on observed frequencies only, is:

```bash
uv run python pca_frequencies.py
```

Cytb and ND2 are separate panels in `figures/pca_aa_freqs.pdf`, `figures/pca_cdn_freqs.pdf`, and `figures/pca_nuc_freqs.pdf`. The fit uses Aves and Actinopteri only. Equilibrium frequencies are projected with the observed-frequency fit.

`prepare_data_for_aa_shift.ipynb` and `gene_strand_forecast.ipynb` are the earlier class-level and strand explorations. `calculate_exp_aa_subst_matrix.py` keeps its own matrix orientation; see [`docs/open_questions.md`](../../docs/open_questions.md).
