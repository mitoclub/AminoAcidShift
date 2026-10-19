# Equilibrium frequencies

Question: given the vertebrate mitochondrial spectrum, which amino acids are expected to become more common, and does gene strand differ among phyla?

Inputs: `MutSpecVertebrates12.csv.gz` from the sibling `192` project (`AAS_NEMU_192`), RefSeq mtDNA in [`data/reference`](../../data/reference), and the nemu FASTA files used by `derive_obs_aa_freqs.py`.

Run `prepare_data_for_aa_shift.ipynb` and `gene_strand_forecast.ipynb` from this directory. `calculate_exp_aa_subst_matrix.py` writes `data/exp_aa_subst_*.csv`. It keeps its own matrix orientation; see [`docs/open_questions.md`](../../docs/open_questions.md).

Outputs: `data/*.csv`, `figures/*.pdf`, `data/gene_strand_counts_by_phylum.csv`, `data/species_list.txt`.
