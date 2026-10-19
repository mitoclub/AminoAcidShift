# Human germline mtDNA

Question: do germline amino-acid substitutions on the human mitochondrial tree follow the spectrum of each gene?

Inputs: [`data/external/human_megatree_genes_spectra.csv`](../../data/external/human_megatree_genes_spectra.csv), gitignored `data/external/raw_human_megatree.tsv`, and [`data/reference/NC_012920.1.gb`](../../data/reference/NC_012920.1.gb).

Run `human_germline.ipynb` from this directory, or `uv run python analyses/03_human_germline/apply_neutral_model_human_megatree.py` from the repository root.

Outputs: `megatree_fit_metrics.csv`, `neutral_model_fit_megatree.pdf`, and the sliding-window, R², amino-acid-flow, and PCA PDFs in this folder.
