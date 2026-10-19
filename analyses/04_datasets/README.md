# Cross-dataset comparison

Question: is the neutral amino-acid prediction strongest in the mtDNA collections where mutation outruns selection?

Inputs: MitoMap, gTEX, and megatree tables in [`data/external`](../../data/external); vertebrate mutations in [`analyses/02_vertebrates/dataset/obs_muts.csv`](../02_vertebrates/dataset/obs_muts.csv); spectra from the sibling projects `192`, `human_mtDNA_megatree`, and `mtdnaMutSpecOfCancers`.

Run `analyse.ipynb` and `vert_species_fit.ipynb` from this directory. The batch fit is `uv run python analyses/04_datasets/apply_neutral_model_all_datasets.py`.

Output: `datasets_fit_metrics.csv`. A fresh run skips the mammal ND1 rows when the vertebrate spectrum file has no ND1 spectrum, and it writes centered-log-ratio correlations. The committed table is from the previous metrics function.
