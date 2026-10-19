# Vertebrate substitutions

Question: what nonsynonymous mtDNA substitutions are observed across vertebrate species, and how should species be sampled for comparison?

Inputs: alignments and observed-mutation tables from `nemu-pipeline` (`AAS_NEMU_PIPELINE`, default `../nemu-pipeline` next to this repository). Two local files are not in git: a Cytb FASTA symlink `MIDORI2_LONGEST_NUC_GB264_Cytb_BLAST.fasta`, and the chordate global-tree workdir that used to be linked as `globa_tree_dataset_wd` (previously `/home/dolphin/nemu_workdir/...`).

Run `collect_nonsyn_substitutions_vertebrates.ipynb` and `prepare_for_comp_sp.ipynb` from this directory. Nucleotide alignment of large Cytb sets is `align_large_nucleotides.sh` (MACSE at `/opt/tools/macse_v2.07.jar`).

Outputs: `dataset/obs_muts.csv`, `dataset/info.csv`, `msa_ncols_distributions.pdf`.
