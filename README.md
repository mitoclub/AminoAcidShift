# AminoAcidShift

The mitochondrial light-strand spectrum is strongly asymmetric: G to A and T to C are in excess. That spectrum predicts a long-term flow from loser amino acids (Leu, Phe, Cys, Trp, Gly, Val) toward gainer amino acids (Pro, His, Gln, Asn, Lys, Thr). This repository tests that prediction where mutation is weak relative to selection and where mutation outruns it.

## Layout

Shared code is the `aminoacidshift` package in [`src/aminoacidshift`](src/aminoacidshift). Each scientific question has its own folder. Run a notebook or script with that folder as the working directory, or run the scripts directly: they locate data through `aminoacidshift.paths`.

| Directory | Question |
| --- | --- |
| [`analyses/01_equilibrium`](analyses/01_equilibrium) | Equilibrium nucleotide, codon, and amino-acid frequencies, and how gene strand is arranged across taxa |
| [`analyses/02_vertebrates`](analyses/02_vertebrates) | Nonsynonymous mtDNA substitutions across vertebrate species |
| [`analyses/03_human_germline`](analyses/03_human_germline) | Human germline mtDNA, gene by gene, from the megatree |
| [`analyses/04_datasets`](analyses/04_datasets) | One neutral-model comparison across mtDNA datasets |
| [`analyses/05_pcawg`](analyses/05_pcawg) | Nuclear somatic missense mutations (PCAWG), the control for that comparison |
| [`analyses/archive`](analyses/archive) | Viral spectrum script kept out of the mtDNA paper |
| [`data/external`](data/external) | Downloaded tables (MitoMap, gTEX, megatree) |
| [`data/reference`](data/reference) | Human and RefSeq mitochondrial GenBank files, species lineages |
| [`site`](site) | Short Quarto pages for the current figures and tables |

Sibling projects are not vendored. Override them with `AAS_NEMU_192`, `AAS_MEGATREE`, `AAS_CANCER_MUTSPEC`, or `AAS_NEMU_PIPELINE` when the default relative paths are wrong. `raw_human_megatree.tsv` is gitignored.

Saved metric tables and PDFs are the current results. Re-running the fit scripts will not reproduce those files exactly: `calc_metrics` now reports centered-log-ratio correlation, and a few corrections are listed in [`docs/open_questions.md`](docs/open_questions.md).

## Setup

```bash
uv sync --extra notebooks --extra site
uv run jupyter lab
```

Python is 3.13. Mutational spectra come from `pymutspec` 0.0.15, which no longer depends on `ete3`.

The count-model helpers in `aminoacidshift.glm` are an adaptation of the E. coli analysis from Cano et al. ([PNAS 2022](https://doi.org/10.1073/pnas.2119720119), [Mutbias2022](https://github.com/alejvcano/Mutbias2022)). Refit that table with `uv run python -m aminoacidshift.glm --input dataEcoli.txt`.

Preview the website, without refitting models, from the repository root:

```bash
uv run python site/build_interactive.py
quarto preview
```
