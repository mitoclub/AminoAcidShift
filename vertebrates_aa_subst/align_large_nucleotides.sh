#!/bin/bash

FASTA=$1
SUFFIX=".fa"

GENCODE=2
THREADS=16
MACSE="java -jar /opt/tools/macse_v2.07.jar"

PREFIX=$(basename $FASTA $SUFFIX)

# NT2AA
$MACSE -prog translateNT2AA -seq $FASTA -gc_def $GENCODE -out_AA $PREFIX.faa
#ALN AA
mafft --thread $THREADS $PREFIX.faa > ${PREFIX}_aln.faa
#AA_ALN --> NT_ALN
$MACSE -prog reportGapsAA2NT -align_AA ${PREFIX}_aln.faa -seq $FASTA -out_NT ${PREFIX}_aln.fna

rm $PREFIX.faa



# $MACSE -prog trimNonHomologousFragments -seq chordates_sample.dedup_aln.faa -out_trim_info output_stats.csv -min_homology_to_keep_seq 0.6 -min_trim_in 40 -min_trim_ext 20
