#!/usr/bin/env bash
set -e
set -x

rm -rf aln2meta_output Ngo*

mkdir aln2meta_output

#presence_absence directory contains gene sequences for presence/absence genes along with manually created metadata and cluster files for those sequences
#aln2meta_input directory contains alignments for genes with known AMR variants along with manually created metadata tsv files for gene in which one representative of each variant is defined
#alignments were created in Seaview v? by translating to amino acid sequences, aligning with Clustal using default parameters and back translating to nucleotides
#ARIBA v2.8.1 was used to create an ARIBA database using the following commands


for x in folP gyrA mtrR parC parE penA ponA porB1b rpoB rpsJ
do
    ariba aln2meta --variant_only aln2meta_input/$x.aln aln2meta_input/$x\_in.tsv coding aln2meta_output/$x
done

for x in 16S 23S
do
    ariba aln2meta --variant_only aln2meta_input/$x.aln aln2meta_input/$x\_in.tsv noncoding aln2meta_output/$x
done


cat aln2meta_output/*.fa presence_absence/*.fa > Ngo_ARIBA.fa
cat aln2meta_output/*.tsv presence_absence/presence_absence.tsv > Ngo_ARIBA.tsv
cat aln2meta_output/*.cluster presence_absence/presence_absence.clusters > Ngo_ARIBA.clusters
ariba prepareref -f Ngo_ARIBA.fa -m Ngo_ARIBA.tsv --cdhit_clusters Ngo_ARIBA.clusters Ngo_ARIBAdb
