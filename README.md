# Split Merge Pipeline for the NAM genome assemblies 

Running Meesh's Split Merge pipeline on the maize NAM genome assemblies. Original code here: https://github.com/HirschLabUMN/Split_genes/tree/master/Split_Merge_Pipeline

Where are the data from? 
Downloaded on January 14,2020 from cyverse
All data are publically available \[somewhere\]

Meesh's pipeline should work as is, but the NAM assemblies require some pre-processing 

## Data pre-processing steps 

We will use the evidence based and abinitio genes to run the split merge pipeline. The provided gff and CDS fasta files require some pre-processing. 

1) Create gene name keys for every gene in a \*mikado.braker.pasa.aed_lt_75.androp_abv.gff file. This is necessary because the evidence based and abinito gene names have different naming schemes and for the split merge pipeline the gene names within a genome need to be numbered relative to the order they are in in a gff file. The new gene names will start with Zm1xyz and count up from 1. The Zm number used in the evidence based gene names will be used here. Zm1xyz instead of Zm0xyz will indicate new names. Between the Zm number and gene number will be a 2 character seperator: ag, bg etc. 
2) Replace all old names with new gene names and write results to gff file using edit_gene_names.py
3) Edit CDS fasta file gene names.
4) Select canonical transcripts from the CDS fasta file. Canonical transcripts are pre-defined, and are not necessary the longest transcript. 
5) Select canonical transcript gene annotations from the gff file (the version with the new gene names) and save just the gene annotation lines. 

## Split-Merge pipeline steps 

1) Run nucmer mummer on fasta files. mummer 4.0 beta was used. 
