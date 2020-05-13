# Split Merge Pipeline for the NAM genome assemblies 

Running Meesh's Split Merge pipeline on the maize NAM genome assemblies. Original code here: https://github.com/HirschLabUMN/Split_genes/tree/master/Split_Merge_Pipeline

Where are the data from? 
Downloaded on January 14,2020 from cyverse

Meesh's pipeline should work as is, but the NAM assemblies require some pre-processing 

## Data pre-processing steps 

We will use the evidence based and abinitio genes to run the split merge pipeline. The provided gff and CDS fasta files require some pre-processing. 

1) Create gene name keys for every gene in a \*mikado.braker.pasa.aed_lt_75.androp_abv.gff file. This is necessary because the evidence based and abinito gene names have different naming schemes and for the split merge pipeline the gene names within a genome need to be numbered relative to the order they are in in a gff file. The new gene names will start with Zm1xyz and count up from 1. The Zm number used in the evidence based gene names will be used here. Zm1xyz instead of Zm0xyz will indicate new names. Between the Zm number and gene number will be a 2 character seperator: ag, bg etc. 
create_gene_name_keys.sh has commands for all 26 assemblies; an example is below. 

```grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/B73_NAMassembly/B73_gene_annotations_20200413/b73.mikado.braker.pasa.aed_lt_75.androp_abv.gff | grep "	gene" | cut -f 9 | sed -e 's/gene://' | tr "=" "\t" | tr ";" "\t" | cut -f 2 | awk '{print $1"\t"NR}' | awk -F "\t" '{printf "%s\tZm10001ag%06i\n", $1,$2}' > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/B73_abinitio_evidence_gene.name.key.txt  ```

2) Replace all old names with new gene names and write results to gff file using edit_gene_names_2020.04annotations.py. With the new annotations the issue with overlapping/nested genes was fixed. 

```python /home/hirschc1/oconnorc/scripts/Split_genes/NAM_genome_preprocessing/NAM_preprocessing_20200413_annotations/edit_gene_names_2020.04annotations.py -g /home/maize/shared/databases/genomes/Zea_mays/B73_NAMassembly/B73_gene_annotations_20200413/b73.mikado.braker.pasa.aed_lt_75.androp_abv.gff -k /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/B73_abinitio_evidence_gene.name.key.txt -o /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/B73_evidence.abinitio_annotation_v2.gff ```

3) Edit CDS fasta file gene names using edit_cds_fasta_names.py. The names of the CDS fasta file also need to be updated. 

``` python edit_cds_fasta_names.py -f /home/maize/shared/databases/genomes/Zea_mays/B73_NAMassembly/B73_gene_annotations_20200413/b73.mikado.braker.pasa.aed_lt_75.androp_abv.cds.fasta -k ~/SplitGenes/NAM_genome_comparisons/gff_files_oldname/B73_abinitio_evidence_gene.name.key.txt -o /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/cds_fasta_files/B73_new.annotations.rename_cds.fasta ```

4) Select canonical transcripts from the CDS fasta file. Canonical transcripts are pre-defined, and are not necessary the longest transcript. Since we have a list of canonical transcripts we 

``` perl ~/scripts/Split_genes/pull_out_canonical_transcript_coded.pl -i /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/cds_fasta_files/B73_new.annotations.rename_cds.fasta -l /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/longest_transcript_files/B73_canonical_transcripts.txt -o /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/longest_transcript_files/B73_canonical_ts.fasta ```

5) Select canonical transcript gene annotations from the gff file (the version with the new gene names) and save just the gene annotation lines. 

``` perl ~/scripts/Split_genes/pull_out_canonical_transcript_coded_from_gff.pl -i /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/B73_evidence.abinitio_annotation_v2.gff -l  /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/longest_transcript_files/B73_canonical_transcripts.txt -o /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/B73_evidence.abinitio_annotation_canonical_ts.gff ```

## Split-Merge pipeline steps 

These are the step as they were run for the NAM genome assembly. There is no need here to worry about finding the longest trancript gene and CDS fasta sequence, as those files were created above using the canonical transcript information. For a full description of the original pipeline see this repository: https://github.com/HirschLabUMN/Split_genes/tree/master/Split_Merge_Pipeline

1) Run nucmer mummer on fasta files. mummer 4.0 beta was used. 

```nucmer --mum -c 1000 -p /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/nucmer_output/B73_B97_c1000 /home/maize/shared/databases/genomes/Zea_mays/B73_NAMassembly/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fasta /home/maize/shared/databases/genomes/Zea_mays/B97_NAMassembly/Zm-B97-REFERENCE-NAM-1.0/Zm-B97-REFERENCE-NAM-1.0.fasta ```

2) Create blast databases from the CDS fasta file. Only use canonical transcripts in the blast database. 

```makeblastdb -in /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/longest_transcript_files/B73_canonical_ts.fasta -out B73_cds_db -dbtype nucl ```

3) Run All_by_All_Blast_COedits.py. This is 99.9% identical to the original script, but a couple of edits needed to be made because of some slight differences in gene name formats in the gff file. 
I used all_by_all_propogator_NAM.sh to create the commands needed to run the all by all blast 650 times (26x25). The extra 'test' lines are added in so the command propogator loop works correctly, and lines with 'test' in them were excluded from any actual runs. 

``` python3 All_by_All_Blast_COedits.py -q /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/B73_evidence.abinitio_annotation_canonical_ts.gff -s /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/B97_evidence.abinitio_annotation_canonical_ts.gff -b /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/blast_dir/B97_cds_db -l /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/longest_transcript_files/B73_canonical_ts.fasta -o /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/all_by_all_output/B73_B97_AllbyAll_res.txt -n /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/nucmer_output/B73_B97_c1000.fil.coords -g bg ```
