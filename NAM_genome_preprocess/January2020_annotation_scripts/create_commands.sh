FILES=$(ls /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/*gff3)
file_base=$(basename -a $FILES | cut -d_ -f1 | sort | uniq)

#create commands to run edit gene name script 

#file_base="B73
#HP301
#NC350"

for i in $file_base
do
#    echo "python /home/hirschc1/oconnorc/scripts/Split_genes/NAM_genome_preprocessing/edit_gene_names_rem_overlapping.py -g <(sort -k 1,1 -k 7,7 -k 4,4n -k 5,5nr -k 3,3r /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/${i}_abinitio_evidence.gff3) -k /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/${i}_abinitio_evidence_gene.name.key.txt -o /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/B73_NC3250_HP301_comparison/${i}_abinitio_evidence_new.gene.names.gff3 -r /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/B73_NC3250_HP301_comparison/${i}_abinitio_evidence_new.gene.names_nested.genes.removed.gff3 -n /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/B73_NC3250_HP301_comparison/${i}_abinitio_evidence_genes_not.removed.txt -f /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/B73_NC3250_HP301_comparison/${i}_abinitio_evidence_genes_to.be.removed.txt" >> $1
#    echo "python /home/hirschc1/oconnorc/scripts/Split_genes/NAM_genome_preprocessing/edit_gene_names_rem_overlapping.py -g /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/${i}_abinitio_evidence.gff3 -k /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/${i}_abinitio_evidence_gene.name.key.txt -o /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/${i}_abinitio_evidence_new.gene.names.all.gff3 -r /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/${i}_abinitio_evidence_new.gene.names.red.gff3 -i /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/${i}_abinitio_evidence_genes.not.removed.txt -l /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/${i}_abinitio_evidence_genes_to.be.removed.txt" >> $1

#echo "bedtools getfasta -fi /home/maize/shared/databases/genomes/Zea_mays/${i}_NAMassembly/Zm-${i}-REFERENCE-NAM-1.0/Zm-${i}-REFERENCE-NAM-1.0.fasta -bed <(grep \"	CDS\" /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/${i}_abinitio_evidence_new.gene.names.red.gff3 | tr \";\" \"\t\" | cut -f 1-8,10 | awk '{print \$1\"\t\"\$2\"\t\"\$9\"\t\"\$4\"\t\"\$5\"\t\"\$6\"\t\"\$7\"\t\"\$8\"\t\"\$9}') -name -s > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/cds_fasta_files/${i}_abinitio_evidence_new.names.red.by.CDS.fasta" >> $1
#this will have to be altered for B73 since it's the 5.0 version 
#sed command that was needed in v 1.0: sed -e 's/^>.\{26,27\};/>/'
#echo "awk -f /home/hirschc1/oconnorc/scripts/Split_genes/NAM_genome_preprocessing/combine_fasta_headers.awk /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/cds_fasta_files/${i}_abinitio_evidence_new.names.red.by.CDS.fasta > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/cds_fasta_files/${i}_abinitio_evidence_new.names.red.by.CDS.combined.fasta" >> $2

#before getting longest transcript, sort gff file by chromosome, strand, start position, reverse alphabetical order and column 9 entry
#why? because parse gff script needs the exons of a gene to come after the gene 
echo "python /home/hirschc1/oconnorc/scripts/Split_genes/COversion_SplitMerge_Pipeline/Parse_GFF_LongestTranscript.py -g <(sort -k 1,1 -k 7,7 -k 4,4n -k 3,3r -k 9,9 /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/${i}_abinitio_evidence_new.gene.names.red.gff3) -o /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/longest_transcript_files/${i}_abinitio_evidence_longest_ts_exons_sorted.gff -i Parent" >> $1

#commands to run get_longest_TSname.py
#echo "python /home/hirschc1/oconnorc/scripts/Split_genes/COversion_SplitMerge_Pipeline/get_longest_TSname.py -g /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/longest_transcript_files/${i}_abinitio_evidence_longest_ts_exons_sorted.gff -o /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/longest_transcript_files/${i}_abinitio_longest_ts.name.txt" >> $1
#create files that have all the nested and non-nested genes for all the genotypes for some quick graphing in R
#echo "awk '{len=\$5-\$4; print \$1\"\t\"\$2\"\t\"\$7\"\t\"len\"\t${i}\tnested_gene\"}' /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/${i}_abinitio_evidence_genes_to.be.removed.txt >> NAM_genomes_nested.genes.txt" >> $1
#echo "awk '{len=\$5-\$4; print \$1\"\t\"\$2\"\t\"\$7\"\t\"len\"\t${i}\tnon_nested_gene\"}' /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files/${i}_abinitio_evidence_genes_not.removed.txt >> NAM_genomes_non.nested.genes.txt" >> $1
done

