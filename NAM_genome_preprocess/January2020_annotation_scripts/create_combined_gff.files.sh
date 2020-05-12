#commands to edit gene names and transcript name formats 
#new ZmXYZ name format comes from the Jan2020 release 
#grep "^##" /home/maize/shared/databases/genomes/Zea_mays/B73_NAMassembly/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/B73_abinitio_evidence.gff3
#cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/B73_NAMassembly/B73_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/B73_NAMassembly/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/B73_abinitio_evidence.gff3
#grep "^##" /home/maize/shared/databases/genomes/Zea_mays/B97_NAMassembly/Zm-B97-REFERENCE-NAM-1.0/Zm-B97-REFERENCE-NAM-1.0_Zm00018a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/B97_abinitio_evidence.gff3
#cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/B97_NAMassembly/B97_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/B97_NAMassembly/Zm-B97-REFERENCE-NAM-1.0/Zm-B97-REFERENCE-NAM-1.0_Zm00018a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/B97_abinitio_evidence.gff3
#grep "^##" /home/maize/shared/databases/genomes/Zea_mays/CML103_NAMassembly/Zm-CML103-REFERENCE-NAM-1.0/Zm-CML103-REFERENCE-NAM-1.0_Zm00021a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML103_abinitio_evidence.gff3
#cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML103_NAMassembly/CML103_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML103_NAMassembly/Zm-CML103-REFERENCE-NAM-1.0/Zm-CML103-REFERENCE-NAM-1.0_Zm00021a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML103_abinitio_evidence.gff3
#grep "^##" /home/maize/shared/databases/genomes/Zea_mays/CML228_NAMassembly/Zm-CML228-REFERENCE-NAM-1.0/Zm-CML228-REFERENCE-NAM-1.0_Zm00022a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML228_abinitio_evidence.gff3
#cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML228_NAMassembly/CML228_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML228_NAMassembly/Zm-CML228-REFERENCE-NAM-1.0/Zm-CML228-REFERENCE-NAM-1.0_Zm00022a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML228_abinitio_evidence.gff3
#grep "^##" /home/maize/shared/databases/genomes/Zea_mays/CML247_NAMassembly/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0_Zm00023a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML247_abinitio_evidence.gff3
#cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML247_NAMassembly/CML247_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML247_NAMassembly/Zm-CML247-REFERENCE-NAM-1.0/Zm-CML247-REFERENCE-NAM-1.0_Zm00023a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML247_abinitio_evidence.gff3
#grep "^##" /home/maize/shared/databases/genomes/Zea_mays/CML277_NAMassembly/Zm-CML277-REFERENCE-NAM-1.0/Zm-CML277-REFERENCE-NAM-1.0_Zm00024a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML277_abinitio_evidence.gff3
#cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML277_NAMassembly/CML277_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML277_NAMassembly/Zm-CML277-REFERENCE-NAM-1.0/Zm-CML277-REFERENCE-NAM-1.0_Zm00024a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML277_abinitio_evidence.gff3
#grep "^##" /home/maize/shared/databases/genomes/Zea_mays/CML322_NAMassembly/Zm-CML322-REFERENCE-NAM-1.0/Zm-CML322-REFERENCE-NAM-1.0_Zm00025a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML322_abinitio_evidence.gff3
#cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML322_NAMassembly/CML322_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML322_NAMassembly/Zm-CML322-REFERENCE-NAM-1.0/Zm-CML322-REFERENCE-NAM-1.0_Zm00025a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML322_abinitio_evidence.gff3
#grep "^##" /home/maize/shared/databases/genomes/Zea_mays/CML333_NAMassembly/Zm-CML333-REFERENCE-NAM-1.0/Zm-CML333-REFERENCE-NAM-1.0_Zm00026a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML333_abinitio_evidence.gff3
#cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML333_NAMassembly/CML333_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML333_NAMassembly/Zm-CML333-REFERENCE-NAM-1.0/Zm-CML333-REFERENCE-NAM-1.0_Zm00026a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML333_abinitio_evidence.gff3
#grep "^##" /home/maize/shared/databases/genomes/Zea_mays/CML52_NAMassembly/Zm-CML52-REFERENCE-NAM-1.0/Zm-CML52-REFERENCE-NAM-1.0_Zm00019a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML52_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML52_NAMassembly/CML52_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML52_NAMassembly/Zm-CML52-REFERENCE-NAM-1.0/Zm-CML52-REFERENCE-NAM-1.0_Zm00019a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML52_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/CML69_NAMassembly/Zm-CML69-REFERENCE-NAM-1.0/Zm-CML69-REFERENCE-NAM-1.0_Zm00020a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML69_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML69_NAMassembly/CML69_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/CML69_NAMassembly/Zm-CML69-REFERENCE-NAM-1.0/Zm-CML69-REFERENCE-NAM-1.0_Zm00020a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/CML69_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/HP301_NAMassembly/Zm-HP301-REFERENCE-NAM-1.0/Zm-HP301-REFERENCE-NAM-1.0_Zm00027a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/HP301_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/HP301_NAMassembly/HP301_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/HP301_NAMassembly/Zm-HP301-REFERENCE-NAM-1.0/Zm-HP301-REFERENCE-NAM-1.0_Zm00027a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/HP301_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/Il14H_NAMassembly/Zm-Il14H-REFERENCE-NAM-1.0/Zm-Il14H-REFERENCE-NAM-1.0_Zm00028a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Il14H_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Il14H_NAMassembly/IL14H_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Il14H_NAMassembly/Zm-Il14H-REFERENCE-NAM-1.0/Zm-Il14H-REFERENCE-NAM-1.0_Zm00028a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Il14H_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/Ki11_NAMassembly/Zm-Ki11-REFERENCE-NAM-1.0/Zm-Ki11-REFERENCE-NAM-1.0_Zm00030a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Ki11_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Ki11_NAMassembly/Ki11_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Ki11_NAMassembly/Zm-Ki11-REFERENCE-NAM-1.0/Zm-Ki11-REFERENCE-NAM-1.0_Zm00030a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Ki11_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/Ki3_NAMassembly/Zm-Ki3-REFERENCE-NAM-1.0/Zm-Ki3-REFERENCE-NAM-1.0_Zm00029a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Ki3_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Ki3_NAMassembly/Ki3_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Ki3_NAMassembly/Zm-Ki3-REFERENCE-NAM-1.0/Zm-Ki3-REFERENCE-NAM-1.0_Zm00029a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Ki3_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/Ky21_NAMassembly/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Ky21_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Ky21_NAMassembly/Ky21_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Ky21_NAMassembly/Zm-Ky21-REFERENCE-NAM-1.0/Zm-Ky21-REFERENCE-NAM-1.0_Zm00031a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Ky21_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/M162W_NAMassembly/Zm-M162W-REFERENCE-NAM-1.0/Zm-M162W-REFERENCE-NAM-1.0_Zm00033a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/M162W_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/M162W_NAMassembly/M162W_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/M162W_NAMassembly/Zm-M162W-REFERENCE-NAM-1.0/Zm-M162W-REFERENCE-NAM-1.0_Zm00033a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/M162W_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/M37W_NAMassembly/Zm-M37W-REFERENCE-NAM-1.0/Zm-M37W-REFERENCE-NAM-1.0_Zm00032a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/M37W_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/M37W_NAMassembly/M37W_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/M37W_NAMassembly/Zm-M37W-REFERENCE-NAM-1.0/Zm-M37W-REFERENCE-NAM-1.0_Zm00032a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/M37W_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/Mo18W_NAMassembly/Zm-Mo18W-REFERENCE-NAM-1.0/Zm-Mo18W-REFERENCE-NAM-1.0_Zm00034a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Mo18W_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Mo18W_NAMassembly/Mo18W_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Mo18W_NAMassembly/Zm-Mo18W-REFERENCE-NAM-1.0/Zm-Mo18W-REFERENCE-NAM-1.0_Zm00034a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Mo18W_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/Ms71_NAMassembly/Zm-Ms71-REFERENCE-NAM-1.0/Zm-Ms71-REFERENCE-NAM-1.0_Zm00035a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Ms71_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/MS71_NAMassembly/Ms71_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Ms71_NAMassembly/Zm-Ms71-REFERENCE-NAM-1.0/Zm-Ms71-REFERENCE-NAM-1.0_Zm00035a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Ms71_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/NC350_NAMassembly/Zm-NC350-REFERENCE-NAM-1.0/Zm-NC350-REFERENCE-NAM-1.0_Zm00036a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/NC350_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/NC350_NAMassembly/NC350_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/NC350_NAMassembly/Zm-NC350-REFERENCE-NAM-1.0/Zm-NC350-REFERENCE-NAM-1.0_Zm00036a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/NC350_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/NC358_NAMassembly/Zm-NC358-REFERENCE-NAM-1.0/Zm-NC358-REFERENCE-NAM-1.0_Zm00037a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/NC358_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/NC358_NAMassembly/NC358_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/NC358_NAMassembly/Zm-NC358-REFERENCE-NAM-1.0/Zm-NC358-REFERENCE-NAM-1.0_Zm00037a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/NC358_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/Oh43_NAMassembly/Zm-Oh43-REFERENCE-NAM-1.0/Zm-Oh43-REFERENCE-NAM-1.0_Zm00039a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Oh43_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Oh43_NAMassembly/Oh43_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Oh43_NAMassembly/Zm-Oh43-REFERENCE-NAM-1.0/Zm-Oh43-REFERENCE-NAM-1.0_Zm00039a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Oh43_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/Oh7B_NAMassembly/Zm-Oh7B-REFERENCE-NAM-1.0/Zm-Oh7B-REFERENCE-NAM-1.0_Zm00038a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Oh7B_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Oh7B_NAMassembly/Oh7b_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Oh7B_NAMassembly/Zm-Oh7B-REFERENCE-NAM-1.0/Zm-Oh7B-REFERENCE-NAM-1.0_Zm00038a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Oh7B_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/P39_NAMassembly/Zm-P39-REFERENCE-NAM-1.0/Zm-P39-REFERENCE-NAM-1.0_Zm00040a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/P39_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/P39_NAMassembly/P39_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/P39_NAMassembly/Zm-P39-REFERENCE-NAM-1.0/Zm-P39-REFERENCE-NAM-1.0_Zm00040a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/P39_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/Tx303_NAMassembly/Zm-Tx303-REFERENCE-NAM-1.0/Zm-Tx303-REFERENCE-NAM-1.0_Zm00041a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Tx303_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Tx303_NAMassembly/Tx303_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Tx303_NAMassembly/Zm-Tx303-REFERENCE-NAM-1.0/Zm-Tx303-REFERENCE-NAM-1.0_Zm00041a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Tx303_abinitio_evidence.gff3
grep "^##" /home/maize/shared/databases/genomes/Zea_mays/Tzi8_NAMassembly/Zm-Tzi8-REFERENCE-NAM-1.0/Zm-Tzi8-REFERENCE-NAM-1.0_Zm00042a.1.gff > /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Tzi8_abinitio_evidence.gff3
cat <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Tzi8_NAMassembly/Tzi8_braker-non-overlapping.gff3) <(grep -v "^#" /home/maize/shared/databases/genomes/Zea_mays/Tzi8_NAMassembly/Zm-Tzi8-REFERENCE-NAM-1.0/Zm-Tzi8-REFERENCE-NAM-1.0_Zm00042a.1.gff) | sort -k 1,1 -k 7,7 -k 4,4n >> /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/gff_files_oldname/Tzi8_abinitio_evidence.gff3