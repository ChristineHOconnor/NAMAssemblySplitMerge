#!/bin/bash
#PBS -l mem=210gb,nodes=1:ppn=6,walltime=60:00:00
#PBS -A hirschc1
#PBS -m abe
#PBS -o /home/hirschc1/oconnorc/OandE
#PBS -e /home/hirschc1/oconnorc/OandE
#PBS -M oconnorc@umn.edu
#PBS -q mangi
#PBS -N NAM_mum_1

#load GNU parallel
module load parallel
#load mummer 4.0
module load mummer/4.0.0.beta2


cd /home/hirschc1/oconnorc/scripts/Split_genes/COversion_SplitMerge_Pipeline/Nucmer_All_by_All
export PARALLEL="--workdir . --env PATH --env LD_LIBRARY_PATH --env LOADEDMODULES --env _LMFILES_ --env MODULE_VERSION --env MODULEPATH --env MODULEVERSION_STACK --env MODULESHOME --env OMP_DYNAMICS --env OMP_MAX_ACTIVE_LEVELS --env OMP_NESTED --env OMP_NUM_THREADS --env OMP_SCHEDULE --env OMP_STACKSIZE --env OMP_THREAD_LIMIT --env OMP_WAIT_POLICY --tmpdir /home/hirschc1/oconnorc/OandE/tmp_files --compress"

sort -u $PBS_NODEFILE > unique_nodelist.txt
#go to where the parallelization command is stored
parallel --jobs 6 --sshloginfile unique_nodelist.txt --workdir $PWD < nucmer_commands_NAM_genomes_rerun1.txt


#/home/hirschc1/oconnorc/software/MUMmer3.23/nucmer --mum -c 1000 -p /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/B73_NC350_HP301_comparison/nucmer_output/B73v4_Mo17_c1000 /home/hirschc1/oconnorc/SplitGenes/second_pipeline_attempt_May2019/fasta_files/Zea_mays.AGPv4.dna.toplevel.fa /home/hirschc1/oconnorc/SplitGenes/second_pipeline_attempt_May2019/fasta_files/mo17_10chr_gbrowse_1August2018.fa

#/home/hirschc1/oconnorc/software/MUMmer3.23/nucmer --mum -c 1000 -p /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/B73_NC350_HP301_comparison/nucmer_output/B73_HP301_c1000_chr1 /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/B73_NC350_HP301_comparison/nucmer_output/Zm-B73-REFERENCE-NAM-5.0_chr1.fasta  /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/B73_NC350_HP301_comparison/nucmer_output/Zm-HP301-REFERENCE-NAM-1.0.fasta_chr1.fasta

#module load mummer/4.0.0.beta2
#nucmer --mum -c 1000 -p /home/hirschc1/oconnorc/SplitGenes/NAM_genome_comparisons/B73v4_Mo17_c1000_mum4.0 /home/hirschc1/oconnorc/SplitGenes/second_pipeline_attempt_May2019/fasta_files/Zea_mays.AGPv4.dna.toplevel.fa /home/hirschc1/oconnorc/SplitGenes/second_pipeline_attempt_May2019/fasta_files/mo17_10chr_gbrowse_1August2018.fa


