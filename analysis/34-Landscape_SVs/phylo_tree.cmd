#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --job-name=tree2
#SBATCH --output=stdout_tree.txt
#SBATCH --error=error_tree.txt
#SBATCH --time=30:00:00
#SBATCH --mem=100GB

#Load RAxML
module load GCC/11.3.0  
module load OpenMPI/4.1.4
module load RAxML/8.2.12-hybrid-avx2

cd /data/projects/punim1637/Aleix/SV/SV_analysis/data/tree

#Run RAxML
raxmlHPC \
-s /scratch/punim1637/Dunstan/SV/vcf/snps/all_samples_snps_sv_filt.fasta \
-m GTRGAMMA \
-n mtb_graph_tree_sv_filt2 \
-o GCA_012923765.1 \
-T 50 \
-p 123456 \
-N 100 \
-b 123456

