#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --job-name=mats
#SBATCH --output=std/stdout_mats.txt
#SBATCH --error=std/error_mats.txt
#SBATCH --time=13-00:00:00
#SBATCH --mem=300GB

module load R/4.1.2

Rscript /data/projects/punim1637/Aleix/SV/SV_analysis/cmd/aim2.1/DR_matrix.R
