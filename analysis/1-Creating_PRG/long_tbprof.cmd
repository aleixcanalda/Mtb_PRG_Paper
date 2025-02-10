#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=quast_ncbi
#SBATCH --output=stdout_quast.txt
#SBATCH --error=error_quast.txt
#SBATCH --time=24:00:00
#SBATCH --array=0-114
#SBATCH --mem=15GB

module load GCCcore/11.3.0
module load minimap2/2.24
module load SAMtools/1.16.1
module load picard/2.25.1-Java-11
module load Pilon/1.23-Java-11
module load GCC/11.3.0
module load HTSlib/1.15.1
module load tabix
module load Anaconda3/2023.07-2

source /apps/easybuild-2022/easybuild/software/Core/Anaconda3/2022.10/etc/profile.d/conda.sh

FILES=($(<$1))
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
OUT_DIR=/data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies_ncbi/${FILE}
OUT_DIR_FILE=/scratch/punim1637/Dunstan/SV/fastq/${FILE}
FASTQ='/scratch/punim1637/Dunstan/SV/fastq'
DB_DIR='/data/projects/punim1637/Aleix/SV/SV_analysis/data'
REF_GEN='/data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa'
VCF='/scratch/punim1637/Dunstan/SV/vcf'
VCF_OUT=${VCF}/${FILE}

eval "$(conda shell.bash hook)"

conda activate /data/gpfs/projects/punim1637/conda_envs/tb-profiler

tb-profiler profile --ram 15 -t 12 --temp $VCF --prefix ${FILE}_lineage --read1 ${OUT_DIR_FILE}_reads_enriched.fq.gz --dir $FASTQ --platform nanopore

conda deactivate
