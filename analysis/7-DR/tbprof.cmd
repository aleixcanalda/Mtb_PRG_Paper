#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --job-name=tbprof
#SBATCH --output=std/stdout_tbprof.txt
#SBATCH --error=std/error_tbprof.txt
#SBATCH --time=01:00:00
#SBATCH --mem=40GB
#SBATCH --array=0-9999

source /apps/easybuild-2022/easybuild/software/Core/Anaconda3/2022.10/etc/profile.d/conda.sh

FILES=($(<$1))
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

cd /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/tbprof/

#while read FILE

#do

conda activate fastq-dl

fastq-dl --cpus 30 -a $FILE

rm /data/gpfs/projects/punim1637/Aleix/SV/SV_analysis/data/DR/tbprof/fastq-run-info.tsv

conda deactivate

conda activate /data/gpfs/projects/punim1637/conda_envs/tb-profiler

tb-profiler profile -1 ${FILE}_1.fastq.gz -2 ${FILE}_2.fastq.gz -t 30 -p $FILE --call_whole_genome

conda deactivate

rm ${FILE}_1.fastq.gz

rm ${FILE}_2.fastq.gz

rm /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/tbprof/bam/${FILE}*.bam*

#done < $1
