#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --job-name=tbprof
#SBATCH --output=std/stdout_tbprof.txt
#SBATCH --error=std/error_tbprof.txt
#SBATCH --time=01:00:00
#SBATCH --mem=40GB
#SBATCH --array=0-3145

source /apps/easybuild-2022/easybuild/software/Core/Anaconda3/2022.10/etc/profile.d/conda.sh
module load snippy/4.6.0
module load BCFtools/1.15.1

FILES=($(<$1))
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

cd /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/vcf/

#while read FILE

#do

conda activate fastq-dl

fastq-dl --cpus 30 -a $FILE

rm /data/gpfs/projects/punim1637/Aleix/SV/SV_analysis/data/DR/tbprof/fastq-run-info.tsv

conda deactivate

snippy --cpus 30 --outdir /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/vcf/${FILE} --ref /data/projects/punim1637/Shared/Reference/MTB/NCBI/H37Rv.fasta --R1 ${FILE}_1.fastq.gz --R2 ${FILE}_2.fastq.gz

rm ${FILE}_1.fastq.gz

rm ${FILE}_2.fastq.gz

mv /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/vcf/${FILE}/snps.filt.vcf /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/vcf/${FILE}.vcf

rm -r /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/vcf/${FILE}/

bcftools view -v snps /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/vcf/${FILE}.vcf -o /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/vcf/${FILE}.snps.vcf

rm /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/vcf/${FILE}.vcf

bcftools view -T ^/data/projects/punim1637/Aleix/SV/SV_analysis/data/repeats/RLC_Regions.Plus.LowPmapK50E4.H37Rv.bed /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/vcf/${FILE}.snps.vcf -o /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/vcf/${FILE}.snps.final.vcf

rm /data/projects/punim1637/Aleix/SV/SV_analysis/data/DR/vcf/${FILE}.snps.vcf
