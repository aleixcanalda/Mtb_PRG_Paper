#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --job-name=rnaseq
#SBATCH --output=std/stdout_rnaseq.txt
#SBATCH --error=std/error_rnaseq.txt
#SBATCH --time=15-00:00:00
#SBATCH --mem=40GB

module load Kraken2/2.1.2
module load fastp
module load MUMmer/4.0.0rc1
module load minimap2/2.26
module load Anaconda3
module load Python/3.10.4
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load GCCcore/11.3.0
module load SAMtools/1.16.1
module load BWA/0.7.17
module load SeqKit/2.3.1

source /apps/easybuild-2022/easybuild/software/Core/Anaconda3/2022.10/etc/profile.d/conda.sh

while read FILE

do

#FILES=($(<$1))
#FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
OUT_DIR='/data/projects/punim1637/Aleix/SV/SV_analysis/data/rv80/'
OUT_DIR_FILE=/data/projects/punim1637/Aleix/SV/SV_analysis/data/rv80/${FILE}
DB_DIR='/data/projects/punim1637/Aleix/SV/SV_analysis/data'
REF_GEN='/data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa'

conda activate fastq-dl

#We get the file 
fastq-dl -o $OUT_DIR -a $FILE --cpus 30
rm ${OUT_DIR}/fastq-run-info.tsv

conda deactivate

#if test -f ${FILE}_1.fastq.gz; then

	#we preprocess the short reads
	fastp -i ${OUT_DIR_FILE}_1.fastq.gz -I ${OUT_DIR_FILE}_2.fastq.gz -o ${OUT_DIR_FILE}_1_trim.fastq.gz -O ${OUT_DIR_FILE}_2_trim.fastq.gz

	rm ${OUT_DIR_FILE}_1.fastq.gz
	rm ${OUT_DIR_FILE}_2.fastq.gz
	rm fastp*

	bwa mem $REF_GEN ${OUT_DIR_FILE}_1_trim.fastq.gz ${OUT_DIR_FILE}_2_trim.fastq.gz > ${OUT_DIR_FILE}.sam

	samtools view -@ 30 -bS ${OUT_DIR_FILE}.sam | samtools sort -@ 30 -o ${OUT_DIR_FILE}.bam
#        samtools markdup -r ${OUT_DIR_FILE}.bam ${OUT_DIR_FILE}_dup.bam
	samtools index ${OUT_DIR_FILE}.bam -@ 30
	rm ${OUT_DIR_FILE}.sam

#else

#	fastp -i ${OUT_DIR_FILE}.fastq.gz -o ${OUT_DIR_FILE}_trim.fastq.gz

#	rm ${OUT_DIR_FILE}.fastq.gz
#	rm fastp*

#	bwa mem $REF_GEN ${OUT_DIR_FILE}_trim.fastq.gz > ${OUT_DIR_FILE}.sam

#	samtools view -bS ${OUT_DIR_FILE}.sam | samtools sort -o ${OUT_DIR_FILE}.bam
#	samtools markdup -r ${OUT_DIR_FILE}.bam ${OUT_DIR_FILE}_dup.bam
#	samtools index ${OUT_DIR_FILE}_dup.bam
#	rm ${OUT_DIR_FILE}.sam
#fi

done < /data/projects/punim1637/Aleix/SV/SV_analysis/data/rv80/samples.txt

conda activate rnaseq

featureCounts -a /data/gpfs/projects/punim1637/Aleix/SV/Manila/RNA-seq/data/Mycobacterium_tuberculosis_H37Rv_gff_v5_modified.gff -o ${OUT_DIR}counts.txt -s 2 -Q 10 ${OUT_DIR}/*.bam -t CDS -g Locus -p

conda deactivate
