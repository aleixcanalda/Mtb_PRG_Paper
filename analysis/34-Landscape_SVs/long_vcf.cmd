#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name=lin
#SBATCH --output=stdout_lin.txt
#SBATCH --error=error_lin.txt
#SBATCH --time=10:00:00
#SBATCH --array=0-13
#SBATCH --mem=20GB

module load GCCcore/11.3.0
module load minimap2/2.24
module load SAMtools/1.16.1
module load picard/2.25.1-Java-11
module load Pilon/1.23-Java-11
module load GCC/11.3.0
module load HTSlib/1.15.1
module load tabix
module load Anaconda3/2023.07-2
module load BCFtools/1.15.1
module load Python/3.10.4

FILES=($(<$1))
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
OUT_DIR='/data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies'
OUT_DIR_FILE=/scratch/punim1637/Dunstan/SV/fastq/${FILE}
FASTQ='/scratch/punim1637/Dunstan/SV/fastq'
DB_DIR='/data/projects/punim1637/Aleix/SV/SV_analysis/data'
REF_GEN='/data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa'
VCF='/scratch/punim1637/Dunstan/SV/vcf/snps'
VCF_OUT=${VCF}/${FILE}

##VCF CREATION

#we map the reads to the reference genome
minimap2 -a -t 8 -x map-pb -o ${OUT_DIR_FILE}.sam ${REF_GEN} ${OUT_DIR_FILE}_reads_enriched.fq.gz

#we convert the alignment (sam) file to BAM
samtools view -S -b ${OUT_DIR_FILE}.sam > ${OUT_DIR_FILE}.bam

rm ${OUT_DIR_FILE}.sam

#sort the bam file
samtools sort \
-o ${OUT_DIR_FILE}.SORTED.bam \
${OUT_DIR_FILE}.bam

rm ${OUT_DIR_FILE}.bam

#we mark duplicates and remove them
java -Xmx20g -jar $EBROOTPICARD/picard.jar MarkDuplicates \
INPUT=${OUT_DIR_FILE}.SORTED.bam \
OUTPUT=${OUT_DIR_FILE}.SORTED.DUP.FREE.bam \
METRICS_FILE=${OUT_DIR_FILE}.dup_metrics \
REMOVE_DUPLICATES=true

rm ${OUT_DIR_FILE}.SORTED.bam
rm ${OUT_DIR_FILE}.dup_metrics

#we mark read groups
java -Xmx20g -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
I=${OUT_DIR_FILE}.SORTED.DUP.FREE.bam \
O=${OUT_DIR_FILE}.SORTED.DUP.FREE.W.RG.bam \
RGID=4 \
RGLB=lib1 \
RGPL=long \
RGPU=unit1 \
RGSM=20

rm ${OUT_DIR_FILE}.SORTED.DUP.FREE.bam

#we index our file, a necessary step for variant calling
samtools index ${OUT_DIR_FILE}.SORTED.DUP.FREE.W.RG.bam

#we variant call
java -Xmx20g -jar $EBROOTPILON/pilon.jar --genome /data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa --bam ${OUT_DIR_FILE}.SORTED.DUP.FREE.W.RG.bam --outdir $VCF --output ${FILE}.PILON.OUTPUT --variant --mindepth 10 --minmq 40 --minqual 20

rm ${OUT_DIR_FILE}.SORTED.DUP.FREE.W.RG.bam*
rm ${VCF_OUT}.PILON.OUTPUT.fasta

##VCF FILTERING

module load GATK/4.3.0.0-Java-11

gatk SelectVariants -V ${VCF_OUT}.PILON.OUTPUT.vcf -O ${VCF_OUT}.PILON.OUTPUT.FILTERED.VARIANT.SNP.vcf --exclude-filtered true -select-type SNP

rm ${VCF_OUT}.PILON.OUTPUT.vcf

module load BCFtools/1.15.1

bgzip ${VCF_OUT}.PILON.OUTPUT.FILTERED.VARIANT.SNP.vcf

bcftools index ${VCF_OUT}.PILON.OUTPUT.FILTERED.VARIANT.SNP.vcf.gz

bcftools filter -i 'QUAL > 20' ${VCF_OUT}.PILON.OUTPUT.FILTERED.VARIANT.SNP.vcf.gz -Oz -o ${VCF_OUT}.PILON.OUTPUT.FILTERED.VARIANT.QUAL20.SNP.vcf.gz

bcftools index ${VCF_OUT}.PILON.OUTPUT.FILTERED.VARIANT.SNP.vcf.gz
