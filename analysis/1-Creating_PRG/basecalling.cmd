#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name=bc
#SBATCH --output=stdout_bc.txt
#SBATCH --error=error_bc.txt
#SBATCH --time=02:00:00
#SBATCH --array=0
#SBATCH --mem=15GB
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu-a100

module load GCCcore/11.3.0
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load foss/2022a
module load SAMtools/1.16.1

#FILES=($(<$1))
#FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
dorado=../Aleix/dorado-0.5.3-linux-x64/bin/dorado
POD5=../SV/fastq/LRS_ALEIX/Aleix_20240307_LRSeq_3/no_sample/20240307_1407_X1_FAU96689_47244a38/pod5_pass/
OUT_DIR=../SV/fastq/LRS_ALEIX/Aleix_20240307_LRSeq_3/no_sample/20240307_1407_X1_FAU96689_47244a38/dorado_fastq
MODEL=../SV/fastq/LRS_ALEIX/Aleix_20240115_LRSeq_1/no_sample/20240116_1246_X1_FAU96707_f2d4a1df/dorado_fastq/dna_r10.4.1_e8.2_400bps_sup@v4.3.0

$dorado basecaller --kit-name EXP-NBD104 $MODEL $POD5 > $OUT_DIR/dorado.bam

$dorado basecaller --kit-name EXP-NBD114 $MODEL $POD5 > $OUT_DIR/dorado2.bam

$dorado demux --output-dir ${OUT_DIR}/demux --no-classify $OUT_DIR/dorado.bam

$dorado demux --output-dir ${OUT_DIR}/demux --no-classify $OUT_DIR/dorado2.bam

samtools fastq -T '*' ${OUT_DIR}/dorado.bam > ${OUT_DIR}/dorado.fq


