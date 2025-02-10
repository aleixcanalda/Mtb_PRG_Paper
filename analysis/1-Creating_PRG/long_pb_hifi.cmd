#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=pbhifi
#SBATCH --output=std/stdout_pbhifi.txt
#SBATCH --error=std/error_pbhifi.txt
#SBATCH --time=04:00:00
#SBATCH --array=0-32
#SBATCH --mem=40GB

module load fastp/0.23.2
module load foss/2022a
module load Anaconda3
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load Kraken2/2.1.2
module load GCCcore/11.3.0
module load minimap2/2.24
module load SeqKit/2.3.1
module load SAMtools/1.16.1
module load BLAST+/2.13.0
module load BLAST/2.14.0-Linux_x86_64
module load BWA/0.7.17
module load bwa-mem2/2.2.1
module load minimap2/2.24
module load SAMtools/1.16.1
module load picard/2.25.1-Java-11
module load Pilon/1.23-Java-11


source /apps/easybuild-2022/easybuild/software/Core/Anaconda3/2022.10/etc/profile.d/conda.sh

FILES=($(<$1))
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
OUT_DIR='/data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies'
OUT_DIR_FILE=/scratch/punim1637/Dunstan/SV/fastq/${FILE}
DB_DIR='/data/projects/punim1637/Aleix/SV/SV_analysis/data'
FASTQ='/scratch/punim1637/Dunstan/SV/fastq'
Hifi='/data/projects/punim1637/Aleix/HiFiAdapterFilt'

#conda activate fastq-dl

#We get the file 
fastq-dl -o $FASTQ -a $FILE
rm ${FASTQ}/fastq-run-info.tsv

conda deactivate

mv ${OUT_DIR_FILE}*.fastq.gz ${OUT_DIR_FILE}_pass.fastq.gz
cd /scratch/punim1637/Dunstan/SV/fastq/
bash ${Hifi}/pbadapterfilt.sh -p ${FILE}_pass -t 8 -o ${FASTQ}
rm ${OUT_DIR_FILE}_pass.contaminant.blastout
rm ${OUT_DIR_FILE}_pass.blocklist
rm ${OUT_DIR_FILE}_pass.stats

conda activate flye

#we decontaminate our samples by extracting human reads and mtb reads
kraken2 --threads 8 --db ${DB_DIR}/HPRC_db/db/ --memory-mapping --output ${OUT_DIR_FILE}_classifications.tsv ${OUT_DIR_FILE}_pass.fastq.gz

awk -F'\t' '$1=="U" {print $2}' ${OUT_DIR_FILE}_classifications.tsv | \
  seqkit grep -f - -o ${OUT_DIR_FILE}_reads_depleted.fq ${OUT_DIR_FILE}_pass.fastq.gz

minimap2 --secondary=no -c -t 8 -x map-hifi -o ${OUT_DIR_FILE}_reads.aln.paf ${DB_DIR}/Mycobacterium_db/Mycobacterium.rep.fna.gz ${OUT_DIR_FILE}_reads_depleted.fq

grep -Ff ${DB_DIR}/Mycobacterium_db/mtb.ids ${OUT_DIR_FILE}_reads.aln.paf | cut -f1 | \
  seqkit grep -f - -o ${OUT_DIR_FILE}_reads_enriched.fq ${OUT_DIR_FILE}_reads_depleted.fq

rm ${OUT_DIR_FILE}_reads_depleted.fq
rm ${OUT_DIR_FILE}_reads.aln.paf
rm ${OUT_DIR_FILE}.filt.fastq.gz
rm ${OUT_DIR_FILE}_classifications.tsv

gzip -f ${OUT_DIR_FILE}_reads_enriched.fq

#We create the assembly
flye --pacbio-hifi ${OUT_DIR_FILE}_reads_enriched.fq.gz --out-dir ${OUT_DIR_FILE} --genome-size 4.4m -t 8 --asm-coverage 40
mv ${OUT_DIR_FILE}/assembly.fasta ${OUT_DIR}/${FILE}.fasta
rm -r ${OUT_DIR_FILE}

conda deactivate

#now we're going to polish the assembly
bwa index ${OUT_DIR}/${FILE}.fasta

minimap2 -a -t 8 -x map-hifi -o ${OUT_DIR_FILE}_reads_enriched.sam ${OUT_DIR}/${FILE}.fasta ${OUT_DIR_FILE}_reads_enriched.fq.gz

samtools view -Sb ${OUT_DIR_FILE}_reads_enriched.sam > ${OUT_DIR_FILE}_reads_enriched.bam
rm ${OUT_DIR_FILE}_reads_enriched.sam
samtools sort  ${OUT_DIR_FILE}_reads_enriched.bam > ${OUT_DIR_FILE}_reads_enriched.sort.bam
rm ${OUT_DIR_FILE}_reads_enriched.bam
samtools index ${OUT_DIR_FILE}_reads_enriched.sort.bam

java -Xmx20G -jar $EBROOTPILON/pilon.jar --genome ${OUT_DIR}/${FILE}.fasta --bam ${OUT_DIR_FILE}_reads_enriched.sort.bam --threads 8 --output ${FILE} --outdir ${OUT_DIR}/${FILE}_pilon

mv ${OUT_DIR}/${FILE}_pilon/${FILE}.fasta ${OUT_DIR}/${FILE}_polished.fasta

rm -r ${OUT_DIR}/${FILE}_pilon
rm ${OUT_DIR_FILE}_reads_enriched.sort.bam
rm ${OUT_DIR_FILE}_reads_enriched.sort.bam.bai
rm ${OUT_DIR}/${FILE}.fasta.*

bash /data/projects/punim1637/Aleix/sequence-stats/src/sequence-stats -f ${OUT_DIR}/${FILE}_polished.fasta 1000 > ${OUT_DIR}/${FILE}_polished_filtered.fasta
rm ${OUT_DIR}/${FILE}_polished.fasta

quast ${OUT_DIR}/${FILE}_polished_filtered.fasta \
        -r /data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa \
        --single ${OUT_DIR_FILE}_reads_enriched.fq.gz -t 18 -g /data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff \
        -o ${OUT_DIR_FILE}_quast

mv ${OUT_DIR_FILE}_quast/report.tsv ${OUT_DIR}/QC/${FILE}_quast_general_report.tsv

rm -r  ${OUT_DIR_FILE}_quast
rm ${OUT_DIR_FILE}.fasta
