#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --job-name=pansr
#SBATCH --output=std/stdout_pansr.txt
#SBATCH --error=std/error_pansr.txt
#SBATCH --time=10-00:00:00
#SBATCH --mem=200GB

module load Kraken2/2.1.2
module load fastp
module load MUMmer/4.0.0rc1
module load foss/2022a
module load Anaconda3
module load Python/3.10.4
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load GCCcore/11.3.0
module load minimap2/2.24
module load SAMtools/1.16.1
module load BWA/0.7.17
module load bwa-mem2/2.2.1
module load manta/1.6.0
module load foss
module load minigraph
module load k8/1.0
module load shovill/1.1.0
module load SeqKit/2.3.1

source /apps/easybuild-2022/easybuild/software/Core/Anaconda3/2022.10/etc/profile.d/conda.sh

while read FILE

do

OUT_DIR='SV/vcf/short_reads/'
OUT_DIR_FILE=SV/vcf/short_reads/${FILE}
PAN_DIR=Aleix/SV/SV_analysis/data/long_read/assemblies
SCRIPT=Aleix/SV/SV_analysis/cmd/aim1.4
DB_DIR='Aleix/SV/SV_analysis/data'
REF_GEN='TB.REF.GENOME.H37RV.ONELINE.fa'

conda activate fastq-dl

#We get the file 
fastq-dl --cpus 50 -o $OUT_DIR -a $FILE
rm ${OUT_DIR}/fastq-run-info.tsv

conda deactivate

#we preprocess the short reads
fastp -w 60 -i ${OUT_DIR_FILE}_1.fastq.gz -I ${OUT_DIR_FILE}_2.fastq.gz -o ${OUT_DIR_FILE}_1_trim.fastq.gz -O ${OUT_DIR_FILE}_2_trim.fastq.gz

rm ${OUT_DIR_FILE}_1.fastq.gz
rm ${OUT_DIR_FILE}_2.fastq.gz
rm fastp*

gunzip ${OUT_DIR_FILE}_2_trim.fastq.gz
gunzip ${OUT_DIR_FILE}_1_trim.fastq.gz

#we decontaminate our samples by extracting human reads and mtb reads
kraken2 --paired --threads 60 --db ${DB_DIR}/HPRC_db/db/ --memory-mapping --output ${OUT_DIR_FILE}_classifications.tsv ${OUT_DIR_FILE}_1_trim.fastq ${OUT_DIR_FILE}_2_trim.fastq

awk -F'\t' '$1=="U" {print $2}' ${OUT_DIR_FILE}_classifications.tsv > ${OUT_DIR_FILE}_ids.txt
seqkit grep -f ${OUT_DIR_FILE}_ids.txt -o ${OUT_DIR_FILE}_1_depleted.fq ${OUT_DIR_FILE}_1_trim.fastq
seqkit grep -f ${OUT_DIR_FILE}_ids.txt -o ${OUT_DIR_FILE}_2_depleted.fq ${OUT_DIR_FILE}_2_trim.fastq

minimap2 --secondary=no -c -t 60 -x sr -o ${OUT_DIR_FILE}_reads.aln.paf ${DB_DIR}/Mycobacterium_db/Mycobacterium.rep.fna.gz ${OUT_DIR_FILE}_1_depleted.fq ${OUT_DIR_FILE}_2_depleted.fq

grep -Ff ${DB_DIR}/Mycobacterium_db/mtb.ids ${OUT_DIR_FILE}_reads.aln.paf | cut -f1 > ${OUT_DIR_FILE}_keep.ids
seqkit grep -f ${OUT_DIR_FILE}_keep.ids -o ${OUT_DIR_FILE}_1_enriched.fq ${OUT_DIR_FILE}_1_depleted.fq
seqkit grep -f ${OUT_DIR_FILE}_keep.ids -o ${OUT_DIR_FILE}_2_enriched.fq ${OUT_DIR_FILE}_2_depleted.fq

rm ${OUT_DIR_FILE}_1_depleted.fq
rm ${OUT_DIR_FILE}_2_depleted.fq
rm ${OUT_DIR_FILE}_reads.aln.paf
rm ${OUT_DIR_FILE}_classifications.tsv
rm ${OUT_DIR_FILE}_keep.ids
rm ${OUT_DIR_FILE}_ids.txt

mv ${OUT_DIR_FILE}_1_enriched.fq ${OUT_DIR_FILE}_1_trim.fastq
mv ${OUT_DIR_FILE}_2_enriched.fq ${OUT_DIR_FILE}_2_trim.fastq

gzip ${OUT_DIR_FILE}_1_trim.fastq
gzip ${OUT_DIR_FILE}_2_trim.fastq

shovill --trim --force --gsize 4.4M --cpus 30 --outdir ${OUT_DIR_FILE}_asm --R1 ${OUT_DIR_FILE}_1_trim.fastq.gz --R2 ${OUT_DIR_FILE}_2_trim.fastq.gz
rm ${OUT_DIR_FILE}_1.fastq.gz
rm ${OUT_DIR_FILE}_2.fastq.gz
mv ${OUT_DIR_FILE}_asm/contigs.fa ${OUT_DIR_FILE}_asm.fasta

rm -r ${OUT_DIR_FILE}_asm/
rm ${OUT_DIR_FILE}_1_trim.fastq.gz
rm ${OUT_DIR_FILE}_2_trim.fastq.gz

module load Python/3.10.4
virtualenv ~/venvs/venv-3.10.4
source ~/venvs/venv-3.10.4/bin/activate

minigraph -l 10000 -d 5000 -t 50 -cxasm --call ${PAN_DIR}/Mtb_PRG_NN.gfa ${OUT_DIR_FILE}_asm.fasta > ${OUT_DIR_FILE}_map2pan_short.bed

paste ${PAN_DIR}/h37rv_NN.bed ${OUT_DIR_FILE}_map2pan_short.bed | k8 /data/projects/punim1637/Aleix/mgutils.js merge - > ${OUT_DIR_FILE}_merged_short.bed

echo NC_000962.3 > ${OUT_DIR_FILE}_samples.txt
echo $FILE >> ${OUT_DIR_FILE}_samples.txt

k8 /data/projects/punim1637/Aleix/mgutils-es6.js merge2vcf -s ${OUT_DIR_FILE}_samples.txt ${OUT_DIR_FILE}_merged_short.bed > ${OUT_DIR_FILE}_pan_short.vcf

rm ${OUT_DIR_FILE}_merged_short.bed
rm ${OUT_DIR_FILE}_samples.txt

#Now we obtain the final vcf versions
miniwalk mod -b ${OUT_DIR_FILE}_map2pan_short.bed -v ${OUT_DIR_FILE}_pan_short.vcf -g ${PAN_DIR}/Mtb_PRG_NN.gfa -o ${OUT_DIR_FILE}_mod_short.vcf -na

rm ${OUT_DIR_FILE}_pan_short.vcf

miniwalk ref -v ${OUT_DIR_FILE}_mod_short.vcf -o ${OUT_DIR_FILE}_ref_short.vcf

rm ${OUT_DIR_FILE}_mod_short.vcf
rm ${OUT_DIR_FILE}_asm.fasta
rm ${OUT_DIR_FILE}_map2pan_short.bed

miniwalk ins2dup -c ${OUT_DIR_FILE}_ref_short.vcf -r /data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa

rm ${OUT_DIR_FILE}_ref_short.vcf.fa

done < $1
