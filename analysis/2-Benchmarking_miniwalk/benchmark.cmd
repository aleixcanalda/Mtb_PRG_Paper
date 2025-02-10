#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=pan
#SBATCH --output=stdout_pan.txt
#SBATCH --error=error_pan.txt
#SBATCH --time=00:40:00
#SBATCH --array=0-16

module load Kraken2/2.1.2
module load fastp
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

for FILE in $(cat gold_ass.txt)

do

OUT_DIR='/data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/gold_assemblies'
OUT_DIR_FILE=/data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/gold_assemblies/${FILE}
PAN_DIR=/data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies
SCRIPT=/data/projects/punim1637/Aleix/SV/SV_analysis/cmd/aim1.4

#First we get the golden assembly's truth variants

minimap2 -a -x asm5 --cs -r2k -t 4 /data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa ${OUT_DIR_FILE}.fasta > ${OUT_DIR_FILE}.sam
samtools sort -m4G -@4 -o ${OUT_DIR_FILE}.bam ${OUT_DIR_FILE}.sam
samtools index ${OUT_DIR_FILE}.bam
rm ${OUT_DIR_FILE}.sam

conda activate svim-asm
svim-asm haploid ${OUT_DIR} ${OUT_DIR_FILE}.bam /data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa
deactivate
mv ${OUT_DIR}/variants.vcf ${OUT_DIR_FILE}.vcf
rm ${OUT_DIR}/sv-lengths.png

#We then get the SVs from Manta

#we preprocess the short reads
fastp -i ${OUT_DIR_FILE}_1.fastq.gz -I ${OUT_DIR_FILE}_2.fastq.gz -o ${OUT_DIR_FILE}_1_trim.fastq.gz -O ${OUT_DIR_FILE}_2_trim.fastq.gz

rm ${OUT_DIR_FILE}_1.fastq.gz
rm ${OUT_DIR_FILE}_2.fastq.gz

#We map the short reads to the reference genome
bwa-mem2.avx mem /data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa ${OUT_DIR_FILE}_1_trim.fastq.gz ${OUT_DIR_FILE}_2_trim.fastq.gz > ${OUT_DIR_FILE}_short.sam

samtools sort -m4G -@4 -o ${OUT_DIR_FILE}_short.bam ${OUT_DIR_FILE}_short.sam
samtools index ${OUT_DIR_FILE}_short.bam

rm ${OUT_DIR_FILE}_short.sam

#we call SVs from the short read data
module load Python/2.7.18

configManta.py --bam=${OUT_DIR_FILE}_short.bam --referenceFasta=/data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa --runDir=${OUT_DIR_FILE}_manta
${OUT_DIR_FILE}_manta/runWorkflow.py

mv ${OUT_DIR_FILE}_manta/results/variants/candidateSV.vcf.gz ${OUT_DIR_FILE}_short_manta.vcf.gz
rm -r ${OUT_DIR_FILE}_manta/
gunzip ${OUT_DIR_FILE}_short_manta.vcf.gz

#Now we get the SVs using the PRG - we need to create an assembly first

shovill --gsize 4.4M --outdir ${OUT_DIR_FILE}_asm --R1 ${OUT_DIR_FILE}_1_trim.fastq.gz --R2 ${OUT_DIR_FILE}_2_trim.fastq.gz

mv ${OUT_DIR_FILE}_asm/contigs.fa ${OUT_DIR_FILE}_asm.fasta

rm -r ${OUT_DIR_FILE}_asm/

module load Python/3.10.4
virtualenv ~/venvs/venv-3.10.4
source ~/venvs/venv-3.10.4/bin/activate

minigraph -cxasm --call ${PAN_DIR}/Mtb_PRG.gfa ${OUT_DIR_FILE}_asm.fasta > ${OUT_DIR_FILE}_map2pan_short.bed

paste ${PAN_DIR}/h37rv.bed ${OUT_DIR_FILE}_map2pan_short.bed | k8 /data/projects/punim1637/Aleix/mgutils.js merge - > ${OUT_DIR_FILE}_merged_short.bed

echo NC_000962.3 > samples.txt
echo $FILE >> samples.txt

k8 /data/projects/punim1637/Aleix/mgutils-es6.js merge2vcf -s samples.txt ${OUT_DIR_FILE}_merged_short.bed > ${OUT_DIR_FILE}_pan_short.vcf

rm ${OUT_DIR_FILE}_merged_short.bed
rm samples.txt

#Now we obtain the final vcf versions and calculate the precision/recall
miniwalk mod -b ${OUT_DIR_FILE}_map2pan_short.bed -v ${OUT_DIR_FILE}_pan_short.vcf -g ${PAN_DIR}/Mtb_PRG.gfa -o ${OUT_DIR_FILE}_mod_short.vcf

rm ${OUT_DIR_FILE}_pan_short.vcf

miniwalk ref -v ${OUT_DIR_FILE}_mod_short.vcf -o ${OUT_DIR_FILE}_ref_short.vcf

rm ${OUT_DIR_FILE}_mod_short.vcf

#Now we get the SVs for the long-read assembly

minigraph -cxasm --call ${PAN_DIR}/Mtb_PRG.gfa ${OUT_DIR_FILE}.fasta > ${OUT_DIR_FILE}_map2pan.bed

paste ${PAN_DIR}/h37rv.bed ${OUT_DIR_FILE}_map2pan.bed | k8 /data/projects/punim1637/Aleix/mgutils.js merge - > ${OUT_DIR_FILE}_merged.bed

echo NC_000962.3 > samples.txt
echo $FILE >> samples.txt

k8 /data/projects/punim1637/Aleix/mgutils-es6.js merge2vcf -s samples.txt ${OUT_DIR_FILE}_merged.bed > ${OUT_DIR_FILE}_pan.vcf

rm ${OUT_DIR_FILE}_merged.bed
rm samples.txt

#Now we obtain the final vcf versions and calculate the precision/recall
miniwalk mod -b ${OUT_DIR_FILE}_map2pan.bed -v ${OUT_DIR_FILE}_pan.vcf -g ${PAN_DIR}/Mtb_PRG.gfa -o ${OUT_DIR_FILE}_mod.vcf

rm ${OUT_DIR_FILE}_pan.vcf

miniwalk ref -v ${OUT_DIR_FILE}_mod.vcf -o ${OUT_DIR_FILE}_ref.vcf

rm ${OUT_DIR_FILE}_mod.vcf

#We calculate the PR

miniwalk bench -t gns -c ${OUT_DIR_FILE}_ref_short.vcf -v ${OUT_DIR_FILE}.vcf -r data/repeats/TR_full_100_h37rv.txt -e ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta >> data/long_read/gold_assemblies/bench/short_read_bench_svimasm.tsv
miniwalk bench -t mns -c ${OUT_DIR_FILE}_short_manta.vcf -v ${OUT_DIR_FILE}.vcf -r data/repeats/TR_full_100_h37rv.txt -e ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta >> data/long_read/gold_assemblies/bench/short_read_bench_svimasm.tsv
miniwalk bench -t mnn -c ${OUT_DIR_FILE}_short_manta.vcf -v ${OUT_DIR_FILE}_ref.vcf -r data/repeats/TR_full_100_h37rv.txt -e ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta >> data/long_read/gold_assemblies/bench/short_read_bench.tsv
miniwalk bench -t gnn -c ${OUT_DIR_FILE}_ref_short.vcf -v ${OUT_DIR_FILE}_ref.vcf -r data/repeats/TR_full_100_h37rv.txt -e ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta >> data/long_read/gold_assemblies/bench/short_read_bench.tsv
miniwalk bench -t mon -c ${OUT_DIR_FILE}_short_manta.vcf -v ${OUT_DIR_FILE}_ref.vcf -r data/repeats/TR_full_100_h37rv.txt -e ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta >> data/long_read/gold_assemblies/bench/pr_manta_histogram.csv
miniwalk bench -t mos -c ${OUT_DIR_FILE}_short_manta.vcf -v ${OUT_DIR_FILE}.vcf -r data/repeats/TR_full_100_h37rv.txt -e ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta >> data/long_read/gold_assemblies/bench/pr_manta_histogram_svimasm.csv
miniwalk bench -t gon -c ${OUT_DIR_FILE}_ref_short.vcf -v ${OUT_DIR_FILE}_ref.vcf -r data/repeats/TR_full_100_h37rv.txt -e ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta >> data/long_read/gold_assemblies/bench/pr_graph_histogram.csv
miniwalk bench -t gos -c ${OUT_DIR_FILE}_ref_short.vcf -v ${OUT_DIR_FILE}.vcf -r data/repeats/TR_full_100_h37rv.txt -e ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta >> data/long_read/gold_assemblies/bench/pr_graph_histogram_svimasm.csv
miniwalk bench -t gns -c ${OUT_DIR_FILE}_ref.vcf -v ${OUT_DIR_FILE}.vcf -r data/repeats/TR_full_100_h37rv.txt -e ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta >> data/long_read/gold_assemblies/bench/svimasm_read_bench.tsv
miniwalk bench -t gos -c ${OUT_DIR_FILE}_ref.vcf -v ${OUT_DIR_FILE}.vcf -r data/repeats/TR_full_100_h37rv.txt -e ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta >> data/long_read/gold_assemblies/bench/pr_svimasm_histogram.csv

#VCFDIST
python3 cmd/aim1.4/replace_asterisk.py -v ${OUT_DIR_FILE}_ref_def.vcf -r ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta -o ${OUT_DIR_FILE}_ref_short2.vcf
gunzip ${OUT_DIR_FILE}_ref_short2.vcf.gz
sed -i 's/GT:GT0/GT/g;s/0:0/1/g' ${OUT_DIR_FILE}_ref_short2.vcf
sed -i 's/1\/1/1/g' ${OUT_DIR_FILE}.vcf

module load vcfdist
vcfdist ${OUT_DIR_FILE}_ref_short2.vcf ${OUT_DIR_FILE}.vcf ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta -l 100000 -t 4 -ct 0.001

cat query.tsv | awk '{if (length($4)>49){print $6,length($4),$8}}' >> data/long_read/gold_assemblies/bench/vcfdist_longass_precision.tsv
cat truth.tsv | awk '{if (length($4)>49){print $6,length($4),$8}}' >> data/long_read/gold_assemblies/bench/vcfdist_longass_recall.tsv;done

python3 cmd/aim1.4/replace_asterisk.py -v ${OUT_DIR_FILE}_ref.vcf -r ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta -o ${OUT_DIR_FILE}_ref2.vcf
gunzip ${OUT_DIR_FILE}_ref2.vcf.gz
sed -i 's/GT:GT0/GT/g;s/0:0/1/g' ${OUT_DIR_FILE}_ref2.vcf
sed -i 's/1\/1/1/g' ${OUT_DIR_FILE}.vcf

module load vcfdist
vcfdist ${OUT_DIR_FILE}_ref2.vcf ${OUT_DIR_FILE}.vcf ../../../Shared/Reference/MTB/NCBI/H37Rv.fasta -l 100000 -t 4 -ct 0.001

cat query.tsv | awk '{if (length($4)>49){print $6,length($4),$8}}' >> data/long_read/gold_assemblies/bench/vcfdist_longass_precision.tsv
cat truth.tsv | awk '{if (length($4)>49){print $6,length($4),$8}}' >> data/long_read/gold_assemblies/bench/vcfdist_longass_recall.tsv;done

done