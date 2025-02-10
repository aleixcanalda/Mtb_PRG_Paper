#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=pan
#SBATCH --output=stdout_pan.txt
#SBATCH --error=error_pan.txt
#SBATCH --time=01:00:00
#SBATCH --array=0-1174
#SBATCH --mem=10GB

module load foss/2022a
module load Anaconda3
module load Python/3.10.4
module load minimap2/2.24
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load GCCcore/11.3.0
module load foss
module load minigraph
module load k8/1.0
module load MUMmer/4.0.0rc1

source /apps/easybuild-2022/easybuild/software/Core/Anaconda3/2022.10/etc/profile.d/conda.sh

FILES=($(<$1))
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
OUT_DIR=/data/projects/punim1637/Aleix/SV/SV_analysis/${FILE}
OUT_DIR_FILE=/data/projects/punim1637/Aleix/SV/SV_analysis/${FILE}
OUT_DIR_NCBI=/data/projects/punim1637/Aleix/SV/SV_analysis/${FILE}
OUT_DIR_FILE_NCBI=/data/projects/punim1637/Aleix/SV/SV_analysis/${FILE}
PAN_DIR=/data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies
SCRIPT=/data/projects/punim1637/Aleix/SV/SV_analysis/cmd/aim1.4
REF=/data/projects/punim1637/Shared/Reference/MTB/NCBI/H37Rv.fasta
READ_DIR=/scratch/punim1637/Dunstan/SV/fastq/${FILE}

virtualenv ~/venvs/venv-3.10.4
source ~/venvs/venv-3.10.4/bin/activate

if [[ ${FILE:0:1} == 'G' ]]; then
	gunzip data/long_read/assemblies_ncbi/${FILE}*.fna.gz
	minigraph -cxasm --call ${PAN_DIR}/Mtb_PRG.gfa data/long_read/assemblies_ncbi/${FILE}*.fna > ${OUT_DIR_FILE_NCBI}_map2pan.bed
	gzip data/long_read/assemblies_ncbi/${FILE}*.fna
	paste ${PAN_DIR}/MTBC0_v1.1.bed ${OUT_DIR_FILE_NCBI}_map2pan.bed | k8 /data/projects/punim1637/Aleix/mgutils.js merge - > ${OUT_DIR_FILE_NCBI}_merged.bed

	echo NC_000962.3 > samples_${FILE}.txt
	echo $FILE >> samples_${FILE}.txt

	k8 /data/projects/punim1637/Aleix/mgutils-es6.js merge2vcf -s samples_${FILE}.txt ${OUT_DIR_FILE_NCBI}_merged.bed > ${OUT_DIR_FILE_NCBI}_pan.vcf
	
	rm ${OUT_DIR_FILE_NCBI}_merged.bed
	rm samples_${FILE}.txt
	
	miniwalk mod -v ${OUT_DIR_FILE_NCBI}_pan.vcf -g ${PAN_DIR}/Mtb_PRG.gfa -o ${OUT_DIR_FILE_NCBI}_mod.vcf -b ${OUT_DIR_FILE_NCBI}_map2pan.bed
	
	rm ${OUT_DIR_FILE_NCBI}_pan.vcf
	rm ${OUT_DIR_FILE_NCBI}_map2pan.bed
	
	miniwalk ref -v ${OUT_DIR_FILE_NCBI}_mod.vcf -o ${OUT_DIR_NCBI}_ref.vcf
	
	rm ${OUT_DIR_FILE_NCBI}_mod.vcf

	echo "${FILE}	$(grep -v '^#' ${OUT_DIR_NCBI}_ref.vcf | wc -l)" >> ancestral_svs.tsv

	rm ${OUT_DIR_NCBI}_ref.vcf*

	

else
	minigraph -cxasm --call ${PAN_DIR}/Mtb_PRG.gfa data/long_read/assemblies/${FILE}_polished_filtered.fasta > ${OUT_DIR_FILE}_map2pan.bed

	paste ${PAN_DIR}/MTBC0_v1.1.bed ${OUT_DIR_FILE}_map2pan.bed | k8 /data/projects/punim1637/Aleix/mgutils.js merge - > ${OUT_DIR_FILE}_merged.bed
	
	echo NC_000962.3 > samples_${FILE}.txt
	echo $FILE >> samples_${FILE}.txt
	
	k8 /data/projects/punim1637/Aleix/mgutils-es6.js merge2vcf -s samples_${FILE}.txt ${OUT_DIR_FILE}_merged.bed > ${OUT_DIR_FILE}_pan.vcf
	
	rm ${OUT_DIR_FILE}_merged.bed
	rm samples_${FILE}.txt
	
	#Now we obtain the final vcf versions
	miniwalk mod -v ${OUT_DIR_FILE}_pan.vcf -g ${PAN_DIR}/Mtb_PRG.gfa -o ${OUT_DIR_FILE}_mod.vcf -b ${OUT_DIR_FILE}_map2pan.bed
	
	rm ${OUT_DIR_FILE}_pan.vcf
	rm ${OUT_DIR_FILE}_map2pan.bed
	
	miniwalk ref -v ${OUT_DIR_FILE}_mod.vcf -o ${OUT_DIR}_ref.vcf
	
	rm ${OUT_DIR_FILE}_mod.vcf

	echo "${FILE}	$(grep -v '^#' ${OUT_DIR}_ref.vcf | wc -l)" >> ancestral_svs.tsv

        rm ${OUT_DIR}_ref.vcf*
fi
