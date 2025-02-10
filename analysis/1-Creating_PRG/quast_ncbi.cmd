#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=quast_ncbi
#SBATCH --output=stdout_quast.txt
#SBATCH --error=error_quast.txt
#SBATCH --time=05:00:00
#SBATCH --array=0-493
#SBATCH --mem=15GB

module load GCC/11.3.0
module load OpenMPI/4.1.4
module load QUAST/5.2.0
module load foss/2022a
module load R/4.2.2
module load SAMtools/1.16.1
module load BEDTools/2.30.0
module load Meryl/1.4
module load GCCcore/11.3.0
module load minimap2/2.24
module load BWA/0.7.17

source /apps/easybuild-2022/easybuild/software/Core/Anaconda3/2022.10/etc/profile.d/conda.sh

export MERQURY=/data/projects/punim1637/Aleix/merqury

FILES=($(<$1))
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
OUT_DIR='/data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies_ncbi'
OUT_DIR_FILE=/scratch/punim1637/Dunstan/SV/fastq/${FILE}
DB_DIR='/data/projects/punim1637/Aleix/SV/SV_analysis/data'
FASTQ='/scratch/punim1637/Dunstan/SV/fastq'

quast ${OUT_DIR}/${FILE}*.fna.gz \
        -r /data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa \
        -o $FASTQ/${FILE}_quast -g /data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/Mycobacterium_tuberculosis_H37Rv_gff_v4.gff

mv $FASTQ/${FILE}_quast/contigs_reports/contigs_report_*.stdout /scratch/punim1637/Dunstan/SV/vcf/${FILE}_misassemblies_report.stdout
mv $FASTQ/${FILE}_quast/report.tsv ${OUT_DIR}/quast/${FILE}_quast_general_report.tsv

rm -r  $FASTQ/${FILE}_quast
