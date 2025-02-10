#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --job-name=pan
#SBATCH --output=stdout_pan.txt
#SBATCH --error=error_pan.txt
#SBATCH --time=04:00:00

module load foss
module load minigraph

FILES='/data/projects/punim1637/Aleix/SV/SV_analysis/final_isolates.txt'

#samtools index /data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa

while read i; do
if [[ -f /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/Mtb_PRG.gfa ]]; then
	if [[ -f /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/${i}_polished_filtered.fasta ]]; then
		minigraph -cxggs -L50 -t12 -l 10000 -d 5000 /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/Mtb_PRG.gfa /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/${i}_polished_filtered.fasta > /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/out.gfa
		mv /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/out.gfa /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/Mtb_PRG.gfa
	else
		gzip /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies_ncbi/${i}*.fna
		minigraph -L50 -cxggs -t12 -l 10000 -d 5000 /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/Mtb_PRG.gfa /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies_ncbi/${i}*.fna.gz > /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/out.gfa
	mv /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/out.gfa /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/Mtb_PRG.gfa
	fi
else
	minigraph -cxggs -L50 -t12 -l 10000 -d 5000 /data/gpfs/projects/punim1637/Shared/Data/Pathogen/Thai/Tutorials/Tutorial_04_Useful_Resources/TB.REF.GENOME.H37RV.ONELINE.fa /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/${i}_polished_filtered.fasta > /data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/assemblies/Mtb_PRG.gfa
fi
done < $FILES

