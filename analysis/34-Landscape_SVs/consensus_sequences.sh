#!/bin/bash

# VCF file
VCF="/scratch/punim1637/Dunstan/SV/vcf/snps/merged_sv_re.vcf.gz"

# Output directory
OUTPUT_DIR="/scratch/punim1637/Dunstan/SV/vcf/snps/snp_sequences"

# Create output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Get the list of sample names from the VCF file
SAMPLES=$(bcftools query -l ${VCF})

# Loop through each sample and extract SNP sequences
for SAMPLE in ${SAMPLES}; do
    echo ">${SAMPLE}" > ${OUTPUT_DIR}/${SAMPLE}.fa
    bcftools query -s ${SAMPLE} -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' /scratch/punim1637/Dunstan/SV/vcf/snps/merged_sv_re.vcf.gz | awk -v OFS="" '{
    if ($5 == "./.") { 
        printf $3 
    } 
    else if (split($4, altArray, ",") == 2) {
        if ($5 == "1/1") {
            printf altArray[1]
        }
        else {
            printf altArray[2]
        }
    }
    else if (split($4, altArray, ",") == 3) {
	if ($5 == "1/1") {
		printf altArray[1]
	}
	else if ($5 == "2/2") {
		printf altArray[2]
	}
	else {
		printf altArray[3]
	}
    }
    else { 
        printf $4 
    } 
}'  >> ${OUTPUT_DIR}/${SAMPLE}.fa
done

echo "SNP sequences generated in ${OUTPUT_DIR}"

# Change to the output directory
cd ${OUTPUT_DIR}

mv GCA_012923765.1.fa ../
# Concatenate all individual FASTA files into one
cat ../GCA_012923765.1.fa *.fa > ../all_samples_snps_sv.fasta

# Return to the original directory
cd ..

echo "All samples concatenated into all_samples_snps.fasta"

