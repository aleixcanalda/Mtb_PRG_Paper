import pysam
import os
import pandas as pd
import sys
import gzip
import re
import argparse
import re
from Bio.Seq import Seq
import warnings
from Bio import Align

aligner = Align.PairwiseAligner()

warnings.simplefilter(action='ignore', category=FutureWarning)

parser = argparse.ArgumentParser(description="This script takes a bed file outputted from minigraph -xasm --call and the vcf from merge2vcf to modify the vcf file to show the SV type and exact positions. Moreover, the gfa file of the pangenome will also be necessary to determine the exact SV position and length.")

parser.add_argument('-v', '--vcf',
                    dest = "vcf",
                    action = "store",
                    required = True,
                    help = "The vcf file output from merge2vcf.")

parser.add_argument('-o', '--output',
                    dest = "out",
                    action = "store",
                    required = True,
                    help = "The new vcf file with the SVs.")

parser.add_argument('-r', '--reference',
                    dest = "ref",
                    action = "store",
                    required = True,
                    help = "The map2pan bed file for finding NAs.")

options = parser.parse_args()

def replace_asterisk_with_ref(vcf_file, ref_genome, output_vcf):
    # Open the VCF file and the reference genome
    vcf = pysam.VariantFile(vcf_file)
    ref = pysam.FastaFile(ref_genome)

    # Check if contigs are defined in the VCF header, add if missing
    if not vcf.header.contigs:
        for seq in ref.references:
            vcf.header.contigs.add(seq, length=ref.get_reference_length(seq))

    # Add necessary INFO fields if they are missing
    if 'END' not in vcf.header.info:
        vcf.header.info.add('END', 1, 'String', 'End position of the variant')
    if 'AN' not in vcf.header.info:
        vcf.header.info.add('AN', 1, 'String', 'Allele number')
    if 'NS' not in vcf.header.info:
        vcf.header.info.add('NS', 1, 'String', 'Number of samples')
    if 'NA' not in vcf.header.info:
        vcf.header.info.add('NA', 1, 'String', 'Number of alleles')
    if 'ALEN' not in vcf.header.info:
        vcf.header.info.add('ALEN', '.', 'Integer', 'Lengths of alleles')
    if 'AC' not in vcf.header.info:
        vcf.header.info.add('AC', '.', 'Integer', 'Allele count')
    if 'VS' not in vcf.header.info:
        vcf.header.info.add('VS', 1, 'String', 'Variant source')
    if 'VE' not in vcf.header.info:
        vcf.header.info.add('VE', 1, 'String', 'Variant end')
    if 'AWALK' not in vcf.header.info:
        vcf.header.info.add('AWALK', '.', 'String', 'Allele walk')

    # Create a new VCF file for the output
    output = pysam.VariantFile(output_vcf, 'w', header=vcf.header)

    # Iterate over each record in the VCF file
    for record in vcf:
        ref_base = ref.fetch(record.chrom, record.pos - 1, record.pos)
        
        if record.ref == '*':
            record.ref = ref_base
        record.alts = tuple(alt if alt != '*' else ref_base for alt in record.alts)
        
        output.write(record)
    
    # Close the files
    vcf.close()
    ref.close()
    output.close()

    # Index the output VCF file
    pysam.tabix_index(output_vcf, preset='vcf')

# Replace the placeholders with your actual file paths
vcf_file = '/data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/gold_assemblies/mtb.N0004_ref2.vcf'
ref_genome = '/data/projects/punim1637/Shared/Reference/MTB/NCBI/H37Rv.fasta'
output_vcf = '/data/projects/punim1637/Aleix/SV/SV_analysis/data/long_read/gold_assemblies/mtb.N0004_ref22.vcf'

replace_asterisk_with_ref(options.vcf,options.ref,options.out)

