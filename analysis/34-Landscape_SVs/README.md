# Code for "Landscape of structural variation across the Mtb genome and phylogeny" and "Pangenomic insights into structural variation in the ESX loci of Mtb"

Description of the code:

* ancestral.cmd: pipeline used for calling SVs, genotyping against the ancestral sequence MTBC0.

* consensus_sequences: python code to extract a consensus sequence in a mult-sample SNP VCF file, useful for creating a phylogeny.

* genotype_sv_longass.cmd: pipeline used for genotyping SVs against the reference sequence H37Rv.

* long_vcf.cmd: calling SNPs from long-read data.

* phylo_tree.cmd: code used to create the phylogenetic tree.

* sv_longass_characterize.R: main file of this section. Used for processing all vcf files, doing various statistical analyses and visualisations.