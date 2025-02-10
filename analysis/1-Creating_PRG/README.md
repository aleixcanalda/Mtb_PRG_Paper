# Code for "Construction of a Mtb pangenome reference graph"

Description of the code:

* basecalling.cmd: code used for basecalling the 25 newly ONT sequenced isolates.

* lineages.R: R code for visualizing the number of sub-lineages in our data.

* long_oxford.cmd: code used for downloading ONT reads, doing quality control, creating an assembly using flye and sending the assembly through QUAST.

* long_pb_clr.cmd: code used for downloading PB-CLR reads, doing quality control, creating an assembly using flye and sending the assembly through QUAST.

* long_pb_hifi.cmd: code used for downloading PB-Hifi reads, doing quality control, creating an assembly using flye and sending the assembly through QUAST

* long_tbprof.cmd: code used to determine the lineages of isolates with long-read data using TB-Profiler.

* ncbi_lineage_assign.cmd: code used to get the lineages of isolates with an NCBI assembly as well as a vcf files with SNPs from the NCBI assemblies.

* pangenome_creation.cmd: code used to create the Mtb-PRG using minigraph.

* quast_ncbi.cmd: code used to get the QUAST metrics for NCBI assemblies.

* quast_output_sed.sh: bash script to make the QUAST output readable in R.

* quast.R: R code used to filter out assemblies and make visualizations.
