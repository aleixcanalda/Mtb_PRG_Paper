# Code for "Discovery of DR-associated SVs in 41,134 Mtb isolates"

Description of the code:

* DR_matrix_slurm.cmd: slurm code to run DR_matrix.R in parallel.

* DR_matrix.R: R script for joining all 44k isolates' VCFs into one large matrix.

* DR.R: R script for finding DR-associated SVs and genes. Main script of this section.

* genotyping.cmd: pipeline for obtaining SVs from short-read Mtb data.

* hail.py: python script for running hail, a tool to do PCA of large datasets.

* hsdm_del_plot.py: python code to plot the hsdM deletion.

* RNA_Seq_DR.R R script doing differential expression analysis in an isolate with vs without hsdM.

* snippy.cmd: pipeline for obtaining VCf files from the 44k isolates and keeping their SNPs.

* tbprof.cmd: pipeline for running TB-Profiler to obtain drug resistance predictions for each isolate.
