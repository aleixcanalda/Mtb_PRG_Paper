# Code for "Genotyping SVs using the Mtb-PRG and miniwalk is more precise than using a traditional short-read SV caller"

Description of the code:

* benchmark.cmd: pipeline used for benchmarking the called SVs against the truth SVs using miniwalk bench and vcfdist for both long- and short-read data.

* pr_bench.R: code used for visualizing the benchmark results.

* replace_asterisk.py: python code for replacing the asterisk in the final miniwalk VCF file, to make the file appropriate for vcfidst use.