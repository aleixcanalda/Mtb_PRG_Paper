# Run with python3.7 and jinja2 <3.1.0


import hail as hl

# Initialize Hail
hl.init()

# Define a custom reference genome
my_genome = hl.ReferenceGenome(
    name='Mycobacterium_tuberculosis',  # Name of the reference genome
    contigs=['NC_000962.3'],            # List of contigs
    lengths={'NC_000962.3': 4411532}    # Length of each contig
)

# Load VCF as MatrixTable
mt = hl.import_vcf("merged.vcf.gz",reference_genome=my_genome,force_bgz=True)

# Annotate with minor allele frequency (MAF)
mt = mt.annotate_rows(maf=hl.min(mt.af, 1 - mt.af))

# Filter by MAF (e.g., > 0.01)
filtered_mt = mt.filter_rows(mt.maf > 0.01)

# Perform PCA on genotypes (using GT)
scores, loadings, eigenvalues = hl.pca(filtered_mt.GT.n_alt_alleles(), k=10, compute_loadings=True)

# Write to a file
with open('scores_hail.tsv', 'w') as f:
    for score in scores:
        f.write(f"{score}\n")

# Save PCA loadings for variants
loadings.export('pca_loadings.tsv')
