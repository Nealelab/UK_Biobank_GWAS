
from hail import *
hc = HailContext()

vds = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.genotype.vds')
kt_variants = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.genotype_snp_qc.kt')
kt_samples = hc.read_table('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.kt')

kt_variants = kt_variants.filter('in_PCA')

vds = vds.filter_variants_table(kt_variants, keep=True)
vds = vds.filter_samples_table(kt_samples, keep=True)

print vds.count()
vds.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.pca_subset.vds', overwrite=True)
