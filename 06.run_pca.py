
from hail import *
hc = HailContext()

vds = hc.read('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.pca_subset.vds')
vds = vds.pca(scores='sa.pca_scores', loadings='va.pca_loadings', k=20)

vds.samples_table().write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.pca_scores.kt', overwrite=True)
vds.variants_table().write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.gwas_samples.pca_loadings.kt', overwrite=True)
