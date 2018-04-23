
from hail import *
hc = HailContext()

vds = hc.import_vcf('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.genotype.vcf.bgz')
vds.write('gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.genotype.vds', overwrite=True)
