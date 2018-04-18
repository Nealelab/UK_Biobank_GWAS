
import hail as hl

mt = hl.read_matrix_table('gs://ukb31063-mega-gwas/qc/ukb31063.genotype.mt')
hl.export_vcf(mt, 'gs://ukb31063-mega-gwas/hail-0.1/qc/ukb31063.genotype.vcf.bgz')
