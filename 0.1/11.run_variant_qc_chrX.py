
from hail import *
hc = HailContext()

kt_samples = hc.read_table('gs://ukb31063-mega-gwas/qc/ukb31063.gwas_samples.females.kt')

vds = hc.import_bgen('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chrX_v3.bgen',
                     sample_file='gs://phenotype_31063/ukb31063.imputed_v3.chrX.sample')
vds = vds.filter_samples_table(kt_samples, keep=True)

print 'nSamples: {:,}'.format(vds.num_samples)

vds = vds.variant_qc()
vds = vds.drop_samples()

vds.write('gs://ukb31063-mega-gwas/qc/ukb31063.imputed_v3.variant_qc.chrX.vds', overwrite=True)
