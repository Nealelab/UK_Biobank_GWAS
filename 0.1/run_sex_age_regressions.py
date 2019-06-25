
import sys
from hail import *

sex = sys.argv[1]           
contig = sys.argv[2] 

if contig == 'autosomes':
    hc = HailContext(min_block_size=512)
elif contig == 'chrX':
    hc = HailContext()

print '#####################'
print '## Starting... ######'
print '## Sex: {}'.format(sex)
print '## Contig: {}'.format(contig)
print '#####################'

vds_variants = hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.gwas_variants.{}.with_qc_annotations.vds'.format(contig))
kt_covariates = hc.read_table('gs://ukb31063-mega-gwas/qc/ukb31063.gwas_covariates.{}.kt'.format(sex))

if contig == 'autosomes':
    import_expr = '{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'
elif contig == 'chrX':
    import_expr = 'X'

vds = (hc.import_bgen('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{}_v3.bgen'.format(import_expr),
                      sample_file='gs://phenotype_31063/ukb31063.{}.sample'.format(contig),
                      tolerance=0.2)
         .annotate_variants_vds(vds_variants, expr='va.AF = vds.qc.AF, va.info = vds.info')
         .filter_variants_expr('isDefined(va.AF)', keep=True)
         .annotate_samples_table(kt_covariates, root='sa.covariates'))

if sex == 'both_sexes':
    vds = vds.linreg3(ys=['sa.covariates.isFemale'],
                      covariates=['sa.covariates.age', 'sa.covariates.age_squared'] + ['sa.covariates.PC{:}'.format(i) for i in xrange(1, 21)],
                      use_dosages=True,
                      variant_block_size=8)
    vds = vds.annotate_variants_expr('va.isFemale_results = va.linreg')

age_covariates = ['sa.covariates.PC{:}'.format(i) for i in xrange(1, 21)]
if sex == 'both_sexes':
    age_covariates += ['sa.covariates.isFemale']

vds = vds.linreg3(ys=['sa.covariates.age'],
                  covariates=age_covariates,
                  use_dosages=True,
                  variant_block_size=8)
vds = vds.annotate_variants_expr('va.age_results = va.linreg')

kt_results = (vds.annotate_variants_expr('va = drop(va, linreg)')
                 .variants_table())

kt_results.write('gs://ukb31063-mega-gwas/results-tables/ukb31063.isFemale_age.{0}.{1}.kt'.format(sex, contig))
