import sys
from hail import *

sex = sys.argv[1]             # {'both_sexes', 'female', 'male'}
pipeline_type = sys.argv[2]   # {'phesant', 'icd10', 'finngen'}
contig = sys.argv[3]          # {'autosomes', 'chrX'}
try:
    pipeline_number = sys.argv[4] # {0, 1, 2, ...}
except:
    pipeline_number = ''

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError('Invalid sex value "{}" - must be one of {"both_sexes", "female", "male"}.'.format(sex))
if pipeline_type not in set(['phesant', 'icd10', 'finngen']):
    raise ValueError('Invalid pipeline type "{}" - must be one of {"phesant", "icd10", "finngen"}.'.format(pipeline_type))
if contig not in set(['autosomes', 'chrX']):
    raise ValueError('Invalid contig "{}" - must be one of {"autosomes", "chrX"}.'.format(contig))

if contig == 'autosomes':
    hc = HailContext(min_block_size=512)
elif contig == 'chrX':
    hc = HailContext()

print '#####################'
print '## Starting... ######'
print '## Sex: {}'.format(sex)
print '## Pipeline type: {}'.format(pipeline_type)
print '## Pipeline number: {}'.format(pipeline_number)
print '## Contig: {}'.format(contig)
print '#####################'

vds_variants = hc.read('gs://ukb31063-mega-gwas/qc/ukb31063.gwas_variants.{}.with_qc_annotations.vds'.format(contig))
kt_covariates = hc.read_table('gs://ukb31063-mega-gwas/qc/ukb31063.gwas_covariates.{}.kt'.format(sex))

def group_phenotypes(kt, block_type):
    codes = [x for x in kt.columns if x != 's']
    if block_type == 'single':
        phenotype_groups = [['single_block', codes]]
    elif block_type == 'individual':
        phenotype_groups = [[c, [c]] for c in codes]
    elif block_type == 'mix':
        phenotype_groups = []
        for code in codes:
            prefix = code.split('_')[0]
            try:
                idx = [x[0] for x in phenotype_groups].index(prefix)
            except ValueError:
                phenotype_groups.append([prefix, [code]])
            else:
                phenotype_groups[idx][1].append(code)
    return phenotype_groups

if pipeline_type in set(['icd10', 'finngen']):
    phenotype_groups = group_phenotypes(kt_phenotype, block_type='single')
elif pipeline_type == 'phesant':
    phenotype_groups = group_phenotypes(kt_phenotype, block_type='mix')
elif pipeline_type == 'curated':
    phenotype_groups = group_phenotypes(kt_phenotype, block_type='individual')  

if contig == 'autosomes':
    import_expr = '{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'
elif contig == 'chrX':
    import_expr = 'X'

if sex == 'both_sexes':
    covariate_expr = ['sa.covariates.age',
                      'sa.covariates.age_squared',
                      'sa.covariates.isFemale',
                      'sa.covariates.age_isFemale',
                      'sa.covariates.age_squared_isFemale'] + \
                     ['sa.covariates.PC{:}'.format(i) for i in xrange(1, 21)]
elif sex == 'female' or sex == 'male':
    covariate_expr = ['sa.covariates.age', 
                      'sa.covariates.age_squared'] + \
                     ['sa.covariates.PC{:}'.format(i) for i in xrange(1, 21)]

vds = (hc.import_bgen('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{}_v3.bgen'.format(import_expr),
                      sample_file='gs://phenotype_31063/ukb31063.{}.sample'.format(contig),
                      tolerance=0.2)
         .annotate_variants_vds(vds_variants, expr='va.AF = vds.qc.AF, va.info = vds.info')
         .filter_variants_expr('isDefined(va.AF)', keep=True)
         .annotate_samples_table(kt_covariates, root='sa.covariates')
         .annotate_samples_table(kt_phenotype, root='sa.phenotypes'))

first = True
for group, group_codes in phenotype_groups:
    vds = vds.linreg3(ys=['sa.phenotypes.`{}`'.format(x) for x in group_codes],
                      covariates=covariate_expr,
                      use_dosages=True,
                      variant_block_size=8)
    if first:
        vds = vds.annotate_variants_expr('va.results = [va.linreg]')
        first = False
    else:
        vds = vds.annotate_variants_expr('va.results = va.results.append(va.linreg)')

kt_results = (vds.annotate_variants_expr('va = drop(va, linreg)')
                 .variants_table())

path = 'gs://ukb31063-mega-gwas/dominance-results-tables/{0}/ukb31063.{1}.{0}.{2}.dominance.pipeline.{3}.kt'.format(pipeline_type, sex, contig, pipeline_number)
    
kt_results.write(path, overwrite=True)
kt = hc.read_table(path)
print 'nVariants: {:,}'.format(kt.count())

print '#####################'
print '## COMPLETE #########'
print '## Sex: {}'.format(sex)
print '## Pipeline type: {}'.format(pipeline_type)
print '## Pipeline number: {}'.format(pipeline_number)
print '## Contig: {}'.format(contig)
print '#####################'
