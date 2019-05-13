
"""
from pprint import pprint
from hail import *
hc = HailContext()

kt = hc.read_table('gs://ukb31063-mega-gwas/phenotype-pipelines/icd10/ukb31063.both_sexes.icd10.pipeline.0.kt')
kt.export('gs://ubk31063-mega-gwas/dominance-testing/ukb31063.both_sexes.icd10.pipeline.0.tsv.bgz')
"""

import hail as hl

"""
ht = hl.import_table('gs://ukb31063-mega-gwas/dominance-testing/ukb31063.both_sexes.icd10.pipeline.0.tsv.bgz', key='s', impute=True, types={'s': hl.tstr})
ht.describe()
ht.write('gs://ukb31063-mega-gwas/dominance-testing/ukb31063.both_sexes.icd10.pipeline.0.ht', overwrite=True)
"""

"""
ht = hl.import_table('gs://phenotype_31063/ukb31063.gwas_covariates.both_sexes.tsv', key='s', impute=True, types={'s': hl.tstr})
ht.describe()
ht.write('gs://ukb31063-mega-gwas/dominance-testing/ukb31063.gwas_covariates.both_sexes.ht', overwrite=True)
"""

ht_phenotypes = hl.read_table('gs://ukb31063-mega-gwas/dominance-testing/ukb31063.both_sexes.icd10.pipeline.0.ht')
ht_covariates = hl.read_table('gs://ukb31063-mega-gwas/dominance-testing/ukb31063.gwas_covariates.both_sexes.ht')  
ht_variants = hl.read_table('gs://ukb31063-mega-gwas/dominance-testing/ukb31063.gwas_variants.ht')
ht_variants = ht_variants.filter(ht_variants['v'].startswith('22:'), keep=True)
ht_variants = ht_variants.annotate(**hl.parse_variant(ht_variants['v']))
ht_variants = ht_variants.key_by('locus', 'alleles')
ht_variants.describe()

mt = hl.import_bgen('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr22_v3.bgen',
                    sample_file='gs://phenotype_31063/ukb31063.autosomes.sample',
                    entry_fields=['GP'])
mt = mt.filter_rows(hl.is_defined(ht_variants[mt.locus, mt.alleles]), keep=True)

# annotate samples with the phenotype information.
mt = mt.annotate_cols(**ht_phenotypes[mt.s])
mt = mt.annotate_cols(**ht_covariates[mt.s])

mt = mt.annotate_rows(p = hl.agg.mean(mt.GP[1] + 2 * mt.GP[2]) / 2.0)
#mt = mt.annotate_rows(p = hl.agg.mean(mt.GT.n_alt_alleles()) / 2.0)
mt = mt.annotate_rows(c = mt.p / (mt.p - 1.0))
mt = mt.annotate_rows(one_over_c = 1 / mt.c)

mt = mt.select_entries(x = mt.GP[0] * mt.c + mt.GP[1] + mt.GP[2] * mt.one_over_c)

# run the regression with our edited dosage information.
phenotypes = [x for x in ht_phenotypes.row if x != 's']
mt = hl.linear_regression(y = [mt[x] for x in phenotypes],
                          x = mt.x, 
                          covariates=[1.0, mt.isFemale, mt.age, mt.age_squared, mt.age_isFemale, mt.age_squared_isFemale] + [mt['PC{:}'.format(i)] for i in range(1, 21)])

ht_results = mt.rows()
for i, phenotype in enumerate(phenotypes):
    ht_export = ht_results.annotate(n=ht_results['linreg'].n,
                                    y_transpose_x=ht_results['linreg'].y_transpose_x[i],
                                    beta=ht_results['linreg'].beta[i],
                                    standard_error=ht_results['linreg'].standard_error[i],
                                    t_stat=ht_results['linreg'].t_stat[i],
                                    p_value=ht_results['linreg'].p_value[i])
    ht_export = ht_export.select('n', 'y_transpose_x', 'beta', 'standard_error', 't_stat', 'p_value')
    ht_export.export('gs://ukb31063-mega-gwas/dominance-testing/{}.dominance.chr22.results.tsv.bgz'.format(phenotype))
    if i == 10:
        break

