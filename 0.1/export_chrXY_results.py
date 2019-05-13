
import sys
from hail import *
hc = HailContext()

sex = sys.argv[1] 
pipeline_type = sys.argv[2]
pipeline_number = sys.argv[3] 

print '#####################'
print '## Starting... ######'
print '## Sex: {}'.format(sex)
print '## Pipeline type: {}'.format(pipeline_type)
print '## Pipeline number: {:}'.format(pipeline_number)
print '#####################'
"""
for i in xrange(7):
    kt = hc.read_table('gs://ukb31063-mega-gwas/phenotype-pipelines/finngen/ukb31063.both_sexes.finngen.pipeline.{:}.kt'.format(i))
    if 'M13_ALGONEURO' in set(kt.columns):
        print 'FOUND IT - PIPELINE {:}'.format(i)
        print kt.columns
        kt_results = hc.read_table('gs://ukb31063-mega-gwas/results-tables/finngen/ukb31063.both_sexes.finngen.chrXY.pipeline.{:}.kt'.format(i))
        print 'n_results_variants: {:,}'.format(kt_results.count())

import sys
sys.exit()
"""

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

if pipeline_type == 'phesant':
    pipeline_path = 'gs://ukb31063-mega-gwas/phenotype-pipelines/final-phesant/ukb31063.{0}.phesant.pipeline.{1}.kt'.format(sex, pipeline_number)
else:
    pipeline_path = 'gs://ukb31063-mega-gwas/phenotype-pipelines/{2}/ukb31063.{0}.{2}.pipeline.{1}.kt'.format(sex, pipeline_number, pipeline_type)

results_chrXY_path = 'gs://ukb31063-mega-gwas/results-tables/{0}/ukb31063.{1}.{0}.chrXY.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)

kt_variants = hc.import_table('gs://ukb31063-mega-gwas/qc/ukb31063.gwas_variants.chrXY.tsv', key='v', types={'v': TVariant()})
kt_pipeline = hc.read_table(pipeline_path)
kt_results = hc.read_table(results_chrXY_path)
kt_results = kt_results.annotate('variant = Variant(v.contig.replace("XY", "X").replace("PAR1", "X"), v.start, v.ref, v.alt())')
kt_results = (kt_results.key_by('variant')
                        .order_by('variant')
                        .drop('v'))
kt_results = kt_results.join(kt_variants, how='inner')
kt_results = kt_results.cache()

if pipeline_type == 'phesant':
    kt_phenosummary = hc.read_table('gs://ukb31063-mega-gwas/phenotype-summaries/phesant/ukb31063.phesant_phenosummaries.{}.kt'.format(sex)).cache()
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='mix')
    print kt_phenosummary.schema
elif pipeline_type == 'icd10':
    kt_phenosummary = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/icd10/ukb31063.{0}.icd10.phenosummary.pipeline.{1:}.tsv'.format(sex, pipeline_number)).cache()
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='single')
elif pipeline_type == 'finngen':
    kt_phenosummary = hc.import_table('gs://ukb31063-mega-gwas/phenotype-summaries/finngen/ukb31063.{0}.finngen.phenosummary.pipeline.{1:}.tsv'.format(sex, pipeline_number)).cache()
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='single')

count = 1
for i, group in enumerate(phenotype_groups):
    for j, code in enumerate(group[1]):
        print 'Exporting {0} ({1:})...'.format(code, count)
        kt_export = kt_results.annotate(['n_complete_samples = va.results[{0:}].nCompleteSamples'.format(i),
                                         'AC = va.results[{0:}].AC'.format(i),
                                         'ytx = va.results[{0:}].ytx[{1:}]'.format(i, j),
                                         'beta = va.results[{0:}].beta[{1:}]'.format(i, j),
                                         'se = va.results[{0:}].se[{1:}]'.format(i, j),
                                         'tstat = va.results[{0:}].tstat[{1:}]'.format(i, j),
                                         'pval = va.results[{0:}].pval[{1:}]'.format(i, j)])
        kt_export = kt_export.annotate('AF = AC/(2.0 * n_complete_samples)')
        kt_export = kt_export.annotate(['minor_allele = if (AF <= 0.5) variant.alt() else variant.ref',
                                         'minor_AF = if (AF <= 0.5) AF else 1.0 - AF'])

        if pipeline_type == 'phesant':
            try:
                n_cases = int(kt_phenosummary.query('n_cases.filter(x => phenotype == "{}").collect()'.format(code))[0])
            except:
                value_counter = kt_pipeline.query('`{}`.map(x => str(x)).counter()'.format(code))
                if len(value_counter) <= 4:
                    min_category_count = min([int(x) for x in value_counter.values()])
                    kt_export = kt_export.annotate('expected_min_category_minor_AC = 2.0 * minor_AF * {:}.toInt()'.format(min_category_count))
                    kt_export = kt_export.annotate('low_confidence_variant = (expected_min_category_minor_AC < 25) || (minor_AF < 0.001)')
                else:
                    kt_export = kt_export.annotate('low_confidence_variant = minor_AF < 0.001')
            else:
                kt_export = kt_export.annotate('expected_case_minor_AC = 2.0 * minor_AF * {:}.toInt()'.format(n_cases))
                kt_export = kt_export.annotate('low_confidence_variant = (expected_case_minor_AC < 25) || (minor_AF < 0.001)')
        elif pipeline_type == 'icd10':
            n_cases = kt_phenosummary.query('n_cases.filter(x => code == "{}").collect()'.format(code))[0]
            kt_export = kt_export.annotate('expected_case_minor_AC = 2.0 * minor_AF * {0:}.toInt()'.format(n_cases))
            kt_export = kt_export.annotate('low_confidence_variant = (expected_case_minor_AC < 25) || (minor_AF < 0.001)')
        elif pipeline_type == 'finngen':
            n_cases = kt_phenosummary.query('n_cases.filter(x => code == "{}").collect()'.format(code))[0]
            kt_export = kt_export.annotate('expected_case_minor_AC = 2.0 * minor_AF * {0:}.toInt()'.format(n_cases))
            kt_export = kt_export.annotate('low_confidence_variant = (expected_case_minor_AC < 25) || (minor_AF < 0.001)')

        if 'expected_case_minor_AC' in kt_export.columns:
            kt_export = kt_export.select(['variant',
                                          'minor_allele',
                                          'minor_AF',
                                          'expected_case_minor_AC',
                                          'low_confidence_variant',
                                          'n_complete_samples',
                                          'AC',
                                          'ytx',
                                          'beta',
                                          'se',
                                          'tstat',
                                          'pval'])
        elif 'expected_min_category_minor_AC' in kt_export.columns:
            kt_export = kt_export.select(['variant',
                                          'minor_allele',
                                          'minor_AF',
                                          'expected_min_category_minor_AC',
                                          'low_confidence_variant',
                                          'n_complete_samples',
                                          'AC',
                                          'ytx',
                                          'beta',
                                          'se',
                                          'tstat',
                                          'pval'])
        else:
            kt_export = kt_export.select(['variant',
                                          'minor_allele',
                                          'minor_AF',
                                          'low_confidence_variant',
                                          'n_complete_samples',
                                          'AC',
                                          'ytx',
                                          'beta',
                                          'se',
                                          'tstat',
                                          'pval'])

        kt_export.export('gs://ukb31063-mega-gwas/chrXY-results-tsvs/{0}.imputed_v3.results.{1}.tsv.gz'.format(code, sex))
        count += 1

print '#####################'
print '## COMPLETED ######'
print '## Sex: {}'.format(sex)
print '## Pipeline type: {}'.format(pipeline_type)
print '## Pipeline number: {}'.format(pipeline_number)
print '#####################'

