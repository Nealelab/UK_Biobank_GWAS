
import sys
from hail import *
hc = HailContext()

sex = sys.argv[1] 
pipeline_number = sys.argv[2]

if sex == 'both_sexes':
    #pipelines = [str(i) for i in range(0, 8)]
    plural_sex = 'both_sexes' 
elif sex == 'female':
    #pipelines = [str(i) for i in range(18, 20)]
    plural_sex = 'females'
elif sex == 'male':
    #pipelines = [str(i) for i in range(12, 19)]
    plural_sex = 'males'     

print '#####################'
print '## Starting... ######'
print '## Sex: {}'.format(sex)
print '## Pipeline type: still_more_phesant'
print '## Pipeline number: {:}'.format(pipeline_number)
print '#####################'

with hadoop_read('gs://ukb31063-mega-gwas/still_more_phesant.to_export.{}.tsv'.format(sex)) as f:
    f.readline()
    to_export = set([x.strip().split('\t')[1] for x in f.readlines() if x.strip().split('\t')[0] == str(pipeline_number)])

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

kt_phenosummary = hc.read_table('gs://ukb31063-mega-gwas/phenotype-summaries/phesant/ukb31063.phesant_phenosummaries.{}.kt'.format(sex)).cache()

pipeline_path = 'gs://ukb31063-mega-gwas/phenotype-pipelines/still-more-phesant/ukb31063.{0}.still_more_phesant.pipeline.{1}.kt'.format(plural_sex, pipeline_number)
results_auto_path = 'gs://ukb31063-mega-gwas/results-tables/still-more-phesant/ukb31063.{0}.still_more_phesant.autosomes.pipeline.{1}.kt'.format(sex, pipeline_number)
results_chrX_path = 'gs://ukb31063-mega-gwas/results-tables/still-more-phesant/ukb31063.{0}.still_more_phesant.chrX.pipeline.{1}.kt'.format(sex, pipeline_number)

kt_pipeline = hc.read_table(pipeline_path)
kt_results_auto = hc.read_table(results_auto_path)
kt_results_chrX = hc.read_table(results_chrX_path)
kt_results = KeyTable.union(kt_results_auto, kt_results_chrX)
kt_results = kt_results.annotate('variant = Variant(v.contig.replace("^0", ""), v.start, v.ref, v.alt())')
kt_results = (kt_results.key_by('variant')
                        .order_by('variant')
                        .filter('isDefined(va.AF)')
                        .drop('v')).cache()

phenotype_groups = group_phenotypes(kt_pipeline, block_type='mix')

count = 1
for i, group in enumerate(phenotype_groups):
    for j, code in enumerate(group[1]):
        if code not in to_export:
            continue
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
        kt_export.export('gs://ukb31063-mega-gwas/export2-results-tsvs/{0}.gwas.imputed_v3.{1}.tsv.bgz'.format(code, sex))
        count += 1

print '#####################'
print '## COMPLETED ######'
print '## Sex: {}'.format(sex)
print '## Pipeline type: still_more_phesant'
print '## Pipeline number: {}'.format(pipeline_number)
print '#####################'
