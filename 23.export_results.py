
import sys
from hail import *
hc = HailContext()
from pprint import pprint

sex = sys.argv[1]             # {'both_sexes', 'female', 'male'}
pipeline_type = sys.argv[2]   # {'phesant', 'icd10', 'curated', 'finngen'}
try:
    pipeline_number = sys.argv[3] # {0, 1, 2, ...}
except:
    pipeline_number = ''

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError('Invalid sex value "{}" - must be one of {"both_sexes", "female", "male"}.'.format(sex))
if pipeline_type not in set(['phesant', 'icd10', 'curated', 'finngen']):
    raise ValueError('Invalid pipeline type "{}" - must be one of {"phesant", "icd10", "curated"}.'.format(pipeline_type))

print '#####################'
print '## Starting... ######'
print '## Sex: {}'.format(sex)
print '## Pipeline type: {}'.format(pipeline_type)
print '## Pipeline number: {}'.format(pipeline_number)
print '#####################'

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

if pipeline_type == 'curated':
    pipeline_path = 'gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/curated/ukb31063.{}.curated.kt'.format(sex)
    results_auto_path = 'gs://ukb31063-mega-gwas/hail-0.1/results-tables/curated/ukb31063.{}.curated.autosomes.kt'.format(sex)
    results_chrX_path = 'gs://ukb31063-mega-gwas/hail-0.1/results-tables/curated/ukb31063.{}.curated.chrX.kt'.format(sex)
else:
    pipeline_path = 'gs://ukb31063-mega-gwas/hail-0.1/phenotype-pipelines/{0}/ukb31063.{1}.{0}.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)
    results_auto_path = 'gs://ukb31063-mega-gwas/hail-0.1/results-tables/{0}/ukb31063.{1}.{0}.autosomes.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)
    results_chrX_path = 'gs://ukb31063-mega-gwas/hail-0.1/results-tables/{0}/ukb31063.{1}.{0}.chrX.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)

kt_pipeline = hc.read_table(pipeline_path)
kt_results_auto = hc.read_table(results_auto_path)
kt_results_chrX = hc.read_table(results_chrX_path)
kt_results = KeyTable.union(kt_results_auto, kt_results_chrX)
kt_results = kt_results.annotate('variant = Variant(v.contig.replace("^0", ""), v.start, v.ref, v.alt())')
kt_results = (kt_results.key_by('variant')
                        .order_by('variant')
                        .filter('isDefined(va.AF)')
                        .drop('v')).cache()

if pipeline_type == 'phesant':
    kt_phenosummary = hc.import_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/phesant/ukb31063.{0}.phenosummary.*.tsv'.format(sex)).cache()
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='mix')
elif pipeline_type == 'icd10':
    kt_phenosummary = hc.import_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/icd10/ukb31063.{0}.icd10.phenosummary.pipeline.{1:}.tsv'.format(sex, pipeline_number)).cache()
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='single')
elif pipeline_type == 'finngen':
    kt_phenosummary = hc.import_table('gs://ukb31063-mega-gwas/hail-0.1/phenotype-summaries/finngen/ukb31063.{0}.finngen.phenosummary.pipeline.{1:}.tsv'.format(sex, pipeline_number)).cache()
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='single')

count = 1
for i, group in enumerate(phenotype_groups):
    for j, code in enumerate(group[1]):
        print 'Exporting {0} ({1:})...'.format(code, count)
        kt_export = kt_results.annotate(['nCompleteSamples = va.results[{0:}].nCompleteSamples'.format(i),
                                         'AC = va.results[{0:}].AC'.format(i),
                                         'ytx = va.results[{0:}].ytx[{1:}]'.format(i, j),
                                         'beta = va.results[{0:}].beta[{1:}]'.format(i, j),
                                         'se = va.results[{0:}].se[{1:}]'.format(i, j),
                                         'tstat = va.results[{0:}].tstat[{1:}]'.format(i, j),
                                         'pval = va.results[{0:}].pval[{1:}]'.format(i, j)])
        if pipeline_type == 'phesant':
            try:
                n_cases = int(kt_phenosummary.query('`N.cases`.filter(x => FieldID == "{}").collect()'.format(code))[0])
            except TypeError:
                pass
            else:
                kt_export = kt_export.annotate('expected_case_minor_AC = if (va.AF <= 0.5) 2.0 * va.AF * {0:}.toInt() else 2.0 * (1.0 - va.AF) * {0:}.toInt()'.format(n_cases))
        elif pipeline_type == 'icd10':
            n_cases = kt_phenosummary.query('n_cases.filter(x => code == "{}").collect()'.format(code))[0]
            kt_export = kt_export.annotate('expected_case_minor_AC = if (va.AF <= 0.5) 2.0 * va.AF * {0:}.toInt() else 2.0 * (1.0 - va.AF) * {0:}.toInt()'.format(n_cases))
        elif pipeline_type == 'finngen':
            n_cases = kt_phenosummary.query('n_cases.filter(x => code == "{}").collect()'.format(code))[0]
            kt_export = kt_export.annotate('expected_case_minor_AC = if (va.AF <= 0.5) 2.0 * va.AF * {0:}.toInt() else 2.0 * (1.0 - va.AF) * {0:}.toInt()'.format(n_cases))
        kt_export = kt_export.drop('va')
        kt_export.export('gs://ukb31063-mega-gwas/hail-0.1/results-tsvs/{0}.imputed_v3.results.{1}.tsv.gz'.format(code, sex))
        count += 1

print '#####################'
print '## COMPLETED ######'
print '## Sex: {}'.format(sex)
print '## Pipeline type: {}'.format(pipeline_type)
print '## Pipeline number: {}'.format(pipeline_number)
print '#####################'
