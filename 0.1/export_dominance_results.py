
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

results_auto_path = 'gs://ukb31063-mega-gwas/dominance-results-tables/{0}/ukb31063.{1}.{0}.autosomes.dominance.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)
results_chrX_path = 'gs://ukb31063-mega-gwas/dominance-results-tables/{0}/ukb31063.{1}.{0}.chrX.dominance.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)

kt_pipeline = hc.read_table(pipeline_path)
kt_results_auto = hc.read_table(results_auto_path)
kt_results_chrX = hc.read_table(results_chrX_path)
kt_results = KeyTable.union(kt_results_auto, kt_results_chrX)
kt_results = kt_results.annotate('variant = Variant(v.contig.replace("^0", ""), v.start, v.ref, v.alt())')
kt_results = (kt_results.key_by('variant')
                        .order_by('variant')
                        .filter('isDefined(va.AF)')
                        .drop('v'))
kt_results = kt_results.repartition(116).cache()

if pipeline_type == 'phesant':
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='mix')
else:
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='single')

count = 1
for i, group in enumerate(phenotype_groups):
    for j, code in enumerate(group[1]):
        print 'Exporting {0} ({1:})...'.format(code, count)
        kt_export = kt_results.annotate(['dominance_AC = va.results[{0:}].AC'.format(i),
                                         'dominance_ytx = va.results[{0:}].ytx[{1:}]'.format(i, j),
                                         'dominance_beta = va.results[{0:}].beta[{1:}]'.format(i, j),
                                         'dominance_se = va.results[{0:}].se[{1:}]'.format(i, j),
                                         'dominance_tstat = va.results[{0:}].tstat[{1:}]'.format(i, j),
                                         'dominance_pval = va.results[{0:}].pval[{1:}]'.format(i, j)])
        kt_export = kt_export.select(['variant',
                                      'dominance_AC',
                                      'dominance_ytx',
                                      'dominance_beta',
                                      'dominance_se',
                                      'dominance_tstat',
                                      'dominance_pval'])
        kt_export.export('gs://ukb31063-mega-gwas/results-tsvs-dominance/{0}.dominance.gwas.imputed_v3.{1}.tsv.bgz'.format(code, sex))
        count += 1

print '#####################'
print '## COMPLETED ######'
print '## Sex: {}'.format(sex)
print '## Pipeline type: {}'.format(pipeline_type)
print '## Pipeline number: {}'.format(pipeline_number)
print '#####################'
