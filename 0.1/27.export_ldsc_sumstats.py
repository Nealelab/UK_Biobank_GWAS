
# example submission:
# $ cluster submit mycluster 26.export_ldsc_sumstats.py --args "both_sexes phesant 1"

# need enough workers to cache in memory for optimal export speed (used 15 for export of GWAS results with 13M variants, 
# could probably get away with 2-4 workers for 1M ldsc variants)

import sys
from hail import *
hc = HailContext()

sex = sys.argv[1]                 # {'both_sexes', 'female', 'male'}
pipeline_type = sys.argv[2]       # {'phesant', 'icd10', 'finngen'}
try:
    pipeline_number = sys.argv[3] # {0, 1, 2, ...}
except:
    pipeline_number = ''

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError('Invalid sex value "{}" - must be one of {"both_sexes", "female", "male"}.'.format(sex))
if pipeline_type not in set(['phesant', 'icd10', 'finngen']):
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

pipeline_path = 'gs://ukb31063-mega-gwas/phenotype-pipelines/{0}/ukb31063.{1}.{0}.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)
results_auto_path = 'gs://ukb31063-mega-gwas/results-tables/{0}/ukb31063.{1}.{0}.autosomes.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)
results_chrX_path = 'gs://ukb31063-mega-gwas/results-tables/{0}/ukb31063.{1}.{0}.chrX.pipeline.{2}.kt'.format(pipeline_type, sex, pipeline_number)

kt_pipeline = hc.read_table(pipeline_path)
kt_results_auto = hc.read_table(results_auto_path)
kt_results_chrX = hc.read_table(results_chrX_path)
kt_results = KeyTable.union(kt_results_auto, kt_results_chrX)
kt_results = kt_results.annotate('variant = Variant(v.contig.replace("^0", ""), v.start, v.ref, v.alt())')
kt_results = (kt_results.key_by('variant')
                        .order_by('variant')
                        .filter('isDefined(va.AF)')
                        .drop('v'))

kt_hm3 = hc.read_table('gs://ukb31063-mega-gwas/ldsc/ld_ref_panel/hm3.r3.b37.auto_bi_af.ukbb_gwas_qcpos.no_mhc.kt')
kt_hm3 = kt_hm3.select('v').key_by('v')
kt_join = kt_results.join(kt_hm3, how='inner').cache()

if pipeline_type == 'phesant':
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='mix')
elif pipeline_type == 'icd10':
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='single')
elif pipeline_type == 'finngen':
    phenotype_groups = group_phenotypes(kt_pipeline, block_type='single')

count = 1
for i, group in enumerate(phenotype_groups):
    for j, code in enumerate(group[1]):
        print 'Exporting LDSC sumstats for trait {0} ({1:})...'.format(code, count)
        kt_export = kt_join.annotate(['SNP = va.rsid',
                                      'A1 = variant.alt()',
                                      'A2 = variant.ref',
                                      'N = va.results[{:}].nCompleteSamples'.format(i),
                                      'Z = va.results[{0:}].tstat[{1:}]'.format(i, j)])
        kt_export = kt_export.select(['SNP', 'A1', 'A2', 'N', 'Z'])
        kt_export.export('gs://ukb31063-mega-gwas/ldsc-sumstats-tsvs/{0}.imputed_v3.ldsc.{1}.tsv.gz'.format(code, sex))
        count += 1

print '#####################'
print '## COMPLETED ######'
print '## Sex: {}'.format(sex)
print '## Pipeline type: {}'.format(pipeline_type)
print '## Pipeline number: {}'.format(pipeline_number)
print '#####################'
