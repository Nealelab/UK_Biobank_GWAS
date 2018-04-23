from __future__ import print_function
import sys
from hail import *

hc = HailContext()

APPLICATION = sys.argv[1]   # "ukb" + 4-digit application number
PIPELINE = sys.argv[2]		# pipeline number, 0, 1, 2, etc.
EXPORT_TYPE = sys.argv[3]   # type of tsv to export, "assoc" or "ldsc"

WORKING_BUCKET = '...'
RESULTS_BUCKET = '...'
EXPORT_BUCKET = '...'

HM3_TABLE = 'hm3.r3.hg19.auto_bi_af.ukbb_gwas_qcpos.kt'

PIPELINE_TABLE = WORKING_BUCKET + APPLICATION + '/pipelines/pipeline' + str(PIPELINE) + '.kt'
RESULTS_TABLE = RESULTS_BUCKET + APPLICATION + '/results' + str(PIPELINE) + '.kt'
HM3_JOIN_TABLE = RESULTS_BUCKET + APPLICATION + '/hm3_joins/' + str(PIPELINE) + '.hm3.join.kt' 

phenotype_groups = []
for field in hc.read_table(PIPELINE_TABLE).columns:
    if field in set(['ID', 'isFemale'] + ['PC{}'.format(x) for x in xrange(1,11)]):
        continue
    prefix = field.split('_', 1)[0]
    try:
        idx = [x[0] for x in phenotype_groups].index(prefix)
    except ValueError:
        phenotype_groups.append([prefix, [field]])
    else:
        phenotype_groups[idx][1].append(field)

kt = hc.read_table(RESULTS_TABLE).annotate('variant = Variant(v.contig.replace("^0", ""), v.start, v.ref, v.alt())').key_by('variant')

if EXPORT_TYPE == 'ldsc':
	kt_hm3 = (
		hc
		.read_table(HM3_TABLE)
		.annotate('variant = Variant(v_ukb.contig.replace("^0", ""), v_ukb.start, v_ukb.ref, v_ukb.alt())')
		.key_by('variant')
		.select('variant')
	)
	(
		kt
		.join(kt_hm3, how='inner')
		.write(HM3_JOIN_TABLE, overwrite=True)
	)
	kt = hc.read_table(HM3_JOIN_TABLE)

n_variants = kt.count()

print('')
print('Application: ', APPLICATION)
print('Pipeline: ', PIPELINE)
print('Export: ', EXPORT_TYPE)
print('nVariants: ', '{:,}'.format(n_variants))
print('Phenotypes: ', phenotype_groups)

for i, group in enumerate(phenotype_groups):
	for j, element in enumerate(group[1]):
		if EXPORT_TYPE == 'assoc':
			(
				kt
				.annotate([
					'rsid = va.rsid',
					'nCompleteSamples = va.results[{}].nCompleteSamples'.format(i),
					'AC = va.results[{}].AC'.format(i),
					'ytx = va.results[{0}].ytx[{1}]'.format(i, j),
					'beta = va.results[{0}].beta[{1}]'.format(i, j),
					'se = va.results[{0}].se[{1}]'.format(i, j),
					'tstat = va.results[{0}].tstat[{1}]'.format(i, j),
					'pval = va.results[{0}].pval[{1}]'.format(i, j)
				])
				.select([
					'variant',
					'rsid',
					'nCompleteSamples',
					'AC',
					'ytx',
					'beta',
					'se',
					'tstat',
					'pval'
				])
				.export(EXPORT_BUCKET + APPLICATION + '/associations/' + element + '.assoc.tsv.gz')
			)
		elif EXPORT_TYPE == 'ldsc':
			(
	            kt
	            .annotate([
	                'SNP = va.rsid',
	                'A1 = v.alt()',
	                'A2 = v.ref',
	                'N = va.results[{}].nCompleteSamples'.format(i),
	                'Z = va.results[{0}].tstat[{1}]'.format(i, j)
	            ])
	            .select([
	                'SNP',
	                'A1',
	                'A2',
	                'Z',
	                'N'
	            ])
	            .export(EXPORT_BUCKET + APPLICATION + '/ldsc/' + element + '.ldsc.tsv.gz')
	        )

