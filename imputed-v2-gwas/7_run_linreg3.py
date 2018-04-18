from __future__ import print_function
from pprint import pprint
import sys
from hail import *

hc = HailContext()

APPLICATION = sys.argv[1]
PIPELINE = sys.argv[2]

WORKING_BUCKET = '...'
RESULTS_BUCKET = '...'
BGEN_BUCKET = '...'

BGEN_FILES = BGEN_BUCKET + 'imputed/ukb_imp_chr*_v2.bgen'
SAMPLE_FILE = WORKING_BUCKET + APPLICATION + '/' + APPLICATION + '.sample'

PIPELINE_TABLE = WORKING_BUCKET + APPLICATION + '/pipelines/pipeline' + str(PIPELINE) + '.kt'
RESULTS_TABLE = RESULTS_BUCKET + APPLICATION + '/results' + str(PIPELINE) + '.kt' 

GWAS_VARIANTS_VDS = WORKING_BUCKET + 'gwas_variants.vds'

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

vds = (
    hc
    .import_bgen(path=BGEN_FILES, sample_file=SAMPLE_FILE, tolerance=0.2, min_partitions=None)
    .annotate_variants_expr('va = {}')
    .annotate_variants_vds(hc.read(GWAS_VARIANTS_VDS).annotate_variants_expr('va.keep = true'), root='va')
    .filter_variants_expr('isDefined(va.keep)')
    .annotate_variants_expr('va = select(va, rsid)')
    .annotate_samples_table(hc.read_table(PIPELINE_TABLE), root='sa')
)

print('')
print('Application: ', APPLICATION)
print('Pipeline: ', PIPELINE)

first = True

for group in phenotype_groups:

    print('Group: ', group[0])
    print('Fields: ', group[1])

    if first:
        expr = 'va.results = [va.linreg]'
        first = False
    else:
        expr = 'va.results = va.results.append(va.linreg)'

    vds = (
        vds
        .linreg3(
            ys=['sa.`{}`'.format(field) for field in group[1]],
            covariates=[
                'sa.isFemale',
                'sa.PC1',
                'sa.PC2',
                'sa.PC3',
                'sa.PC4',
                'sa.PC5',
                'sa.PC6',
                'sa.PC7',
                'sa.PC8',
                'sa.PC9',
                'sa.PC10'
            ],
            use_dosages=True,
            variant_block_size=8
        )
        .annotate_variants_expr(expr)
    )

(
    vds
    .annotate_variants_expr('va = drop(va, linreg)')
    .variants_table()
    .write(RESULTS_TABLE, overwrite=True)
)
