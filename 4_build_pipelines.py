from __future__ import print_function
from pprint import pprint
import sys
from hail import *

hc = HailContext()

APPLICATION = sys.argv[1]

BUCKET = '...'

PREFIX = BUCKET + APPLICATION + '/' + APPLICATION + '_'

SAMPLE_QC_TABLE = PREFIX + 'qc.kt'
PHESANT_FILE = PREFIX + 'output.tsv'

kt = (
    hc
    .import_table(PHESANT_FILE, impute=True, types={'userId': TString()})
    .rename({'userId': 'ID'})
    .key_by('ID')
    .join(hc.read_table(SAMPLE_QC_TABLE), how='inner')
    .cache()
)

phenotype_groups = []
for field in kt.columns:
    if field in set(['ID', 'isFemale', 'age', 'sex'] + ['PC{}'.format(x) for x in xrange(1,11)]):
        continue
    prefix = field.split('_', 1)[0]
    try:
        idx = [x[0] for x in phenotype_groups].index(prefix)
    except ValueError:
        phenotype_groups.append([prefix, [field]])
    else:
        phenotype_groups[idx][1].append(field)

max_phenotypes = 110
max_linreg = 37
linreg_calls = 0
effective_phenotypes_run = 0
pipelines = []
next_pipeline = []

for group, fields in phenotype_groups:
    blocks = [fields[i:(i+max_phenotypes)] for i in xrange(0, len(fields), max_phenotypes)]
    for block in blocks:
        effective_phenotypes_run += 2 + len(block)
        linreg_calls += 1
        if (effective_phenotypes_run > max_phenotypes) or (linreg_calls > max_linreg):
            pipelines.append(next_pipeline)
            next_pipeline = []
            effective_phenotypes_run = 2 + len(block)
            linreg_calls = 1
        next_pipeline.append(block)
pipelines.append(next_pipeline)

for i in xrange(len(pipelines)):

    fields = [field for linreg_call in pipelines[i] for field in linreg_call]

    (
        kt
        .select(['ID', 'isFemale'] + ['PC{}'.format(x) for x in xrange(1,11)] + fields)
        .write(BUCKET + APPLICATION + '/pipelines/pipeline' + str(i) + '.kt', overwrite=True)
    )

    kt_check = hc.read_table(BUCKET + APPLICATION + '/pipelines/pipeline' + str(i) + '.kt')
    n_samples = kt_check.count()

    print('')
    print('Pipeline: ', i)
    print('Fields: ', fields)
    print('nSamples: ', '{:,}'.format(n_samples))
    pprint(kt_check.schema)
